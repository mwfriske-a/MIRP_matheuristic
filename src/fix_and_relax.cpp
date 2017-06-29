#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <ilcplex/ilocplex.h>
#include "../include/util.h"
#include "../include/model.h"
#define NDEBUG
#include <assert.h>
#define UselessVariables
//~ #define NSpotMarket
//~ #define NTravelAtCapacity
//~ #define NTravelEmpty
//~ #define NBerthLimit
#define NBranching
#define NValidInequalities
#define NOperateAndGo
#define NOperateAndGo2
#define NBetas
#define NNoRevisits	
#define N2PortNorevisit
//~ #define NFixZvar
#define NRollingHorizon

ILOSTLBEGIN

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloNumArray> NumNumMatrix;

using namespace std;
using mirp::Model;
using mirp::Instance;

#define zero -0.0001

/* Build a model for the first iteration of fix and relax
 * All variables are created (including that are in end block)
 * Just 1st interval have integer variables x and z
 * Variables that are in the end block are not included in any constraint and objective function
 * Constraints of arcs/nodes in the and block are created, but are empty 
 */
void Model::buildFixAndRelaxModel(IloEnv& env, Instance inst, const double& nIntervals, const int& endBlock){
	int j,t,v,a;
	int timePerInterval = inst.t/nIntervals;  //Also starting index of relaxed block
	int J = inst.numTotalPorts;
	int T = inst.t;
	double intPart;
	int tOEB = inst.t - (timePerInterval*max(0.0,endBlock-modf(nIntervals, &intPart))) ; //Time periods out of End Block (index tOEB is the first in the end block(
	int V = inst.speed.getSize(); //# of vessels
	unordered_map<string,double>::const_iterator it;
	
	//Init Variables, converters and storage of variable values
	#ifndef NSpotMarket
	alpha = NumVarMatrix(env, J);	
	#endif
	#ifndef NBetas
	beta = NumVarMatrix(env, J);
	#endif
	sP = NumVarMatrix(env,J);
	f = IloArray<NumVarMatrix>(env, V);
	sV = NumVarMatrix(env,V);	
	x = IloArray<IloBoolVarArray>(env, V);
	z = IloArray<IloArray<IloBoolVarArray> >(env,V);
	
	//Additional variables for branching
	#ifndef NBranching	
	y = IloIntVarArray(env, J, 0, IloInfinity);	
	#endif
	
	//For implementation OperateAndGo2
	#ifndef NOperateAndGo2
	w = IloArray<IloBoolVarArray>(env,V);
	//~ convertW = IloArray<IloArray<IloConversion> >(env, V);	
	wValue = IloArray<IloNumArray>(env, V);
	#endif
	
	convertX = IloArray<IloArray<IloConversion> >(env, V);
	convertZ = IloArray<IloArray<IloArray<IloConversion> > >(env, V);	
	
	xValue = IloArray<IloNumArray>(env, V);
	zValue = IloArray<IloArray<IloNumArray> >(env, V);
	
	for(v=0;v<V;v++){
		f[v] = NumVarMatrix(env,J);
		sV[v] = IloNumVarArray(env, T+1);
		z[v] = IloArray<IloBoolVarArray>(env,J);		
		x[v] = IloBoolVarArray(env, inst.arcs[v].size());
		zValue[v] = IloArray<IloNumArray>(env, J);
		xValue[v] = IloNumArray(env, inst.arcs[v].size());
		convertZ[v] = IloArray<IloArray<IloConversion>> (env, J);	
		convertX[v] = IloArray<IloConversion>(env, inst.arcs[v].size());
		
		#ifndef NOperateAndGo2
		w[v] = IloBoolVarArray(env,T);
		//~ convertW[v] = IloArray<IloConversion> (env, T);	
		wValue[v] = IloNumArray(env, T);
		#endif
		
		for(a=0; a<x[v].getSize(); a++){
			stringstream ss;			
			ss << "x_(" << inst.arcs[v][a] << ")," << v;
			x[v][a].setName(ss.str().c_str());
			#ifdef NRollingHorizon
			if(inst.shouldRelaxArc(timePerInterval, v, a)){ //Kepp integrality only in the first interval
				convertX[v][a] = IloConversion(env,x[v][a], ILOFLOAT); 
				model.add(convertX[v][a]);
			}
			#endif
		}
		for(t=0;t<T+1;t++){
			stringstream ss, ss1;			
			ss << "supplyOnVessel_" << v << "," << t;
			
			if (t > inst.firstTimeAv[v])
				sV[v][t] = IloNumVar(env, 0, inst.q_v[v], ss.str().c_str());
			else
				sV[v][t] = IloNumVar(env, inst.s_v0[v], inst.s_v0[v], ss.str().c_str());
			#ifndef NOperateAndGo2			
			if(t<T){
				//~ convertW[v][t] = IloConversion (env, w[v][t], ILOFLOAT);	//Relaxing w variables
				ss1 << "w_"<<v<<","<<t;
				w[v][t].setName(ss1.str().c_str());				
				//~ model.add(convertW[v][t]);
			}
			#endif			
		}		
		
		for(j=0;j<J;j++){
			f[v][j] = IloNumVarArray(env, T);
			z[v][j] = IloBoolVarArray(env, T);
			zValue[v][j] = IloNumArray(env, T);
			convertZ[v][j] = IloArray<IloConversion> (env, T);				
			for(t=0;t<T;t++){
				#ifdef NRollingHorizon
				if (t >= timePerInterval){
					convertZ[v][j][t] = IloConversion (env, z[v][j][t], ILOFLOAT);	//Relaxing z variables
					model.add(convertZ[v][j][t]);										
				}
				#endif
				stringstream ss;
				ss << "f_(" << j << "," << t << ")," << v;
				int maxAmountOperation = min(inst.q_v[v], inst.f_max_jt[j][t]);
				f[v][j][t] = IloNumVar(env, 0, maxAmountOperation, ss.str().c_str());				
				
				ss.str(string());
				ss << "z_(" << j << "," << t << ")," << v;
				z[v][j][t].setName(ss.str().c_str());	
				#ifndef UselessVariables
				if(inst.inArcs[v][j][t].getSize() < 1){//Only if there is no entering arc in the node
					f[v][j][t].setBounds(0,0);
					z[v][j][t].setBounds(0,0);
				}
				#endif							
			}			
		}
	}
	
	for(j=0;j<J;j++){
		#ifndef NSpotMarket
		alpha[j] = IloNumVarArray(env,0,inst.alp_max_jt[j]);
		#endif
		#ifndef NBetas
		beta[j] = IloNumVarArray(env,T,0,IloInfinity);
		#endif
		
		sP[j] = IloNumVarArray(env, inst.sMin_jt[j], inst.sMax_jt[j]);
		sP[j].add(IloNumVar(env, inst.sMin_jt[j][0], inst.sMax_jt[j][0])); //Add more one variable to represent inventory at source node (note: using index 0 instead index T does not interfer in the current tested instances)
		for(t=0;t<T+1;t++){
			stringstream ss,ss1;
			if(t<T){
				#ifndef NSpotMarket			
				ss << "alpha_(" << j << "," << t << ")";				
				alpha[j][t].setName(ss.str().c_str());							
				#endif
				#ifndef NBetas
				ss1 << "beta_(" << j << "," << t << ")";
				beta[j][t].setName(ss1.str().c_str());
				#endif
			}
			ss.str(string());
			ss << "sP_(" << j << "," << t << ")";			
			sP[j][t].setName(ss.str().c_str());
			//~ sP[j][t].setBounds(-IloInfinity,IloInfinity); //There is no infeasibility if port inventory is unlimited
			
		}				
		//Initial inventory at each port
		sP[j][0].setBounds(inst.s_j0[j],inst.s_j0[j]); 
		stringstream ss;
		ss << "sP_(" << j << "," << 0 << ")";
		sP[j][0].setName(ss.str().c_str());
	}
	
	///Objective function - Just add variables z and x that are in out of end block (EOB)
	IloExpr expr(env);
	IloExpr expr1(env);
	for(v=0;v<V;v++){
		for(j=0;j<J;j++){
			for(t=0;t<tOEB;t++){
				#ifndef UselessVariables				
				if (inst.inArcs[v][j][t].getSize() > 0){ //Only if there is a intering arc in the node
				#endif
					expr += inst.r_jt[j][t]*f[v][j][t];				
					expr1 += (t*inst.perPeriodRewardForFinishingEarly)*z[v][j][t];					
				#ifndef UselessVariables
				}
				#endif
			}		
		}		
		for (a=0;a<inst.arcs[v].size();a++){
			//Add only the arc in the model that are not in the end block
			if (inst.isInModel(tOEB, v, a)){													
				it = inst.c_va[v].find(inst.arcs[v][a]);				
				if (it == inst.c_va[v].end()){
					cout << "Not found key for arc " << inst.arcs[v][a] << endl;
					exit(1);
				}				
				expr1 +=  it->second * x[v][a];						
			}			
		}
	}	
	for(j=0;j<J;j++){
		for(t=0;t<tOEB;t++){			
			#ifndef NSpotMarket
			expr1 += inst.p_jt[j][t]*alpha[j][t];
			#endif
			#ifndef NBetas
			expr1 += 1000*beta[j][t];
			#endif
		}		
	}
	//~ obj.setSense(IloObjective::Maximize);
	obj.setExpr(expr1-expr);	
	
	model.add(obj);
	
	///Constraints
	//Flow balance in port-time nodes, source and sink
	flowBalance = IloArray<IloArray<IloRangeArray> >(env,V);
	for(v=0;v<V;v++){
		flowBalance[v] = IloArray<IloRangeArray>(env,J+1);
		for(j=0;j<J+1;j++){ //J+1 used for access the information of source and sink node
			if(j<J){ //If is a real port
				flowBalance[v][j] = IloRangeArray(env,T, 0, 0);
				for(t=0;t<T;t++){	
					#ifndef UselessVariables
					if (inst.inArcs[v][j][t].getSize() > 0){ //Only if there is a intering arc in the node
					#endif				
						expr.clear();
						//~ cout << "Out arcs v["<<v<<"]["<<j<<"]["<<t<<"] = ";
						for (a=0;a<inst.outArcs[v][j][t].getSize();a++){						
							int idA = inst.outArcs[v][j][t][a];
							//Verify if arc and node are in the model
							if ( inst.isInModel(tOEB, v, idA) && t < tOEB){
								expr += x[v][idA];						
								//~ cout << inst.arcs[v][idA] << " " ;
							}
						}
						//~ cout << endl;
						for (a=0;a<inst.inArcs[v][j][t].getSize();a++){
							int idA = inst.inArcs[v][j][t][a];						
							//Verify if arc and node is in model
							if ( inst.isInModel(tOEB, v, idA) && t < tOEB)
								expr += -x[v][idA];		
						}
						flowBalance[v][j][t].setExpr(expr);						
					#ifndef UselessVariables
					}
					#endif							
					stringstream ss;
					ss << "flowBalance_(" << j << "," << t << ")," << v;					
					flowBalance[v][j][t].setName(ss.str().c_str());	
														
				}
				model.add(flowBalance[v][j]);
			}else{
				flowBalance[v][j] = IloRangeArray(env,2, 1, 1); // Two index, 0 for the source and 1 for sink node
				expr.clear();
				expr1.clear();
				stringstream ss;
				stringstream ss1;
				ss << "flowBalanceSource" << "," << v;
				ss1 << "flowBalanceSink" << "," << v;
				for (a=0;a<inst.outArcs[v][J][0].getSize();a++){
					//summing for source node
					int idA = inst.outArcs[v][J][0][a];					
					//~ if ( inst.isInModel(tOEB, v, idA) && t < tOEB)//Verify if arc and node is in model
						expr += x[v][idA];										
				}
				for (a=0;a<inst.inArcs[v][J][0].getSize();a++){
					//summing for sink node
					int idA = inst.inArcs[v][J][0][a];					
					if (inst.isInModel(tOEB, v, idA))
						expr1 += x[v][idA];		
				}
				//Flow Balance source node
				flowBalance[v][j][0].setExpr(expr);
				flowBalance[v][j][0].setName(ss.str().c_str());					
				//Flow Balance sink node
				flowBalance[v][j][1].setExpr(expr1);
				flowBalance[v][j][1].setName(ss1.str().c_str());
				
				//~ model.add(flowBalance[v][j]); //Both source and sink node flow				
				model.add(flowBalance[v][j][0]); //Only source node flow TODO Evaluate if it is necessary.
			}
		}
	}
	//Inventory balance ports 
	portInvBalance = IloArray<IloRangeArray>(env,J);
	for(j=0;j<J;j++){
		portInvBalance[j] = IloRangeArray(env,T);
		for(t=0;t<T;t++){ 
			expr.clear();
			if (t < tOEB){
				for(v=0;v<V;v++){
					#ifndef UselessVariables
					if (inst.inArcs[v][j][t].getSize() > 0) //Only if there is a intering arc in the node
					#endif
					expr += -f[v][j][t];
				}
				#ifndef NSpotMarket			
				expr += -alpha[j][t];
				#endif
				#ifndef NBetas
				expr += -beta[j][t];
				#endif
				stringstream ss;
				ss << "invBalancePort_(" << j << "," << t << ")";				
				portInvBalance[j][t] = IloRange(env, inst.delta[j]*inst.d_jt[j][t], 
				sP[j][t+1] - sP[j][t] - inst.delta[j] * expr,
				inst.delta[j]*inst.d_jt[j][t], 
				ss.str().c_str());
			}else{ //The constraint will no be considered (it is in end block, and is empty)
				stringstream ss;
				ss << "invBalancePort_(" << j << "," << t << ")";				
				portInvBalance[j][t] = IloRange(env, 0, expr, 0, ss.str().c_str());
			}
		}
		//Fix Inventory	in the first time period		
		//~ stringstream ss;
		//~ ss << "invBalancePort_(" << j << "," << 0 << ")";
		//~ portInvBalance[j][0] = IloRange(env, inst.s_j0[j], sP[j][0], inst.s_j0[j], ss.str().c_str());	
		model.add(portInvBalance[j]);
	}
	
	//Vessels inventory balance 
	vesselInvBalance = IloArray<IloRangeArray>(env,V);
	for(v=0;v<V;v++){
		vesselInvBalance[v] = IloRangeArray(env,T);			
		for(t=0;t<T;t++){
			expr.clear();
			if (t<tOEB){
				expr += sV[v][t];
				for(j=0;j<J;j++){
					#ifndef UselessVariables
					if (inst.inArcs[v][j][t].getSize() > 0) //Only if there is a intering arc in the node
					#endif
					expr += inst.delta[j]*f[v][j][t];
				}
			}
			stringstream ss;
			ss << "invBalanceVessel_" << v << "," << t;
			if (t < tOEB)
				vesselInvBalance[v][t] = IloRange(env, 0, -sV[v][t+1] + expr, 0, ss.str().c_str());			
			else
				vesselInvBalance[v][t] = IloRange(env, 0, expr, 0, ss.str().c_str());			
		}		
		//~ stringstream ss;
		//~ ss << "invBalanceVessel_" << v << "," << 0;
		//~ vesselInvBalance[v][0] = IloRange(env, inst.s_v0[v], sV[v][0], inst.s_v0[v], ss.str().c_str());
		
		model.add(vesselInvBalance[v]);
	}
	
	//Berth Limit
	#ifndef NBerthLimit
	berthLimit = IloArray<IloRangeArray>(env, J);
	for(j=0;j<J;j++){
		berthLimit[j] = IloRangeArray(env,T);
		for(t=0;t<T;t++){
			expr.clear();
			stringstream ss;
			ss << "berthLimit_(" << j << "," << t << ")";
			if(t < tOEB){
				for(v=0;v<V;v++){
					#ifndef UselessVariables
					if (inst.inArcs[v][j][t].getSize() > 0) //Only if there is a intering arc in the node
					#endif
					expr += z[v][j][t];				
				}			
				berthLimit[j][t] = IloRange(env,-IloInfinity, expr, inst.b_j[j], ss.str().c_str());						
			}else
				berthLimit[j][t] = IloRange(env,0, expr, 0, ss.str().c_str());						
		}
		model.add(berthLimit[j]);		
	}
	#endif
	
	//Vessels only attempt to load/discharge at node only if the vessel is actually at the node
	atemptToOperate = IloArray<IloArray<IloRangeArray> > (env, J);
	for(j=0;j<J;j++){
		atemptToOperate[j] = IloArray<IloRangeArray>(env, T);
		for(t=0;t<T;t++){
			atemptToOperate[j][t] = IloRangeArray(env, V);
			for(v=0;v<V;v++){
				expr.clear();
				stringstream ss;
				ss << "zOnlyIfSumx_(" << j << "," << t << ")," << v;
				#ifndef UselessVariables
				if (inst.inArcs[v][j][t].getSize() > 0){ //Only if there is a intering arc in the node
				#endif
					if (t < tOEB){					
						for(a=0;a<inst.inArcs[v][j][t].getSize();a++){
							int idA = inst.inArcs[v][j][t][a];
							expr += x[v][idA];
						}					
						atemptToOperate[j][t][v] = IloRange(env, -IloInfinity, z[v][j][t] - expr, 0, ss.str().c_str());				
					}else //Just descondirer the restriction for the port time and vessel, as there is no entering arc ((j,t),(i,t+x)) from a t+x < t
						atemptToOperate[j][t][v] = IloRange(env, 0, 0, ss.str().c_str());				
				#ifndef UselessVariables
				}else{
					expr.clear();
					atemptToOperate[j][t][v] = IloRange(env, 0, z[v][j][t], 0, ss.str().c_str());				
				}
				#endif
			}
			model.add(atemptToOperate[j][t]);			
		}
	}
	
	//Vessels must travel at capacity and empty 
	travelAtCapacity = IloArray<IloRangeArray> (env, V);
	travelEmpty = IloArray<IloRangeArray> (env, V);
	for(v=0;v<V;v++){
		travelAtCapacity[v] = IloRangeArray(env,inst.arcs[v].size()); ///Altough just travel and sink arcs are incluede in this constraint set, the size is on the number of arcs
		travelEmpty[v] = IloRangeArray(env, inst.arcs[v].size());
		for(a=0;a<inst.arcs[v].size();a++){			
			int timeJ1, timeJ2;	
			stringstream ss, ss1;
			ss << "travelAtCap_(" << inst.arcs[v][a] << ")," << v;
			ss1 << "travelEmpty_(" << inst.arcs[v][a] << ")," << v;			
			int arcType = inst.travelArcType(inst.arcs[v][a], timeJ1, timeJ2);
			if ( inst.isInModel(tOEB, v, a) ){ 
				travelAtCapacity[v][a] = IloRange(env, -IloInfinity, 0, ss.str().c_str());		
				travelEmpty[v][a] = IloRange(env, -IloInfinity, inst.q_v[v], ss1.str().c_str());		
				switch (arcType)
				{
					case 1:
						#ifndef NTravelAtCapacity
						 //timeJ1 is incresead in 1 because in vector sV time 0 is considered as source node										
						travelAtCapacity[v][a].setExpr(-sV[v][timeJ1+1] + inst.q_v[v]*x[v][a]);	
						//~ cout << "Travel at capacity " << v << " arc " << inst.arcs[v][a] << endl;														
						#endif
						break;
					case 2:
						#ifndef NTravelEmpty
						travelEmpty[v][a].setExpr(inst.q_v[v]*x[v][a] + sV[v][timeJ1+1]);						
						#endif					
						break;
					case 3:
						#ifndef NTravelAtCapacity
						travelAtCapacity[v][a].setExpr(-sV[v][timeJ1+1] + inst.q_v[v]*x[v][a]);																	
						#endif					
						break;
					case 4:
						#ifndef NTravelEmpty
						travelEmpty[v][a].setExpr(inst.q_v[v]*x[v][a] + sV[v][timeJ1+1]);						
						#endif
						break;
					default:
						break;
				}
			}else{
				travelAtCapacity[v][a] = IloRange(env, 0, expr, 0, ss.str().c_str());				
				travelEmpty[v][a] = IloRange(env, 0, expr, 0, ss1.str().c_str());				
			}
		}		
		model.add(travelAtCapacity[v]);
		model.add(travelEmpty[v]);
	}
	
	//Cumulative amount of product that can be purchased from or sold to a spot market by each port
	#ifndef NSpotMarket
	cumSlack = IloRangeArray(env, J);
	for(j=0;j<J;j++){
		stringstream ss;
		ss << "cumSlack_" << j;
		expr.clear();
		for(t=0;t<T;t++){
			if (t < tOEB)
				expr += alpha[j][t];							
		}
		cumSlack[j] = IloRange(env, expr, inst.alp_max_j[j], ss.str().c_str());
	}
	model.add(cumSlack);
	#endif
	
	//Amount loaded/discharged must be in the pre-specified interval 
	operationLowerLimit = IloArray<IloArray<IloRangeArray> > (env, V);
	operationUpperLimit = IloArray<IloArray<IloRangeArray> > (env, V);
	for(v=0;v<V;v++){
		operationLowerLimit[v] = IloArray<IloRangeArray>(env, J);
		operationUpperLimit[v] = IloArray<IloRangeArray>(env, J);
		for(j=0;j<J;j++){
			operationLowerLimit[v][j] = IloRangeArray(env, T);
			operationUpperLimit[v][j] = IloRangeArray(env, T);
			for(t=0;t<T;t++){
				expr.clear();
				stringstream ss, ss1;
				ss << "fjmin_(" << j << "," << t << ")," << v;
				ss1 << "fjmax_(" << j << "," << t << ")," << v;
				#ifndef UselessVariables
				if (inst.inArcs[v][j][t].getSize() > 0){ //Only if there is a entering arc in the node
				#endif
					if (t < tOEB){
						IloNum minfMax = min(inst.f_max_jt[j][t], inst.q_v[v]);
						operationLowerLimit[v][j][t] = IloRange(env, 0, f[v][j][t] - inst.f_min_jt[j][t]*z[v][j][t] , IloInfinity, ss.str().c_str());				
						operationUpperLimit[v][j][t] = IloRange(env, -IloInfinity, f[v][j][t] - minfMax*z[v][j][t],	0, ss1.str().c_str());				
					}else{
						operationLowerLimit[v][j][t] = IloRange(env, 0, expr, 0, ss.str().c_str());				
						operationUpperLimit[v][j][t] = IloRange(env, 0, expr, 0, ss1.str().c_str());			
					}	
				#ifndef UselessVariables
				}else{
					operationLowerLimit[v][j][t] = IloRange(env, 0, f[v][j][t], 0, ss.str().c_str());				
					operationUpperLimit[v][j][t] = IloRange(env, 0, 0, ss1.str().c_str());							
				}
				#endif
				
			}
			model.add(operationLowerLimit[v][j]);
			model.add(operationUpperLimit[v][j]);
		}		
	}
	
	//Additional constraints definition - y_j = number of visits at each port (entering arcs - no including waiting)
	#ifndef NBranching
	y_Sum = IloRangeArray(env, J,0,0);
	for(j=0;j<J;j++){
		expr.clear();
		for(v=0;v<V;v++){			
			for(t=0;t<T;t++){
				for (a=0;a<inst.inArcs[v][j][t].getSize();a++){
					int idA = inst.inArcs[v][j][t][a];						
					int j1,j2,t1,t2,type;
					type = inst.getArcType(inst.arcs[v][idA], t1, t2, j1, j2);
					if(type == 1 || type == 2 || (type ==-1 && j1 != j2)){
						//Verify if arc is in model
						if ( inst.isInModel(tOEB, v, idA))
							expr += x[v][idA];		
					}
				}			
			}			
		}
		stringstream ss;
		ss << "y_sum_" << j;
		y_Sum[j] = IloRange(env, 0, expr-y[j], 0, ss.str().c_str());		
		model.add(y_Sum[j]);		
	}	
	#endif
	#ifndef NValidInequalities
	//Valid inequalities 
	//Minimum number of visits
	minVisits = IloArray<IloRangeArray>(env, J);
	for(j=0;j<J;j++){
		minVisits[j] = IloRangeArray(env, T, -IloInfinity, 0);
		for(t=0;t<T;t++){
			stringstream ss;
			ss << "minVisit_" << j << ","<<t;
			minVisits[j][t].setName(ss.str().c_str());			
			if(t<tOEB){ // If is in the model
				expr.clear();
				double sum_d_jt=0, lb_RHS=0, denominator=0, max_dj=0;				
				for(int t1=0;t1<=t;t1++){
					sum_d_jt += inst.d_jt[j][t1];
					double dj = inst.sMax_jt[j][t1] + inst.d_jt[j][t1] - inst.sMin_jt[j][t1];
					if (max_dj < dj)
						max_dj = dj;
					for(v=0;v<V;v++)
						expr += -z[v][j][t1];					
				}
				denominator = min(inst.f_max_jt[j][t], max_dj);
				
				//Define the numerator according to the type of port
				if (inst.typePort[j] == 0) //Loading port
					sum_d_jt += inst.s_j0[j] - inst.sMax_jt[j][t];
				else //Discharging port
					sum_d_jt += inst.sMin_jt[j][t] - inst.s_j0[j];				
				lb_RHS = ceil(sum_d_jt/denominator);				
				//~ cout << j << "," << t << "=" << lb_RHS << endl;
				minVisits[j][t].setUB(-lb_RHS);
				minVisits[j][t].setExpr(expr);
				
				//~ cout << ss.str() << " = " << lb_RHS << endl;
			}
		}
		model.add(minVisits[j]);
	}
	#endif
	#ifndef NOperateAndGo
	//Forces a vessel go to another region after operates if a vessel can load(unload) at a port just in 1 time period
	operateAndGo = IloArray<IloArray<IloRangeArray> >(env, V);
	for (v=0;v<V;v++){
		operateAndGo[v] = IloArray<IloRangeArray>  (env, J);
		for(j=0;j<J;j++){						
			operateAndGo[v][j] = IloRangeArray(env,T);
			for(t=0;t<T;t++){		
				stringstream ss;
				ss << "operateAndGo_"<<v<<"_"<<j<<","<<t;
				expr.clear();		
				if (t<tOEB){  //If is in model										
					if (inst.q_v[v] < inst.f_max_jt[j][t]){ //Only if a vessel can load(unload) in a port in just 1 time period						
						for(a=0;a<inst.outRegionArcs[v][j][t].getSize();a++){
							int idA = inst.outRegionArcs[v][j][t][a];
							//Verify if arc and node are in the model
							if (inst.isInModel(tOEB, v, idA))
								expr += x[v][idA];							
						}
						expr += -z[v][j][t];						
					}
				}
				operateAndGo[v][j][t] = IloRange(env, 0, expr, 0,ss.str().c_str());
			}
			model.add(operateAndGo[v][j]);
		}
	}
	#endif
	
	#ifndef NOperateAndGo2
	operateAndGo1 = IloArray<IloArray<IloRangeArray> >(env,V);	
	sideIf = IloArray<IloArray<IloArray<IloIfThen> > >(env, V);
	for (v=0;v<V;v++){
		operateAndGo1[v] = IloArray<IloRangeArray>(env, J);		
		//Side constraints			
		sideIf[v] = IloArray<IloArray<IloIfThen> >(env, J);		
		for(j=0;j<J;j++){
			sideIf[v][j] = IloArray<IloIfThen>(env,T);			
			operateAndGo1[v][j] = IloRangeArray(env,T, -IloInfinity, 0);
			for(t=0;t<T;t++){
				stringstream ss,ss1, ss2;
				ss << "operateAndGo1_"<<v<<"_"<<j<<","<<t;				
				ss2 << "sideOpG1_"<<v<<"_"<<j<<","<<t;
				expr.clear();
				if (inst.typePort[j] == 0)
					sideIf[v][j][t] = IloIfThen(env, (z[v][j][t] == 1) && (inst.q_v[v] - sV[v][t+1] <= inst.f_min_r[inst.idRegion[j]]), w[v][t] == 1, ss2.str().c_str());
					//~ sideIf[v][j][t] = IloIfThen(env, (z[v][j][t] == 1) && (inst.q_v[v] - sV[v][t+1] <= 0.001), w[v][t] == 1, ss2.str().c_str());
				else
					sideIf[v][j][t] = IloIfThen(env, (z[v][j][t] == 1) && (sV[v][t+1] <= inst.f_min_r[inst.idRegion[j]]), w[v][t] == 1, ss2.str().c_str());
					//~ sideIf[v][j][t] = IloIfThen(env, (z[v][j][t] == 1) && (sV[v][t+1] <= 0.001), w[v][t] == 1, ss2.str().c_str());
				if(t<tOEB){
					for(a=0;a<inst.outRegionArcs[v][j][t].getSize();a++){
						int idA = inst.outRegionArcs[v][j][t][a];
						//Verify if arc and node are in the model
						if (inst.isInModel(tOEB, v, idA))
							expr += -x[v][idA];							
					}				
					operateAndGo1[v][j][t].setExpr(expr + w[v][t]);					
					model.add(sideIf[v][j][t]);
				}else{
					operateAndGo1[v][j][t].setBounds(0,0);					
				}					
				operateAndGo1[v][j][t].setName(ss.str().c_str());				
			}
			model.add(operateAndGo1[v][j]);				
		}				
	}
	#endif	
	
	#ifndef NNoRevisits
	noRevisit = IloArray<IloArray<IloRangeArray> >(env,V);
	for(v=0;v<V;v++){
		noRevisit[v] = IloArray<IloRangeArray>(env,J);
		for(j=0;j<J;j++){ 			
			noRevisit[v][j] = IloRangeArray(env,T, 0, 1);
			for(t=0;t<T;t++){				
				expr.clear();										
				if(t<tOEB){ //Only for nodes in model
					for (a=0;a<inst.outArcs[v][j][t].getSize();a++){
						int idA = inst.outArcs[v][j][t][a];
						int j1,j2,t1,t2,type;
						type = inst.getArcType(inst.arcs[v][idA], t1, t2, j1, j2);
						if(type == -1 && j1 != j2 && inst.idRegion[j1] == inst.idRegion[j2]){ //If is intra-regional arc
							//Verify if arc is in model
							if (inst.isInModel(tOEB, v, idA)){
								expr += x[v][idA];
							}
						}
					}
					int forbidTimeArrive = t+2*inst.maxTimeIntraReg[v][j];					
					for(int ta=t; ta<forbidTimeArrive;ta++){ //For each next nodes of port j 
						for (a=0;a<inst.inArcs[v][j][ta].getSize();a++){						
							int idA = inst.outArcs[v][j][ta][a];
							int j1,j2,t1,t2,type;
							type = inst.getArcType(inst.arcs[v][idA], t1, t2, j1, j2);
							if(type == -1 && j1 != j2 && inst.idRegion[j1] == inst.idRegion[j2]){ //If is intra-regional arc
								//Verify if arc is in the model
								if (inst.isInModel(tOEB, v, idA)){
									expr += x[v][idA];						
								}
							}						
						}
					}				
				}else
					noRevisit[v][j][t].setBounds(0,0);
				noRevisit[v][j][t].setExpr(expr);												
				stringstream ss;
				ss << "noRevisit_(" << j << "," << t << ")," << v;					
				noRevisit[v][j][t].setName(ss.str().c_str());				
			}
			model.add(noRevisit[v][j]);
		}
	}
	#endif
	
	#ifndef N2PortNorevisit
	twoPortsNoRevisit = IloArray<IloRangeArray>(env, inst.loadReg + inst.discReg);	
	//First load regions 
	for(int tpR=0;tpR<=1;tpR++){
		for(int reg=0;reg<inst.identifyPort[tpR].getSize();reg++){
			int r; //Unique identifier for each region (first load, after discharge)
			if(tpR==0)
				r = reg;
			else
				r = reg+inst.loadReg;			
			twoPortsNoRevisit[r] = IloRangeArray(env, V, -IloInfinity, 0);						
			for(v=0;v<V;v++){
				expr.clear();
				for(int port=0;port < inst.identifyPort[tpR][reg].getSize();port++){ //For each port of region reg of type tpR										
					j = inst.identifyPort[tpR][reg][port];					
					for(t=0;t<tOEB;t++){ //Iterating until end block
						for (a=0;a<inst.outArcs[v][j][t].getSize();a++){							
							int idA = inst.outArcs[v][j][t][a];
							int j1,j2,t1,t2,type;
							type = inst.getArcType(inst.arcs[v][idA], t1, t2, j1, j2);
							if( (type == -1) && (inst.idRegion[j1] == inst.idRegion[j2]) && (j1 != j2)){ //Intra-regional arcs (no considering waiting arcs);
								if (inst.isInModel(tOEB, v, idA)){
									expr += x[v][idA];						
								}
							}
						}						
						for (a=0;a<inst.inRegionArcs[v][j][t].getSize();a++){ //Arcs arriving at region through node v_j,t
							int idA = inst.inRegionArcs[v][j][t][a];
							if (inst.isInModel(tOEB, v, idA)){
								expr += -x[v][idA];						
							}						
						}						
					}					
				}
				stringstream ss;
				ss << "twoPorts_noRevisit_" << v << "("<<r<<")";				
				twoPortsNoRevisit[r][v].setExpr(expr);
				twoPortsNoRevisit[r][v].setName(ss.str().c_str());				
			}
			model.add(twoPortsNoRevisit[r]);
		}
	}	
	#endif
}
/* Param p = 0 If algorithm can extract var solution values; 1 otherwise*/
void Model::fixSolution(IloEnv& env, Instance inst, const int& t3S, const int& t3F,const int& p){	
	if(p == 0){
		//Get the values 
		//~ getSolVals(env, inst);		//Get all values
		getSolValsW(env, inst, t3S, t3F);	//Get only values of the interval
	}
	//Fix solution - If p == 1, consider the values from the previous p==0
	for(int v=0; v<inst.speed.getSize(); v++){
		#ifndef NFixZvar
		for(int j=0; j<inst.numTotalPorts; j++){							
			for(int t=t3S; t<t3F; t++){										
				z[v][j][t].setBounds(round(zValue[v][j][t]), round(zValue[v][j][t]));					
			}
		}
		#endif
		//Fixing X variables
		for(int a=0; a<x[v].getSize(); a++){
			int t1,t2, arcType;		
			arcType = inst.travelArcType(inst.arcs[v][a],t1,t2);								
			
			if( ( (arcType != 3) && (arcType != 4) && // Never fix a sink ark
			//~ if(  ((arcType == 3 || arcType == 4) && (t1 >= t3S && t1 < t3F) ) || // Is a sink arc in the interval
			(t2 >= t3S && t2 < t3F) ) ){  // Arc times between fixing times	(and travessing fixed -> new fixed)
				x[v][a].setBounds(round(xValue[v][a]), round(xValue[v][a]));
		   }			   
		}
		#ifndef NOperateAndGo2
		//~ for(int t=t3S; t<t3F; t++){
			//~ w[v][t].setBounds(round(wValue[v][t]), round(wValue[v][t]));
		//~ }
		#endif
	}
}

void Model::decreaseEndBlock (IloEnv& env, Instance inst, const double& nIntervals, const int& t2S, const int& t2F){
	//t2S = tOEB in the first call of this method
	int i,j,t,v,a;
	int timePerInterval = inst.t/nIntervals;
	int J = inst.numTotalPorts;
	int T = inst.t;	
	int V = inst.speed.getSize(); //# of vessels
	unordered_map<string,double>::const_iterator it;
	
	///Objective function
	IloExpr current = cplex.getObjective().getExpr(); //Gets the current objective function
	IloExpr expr(env);
	IloExpr expr1(env);	
	for(v=0;v<V;v++){
		for(j=0;j<J;j++){
			for(t=t2S;t<t2F;t++){
				#ifndef UselessVariables
				if (inst.inArcs[v][j][t].getSize() > 0){ //Only if there is a intering arc in the node
				#endif
					expr += inst.r_jt[j][t]*f[v][j][t];				
					expr1 += (t*inst.perPeriodRewardForFinishingEarly)*z[v][j][t];					
				#ifndef UselessVariables
				}
				#endif
			}		
		}		
		for (a=0;a<inst.arcs[v].size();a++){
			//Add only the arc in the model that are not in the end block
			if (inst.shouldAddArc(t2S, t2F, v, a)){
				int t1,t2;
				it = inst.c_va[v].find(inst.arcs[v][a]);
				#ifndef NDEBUG
				int arcType = inst.travelArcType(inst.arcs[v][a], t1,t2);
				if (arcType == 0){
					cout << "Erro - source node was not in model in the first iteration" << endl;
					exit(1);				
				}				
				if (it == inst.c_va[v].end()){
					cout << "Not found key for arc " << inst.arcs[v][a] << endl;
					exit(1);
				}				
				#endif
				expr1 +=  it->second * x[v][a];						
			}			
		}
	}
	
	for(j=0;j<J;j++){
		for(t=t2S;t<t2F;t++){			
			#ifndef NSpotMarket
				expr1 += inst.p_jt[j][t]*alpha[j][t]; 
			#endif
			#ifndef NBetas
				expr1 += 1000*beta[j][t];		
			#endif			
		}		
	}
	obj.setExpr(current + expr1 - expr);	
	
	//Constraints
	//Flow balance in port-time nodes, source and sink
	for(v=0;v<V;v++){
		for(j=0;j<J+1;j++){ //J+1 used for access the information of source and sink node
			if(j<J){ //If is a real port				
				for(t=0;t<t2F;t++){ //May have arcs in constraints already in the model	(TODO Improve for not initalize t=0 ever)
					#ifndef UselessVariables
					if (inst.inArcs[v][j][t].getSize() > 0){ //Only if there is a intering arc in the node
					#endif
						//Get the current expr (empty or not) and updates it with new arcs 
						IloExpr current = flowBalance[v][j][t].getExpr();
						expr.clear();
						for (a=0;a<inst.outArcs[v][j][t].getSize();a++){						
							int idA = inst.outArcs[v][j][t][a];
							//Verify if the arc should be added 							
							if (inst.shouldAddArc(t2S, t2F, v, idA)) 
								expr += x[v][idA];						
						}
						for (a=0;a<inst.inArcs[v][j][t].getSize();a++){
							int idA = inst.inArcs[v][j][t][a];						
							//Verify if the arc should be added 
							if (inst.shouldAddArc(t2S, t2F, v, idA)) 
								expr += -x[v][idA];		
						}										
						flowBalance[v][j][t].setExpr(current + expr);
					#ifndef UselessVariables
					}
					#endif					
				}
				
			}else{ //Is source (0 - assuming it does not need to be updated) or sink (1) node								
				expr1.clear();				
				for (a=0;a<inst.inArcs[v][J][0].getSize();a++){
					//summing for sink nodes
					int idA = inst.inArcs[v][J][0][a];					
						if ( inst.shouldAddArc(t2S, t2F, v, idA))
						expr1 += x[v][idA];		
				}
				
				//Flow Balance sink node (need to be add?)
				expr.clear();
				expr = flowBalance[v][j][1].getExpr();
				flowBalance[v][j][1].setExpr(expr + expr1);								
				//~ model.add(flowBalance[v][j]); //Both source and sink node flow				
				//~ model.add(flowBalance[v][j][0]); //Only source node flow 
			}
		}
	}
	
	//Inventory balance ports
	for(j=0;j<J;j++){		
		for(t=t2S;t<t2F;t++){
			expr.clear();			
			for(v=0;v<V;v++){	
				#ifndef UselessVariables
				if (inst.inArcs[v][j][t].getSize() > 0) //Only if there is a intering arc in the node
				#endif			
					expr += -f[v][j][t];				
			}
			#ifndef NSpotMarket			
			expr += -alpha[j][t]; 
			#endif
			#ifndef NBetas
			expr +=  -beta[j][t];
			#endif
			stringstream ss;
			ss << "invBalancePort_(" << j << "," << t << ")";				
			portInvBalance[j][t].setBounds(inst.delta[j]*inst.d_jt[j][t], inst.delta[j]*inst.d_jt[j][t]);
			portInvBalance[j][t].setExpr(sP[j][t+1] - sP[j][t] - inst.delta[j] * expr);				
		}
		//~ model.add(portInvBalance[j]);
	}
	//Vessels inventory balance 	
	for(v=0;v<V;v++){		
		for(t=t2S;t<t2F;t++){
			expr.clear();	
			expr += sV[v][t];
			for(j=0;j<J;j++){
				#ifndef UselessVariables
				if (inst.inArcs[v][j][t].getSize() > 0) //Only if there is a intering arc in the node
				#endif
					expr += inst.delta[j]*f[v][j][t];
			}			
			vesselInvBalance[v][t].setExpr(-sV[v][t+1] + expr);			
		}				
		//~ model.add(vesselInvBalance[v]);
	}
	
	//Berth Limit
	#ifndef NBerthLimit	
	for(j=0;j<J;j++){		
		for(t=t2S;t<t2F;t++){
			expr.clear();			
			berthLimit[j][t].setBounds(-IloInfinity, inst.b_j[j]);
			for(v=0;v<V;v++){
				#ifndef UselessVariables
				if (inst.inArcs[v][j][t].getSize() > 0) //Only if there is a intering arc in the node
				#endif
					expr += z[v][j][t];				
			}			
			berthLimit[j][t].setExpr(expr);			
		}
		//~ model.add(berthLimit[j]);		
	}
	#endif
	
	//Vessels only attempt to load/discharge at node only if the vessel is actually at the node	
	for(j=0;j<J;j++){		
		for(t=t2S;t<t2F;t++){			
			for(v=0;v<V;v++){
				#ifndef UselessVariables
				if (inst.inArcs[v][j][t].getSize() > 0){ //Only if there is a intering arc in the node
				#endif
					expr.clear();									
					for(a=0;a<inst.inArcs[v][j][t].getSize();a++){
						int idA = inst.inArcs[v][j][t][a];
						expr += x[v][idA];
					}					
					atemptToOperate[j][t][v].setBounds(-IloInfinity, 0);
					atemptToOperate[j][t][v].setExpr(z[v][j][t] - expr);			
				#ifndef UselessVariables
				}
				#endif
			}
			//~ model.add(atemptToOperate[j][t]);			
		}
	}
	
	//Vessels must travel at capacity and empty	
	for(v=0;v<V;v++){				
		for(a=0;a<inst.arcs[v].size();a++){			
			int timeJ1, timeJ2;				
			int arcType = inst.travelArcType(inst.arcs[v][a], timeJ1, timeJ2);						
			switch (arcType)
			{
				case 1:
					#ifndef NTravelAtCapacity	
					//~ cout << "Travel Arc V " << v << ": " <<  inst.arcs[v][a] << " T1 " << timeJ1 << " T2 " << timeJ2 << endl;							 
					if (timeJ2 >= t2S && timeJ2 < t2F){ 
						travelAtCapacity[v][a].setBounds(-IloInfinity, 0);
						travelAtCapacity[v][a].setExpr(-sV[v][timeJ1+1] + inst.q_v[v]*x[v][a]);						
						//~ cout << "Travel at capacity " << v << " arc " << inst.arcs[v][a] << endl;			
					}
					#endif
					break;
				case 2:
					#ifndef NTravelEmpty
					if (timeJ2 >= t2S && timeJ2 < t2F){ 
						travelEmpty[v][a].setBounds(-IloInfinity, inst.q_v[v]);
						travelEmpty[v][a].setExpr(inst.q_v[v]*x[v][a] + sV[v][timeJ1+1]);						
						//~ cout << "Travel empty " << v << " arc " << inst.arcs[v][a] << endl;			
					}
					#endif					
					break;
				case 3:
					#ifndef NTravelAtCapacity
					if (timeJ1 >= t2S && timeJ1 < t2F){
						travelAtCapacity[v][a].setBounds(-IloInfinity, 0);
						travelAtCapacity[v][a].setExpr(-sV[v][timeJ1+1] + inst.q_v[v]*x[v][a]);
						//~ cout << "Travel at capacity " << v << " arc " << inst.arcs[v][a] << endl;			
					}
					#endif					
					break;
				case 4:
					#ifndef NTravelEmpty
					if (timeJ1 >= t2S && timeJ2 < t2F){ 
						travelEmpty[v][a].setBounds(-IloInfinity, inst.q_v[v]);
						travelEmpty[v][a].setExpr(inst.q_v[v]*x[v][a] + sV[v][timeJ1+1]);						
						//~ cout << "Travel empty " << v << " arc " << inst.arcs[v][a] << endl;			
					}
					#endif
					break;
				default:
					break;
			}
		}		
	}
	
	//Cumulative amount of product that can be purchased from or sold to a spot market by each port
	#ifndef NSpotMarket	
	for(j=0;j<J;j++){		
		expr.clear();
		expr = cumSlack[j].getExpr();
		for(t=t2S;t<t2F;t++){			
			expr += alpha[j][t];							
		}
		cumSlack[j].setExpr(expr);
	}
	model.add(cumSlack);
	#endif
	
	//Amount loaded/discharged must be in the pre-specified interval 	
	for(v=0;v<V;v++){		
		for(j=0;j<J;j++){			
			for(t=t2S;t<t2F;t++){
				#ifndef UselessVariables
				if (inst.inArcs[v][j][t].getSize() > 0){ //Only if there is a intering arc in the node
				#endif								
					operationLowerLimit[v][j][t].setBounds(0, IloInfinity);
					operationLowerLimit[v][j][t].setExpr(f[v][j][t] - inst.f_min_jt[j][t]*z[v][j][t]);								 	
					
					IloNum minfMax = min(inst.f_max_jt[j][t], inst.q_v[v]);
					operationUpperLimit[v][j][t].setBounds(-IloInfinity, 0);
					operationUpperLimit[v][j][t].setExpr( f[v][j][t] - minfMax*z[v][j][t]); 									
				#ifndef UselessVariables
				}
				#endif
			}
		}
	}
	
	//Additional constraint for branching	
	#ifndef NBranching	
	for(j=0;j<J;j++){
		expr.clear();
		expr = y_Sum[j].getExpr();
		for(v=0;v<V;v++){			
			for(t=0;t<t2F;t++){
				for (a=0;a<inst.inArcs[v][j][t].getSize();a++){
					int idA = inst.inArcs[v][j][t][a];						
					int j1,j2,t1,t2,type;
					type = inst.getArcType(inst.arcs[v][idA], t1, t2, j1, j2);
					if(type == 1 || type == 2 || (type ==-1 && j1 != j2)){						
						//Verify if arc is in model
						if (inst.shouldAddArc(t2S, t2F, v, idA)) 
							expr += x[v][idA];		
					}					
				}			
			}			
		}		
		y_Sum[j].setExpr(expr);				
	}	
	#endif
	
	#ifndef NValidInequalities
	//Valid inequalities 
	//Minimum number of visits	
	for(j=0;j<J;j++){		
		for(t=t2S;t<t2F;t++){				
			expr.clear();
			double sum_d_jt=0, lb_RHS=0, denominator=0, max_dj=0;				
			for(int t1=0;t1<=t;t1++){
				sum_d_jt += inst.d_jt[j][t1];
				double dj = inst.sMax_jt[j][t1] + inst.d_jt[j][t1] - inst.sMin_jt[j][t1];
				if (max_dj < dj)
					max_dj = dj;
				for(v=0;v<V;v++)
					expr += -z[v][j][t1];					
			}
			denominator = min(inst.f_max_jt[j][t], max_dj);
			
			//Define the numerator according to the type of port
			if (inst.typePort[j] == 0) //Loading port
				sum_d_jt += inst.s_j0[j] - inst.sMax_jt[j][t];
			else //Discharging port
				sum_d_jt += inst.sMin_jt[j][t] - inst.s_j0[j];				
			lb_RHS = ceil(sum_d_jt/denominator);				
			//~ cout << j << "," << t << "=" << lb_RHS << endl;
			minVisits[j][t].setUB(-lb_RHS);
			minVisits[j][t].setExpr(expr);			
		}		
	}
	#endif
	#ifndef NOperateAndGo
	//Forces a vessel go for another region after operates (and if have no enough capacity to operate)	
	for (v=0;v<V;v++){		
		for(j=0;j<J;j++){									
			for(t=0;t<t2F;t++){
				expr.clear();				
				if (inst.q_v[v] < inst.f_max_jt[j][t]){ //Only if a vessel can load(unload) in a port in just 1 time period						
					expr = operateAndGo[v][j][t].getExpr();
					for(a=0;a<inst.outRegionArcs[v][j][t].getSize();a++){
						int idA = inst.outRegionArcs[v][j][t][a];
						//Verify if arc is in the model
						if (inst.shouldAddArc(t2S, t2F, v, idA))
							expr += x[v][idA];
					}
					if(t>=t2S) //If expr for the constraint is empty - needs add variable z
						expr += -z[v][j][t];
				}
				operateAndGo[v][j][t].setExpr(expr);
			}		
		}
	}
	#endif	
	#ifndef NOperateAndGo2	
	for (v=0;v<V;v++){		
		for(j=0;j<J;j++){			
			for(t=0;t<t2F;t++){				
				expr.clear();				
				expr = operateAndGo1[v][j][t].getExpr();				
				for(a=0;a<inst.outRegionArcs[v][j][t].getSize();a++){
					int idA = inst.outRegionArcs[v][j][t][a];
					//Verify if arc is in the model
					if (inst.shouldAddArc(t2S, t2F, v, idA))
						expr += -x[v][idA];												
				}
				if(t>=t2S){//If expr for the constraint was empty in the previous iteration
					expr += w[v][t];
					//Side constraints		
					model.add(sideIf[v][j][t]);
					operateAndGo1[v][j][t].setBounds(-IloInfinity,0);					
				}					
				operateAndGo1[v][j][t].setExpr(expr);				
			}
		}				
	}
	#endif
	
	#ifndef NNoRevisits	
	for(v=0;v<V;v++){		
		for(j=0;j<J;j++){ 						
			for(t=0;t<t2F;t++){				
				expr.clear();
				expr = noRevisit[v][j][t].getExpr();
				if(t>=t2S) //When constraint was recently added
					noRevisit[v][j][t].setBounds(0,1);				
				for (a=0;a<inst.outArcs[v][j][t].getSize();a++){
					int idA = inst.outArcs[v][j][t][a];
					int j1,j2,t1,t2,type;
					type = inst.getArcType(inst.arcs[v][idA], t1, t2, j1, j2);
					if(type == -1 && j1 != j2 && inst.idRegion[j1] == inst.idRegion[j2]){ //If is intra-regional arc						
						if (inst.shouldAddArc(t2S,t2F,v,idA)){
							expr += x[v][idA];
						}
					}
				}
				int forbidTimeArrive = t+2*inst.maxTimeIntraReg[v][j];					
				for(int ta=t; ta<forbidTimeArrive;ta++){ //For each next nodes of port j 
					for (a=0;a<inst.inArcs[v][j][ta].getSize();a++){						
						int idA = inst.outArcs[v][j][ta][a];
						int j1,j2,t1,t2,type;
						type = inst.getArcType(inst.arcs[v][idA], t1, t2, j1, j2);
						if(type == -1 && j1 != j2 && inst.idRegion[j1] == inst.idRegion[j2]){ //If is intra-regional arc							
							if (inst.shouldAddArc(t2S,t2F,v,idA)){
								expr += x[v][idA];						
							}
						}						
					}
				}
				noRevisit[v][j][t].setExpr(expr);																
			}			
		}
	}
	#endif
	
	#ifndef N2PortNorevisit	
	//First load regions 
	for(int tpR=0;tpR<=1;tpR++){
		for(int reg=0;reg<inst.identifyPort[tpR].getSize();reg++){
			int r; //Unique identifier for each region (first load, after discharge)
			if(tpR==0)
				r = reg;
			else
				r = reg+inst.loadReg;			
			for(v=0;v<V;v++){
				expr.clear();
				expr = twoPortsNoRevisit[r][v].getExpr();
				for(int port=0;port<inst.identifyPort[tpR][reg].getSize();port++){ //For each port of region reg of type tpR
					j = inst.identifyPort[tpR][reg][port];
					for(t=0;t<t2F;t++){ 
						for (a=0;a<inst.outArcs[v][j][t].getSize();a++){
							int idA = inst.outArcs[v][j][t][a];
							int j1,j2,t1,t2,type;							
							type = inst.getArcType(inst.arcs[v][idA], t1, t2, j1, j2);
							if( (type == -1) && (inst.idRegion[j1] == inst.idRegion[j2]) && (j1 != j2)){ //Intra-regional arcs (no considering waiting arcs);
								if (inst.shouldAddArc(t2S,t2F,v, idA)){
									expr += x[v][idA];						
								}
							}
						}
						for (a=0;a<inst.inRegionArcs[v][j][t].getSize();a++){ //Arcs arriving at region through node v_j,t
							int idA = inst.inRegionArcs[v][j][t][a];
							if (inst.shouldAddArc(t2S,t2F,v, idA)){
								expr += -x[v][idA];						
							}						
						}
					}
				}				
				twoPortsNoRevisit[r][v].setExpr(expr);
			}			
		}
	}	
	#endif
}
void Model::reIntegralize(IloEnv& env, Instance inst, const int& t1S, const int& t1F){
	int j,t,v,a;	
	int J = inst.numTotalPorts;
	int T = inst.t;	
	int V = inst.speed.getSize(); 
	for(v=0;v<V;v++){
		for(a=0; a<x[v].getSize(); a++){		
			if (inst.shouldIntegralizeArc(v, a, t1S, t1F)){
				model.remove(convertX[v][a]);
			}
		}		
		for(j=0;j<J;j++){
			for(t=t1S;t<t1F;t++){							
				model.remove(convertZ[v][j][t]);				
			}
		}
		#ifndef NOperateAndGo2
		//~ for(t=t1S;t<t1F;t++){							
			//~ model.remove(convertW[v][t]);
		//~ }
		#endif
	}
}

void Model::fixAllSolution(IloEnv& env,const Instance& inst){
	if((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){				
		getSolVals(env, inst);			
		//Fix solution		
		for(int v=0; v<inst.speed.getSize(); v++){
			#ifndef NFixZvar
			for(int j=0; j<inst.numTotalPorts; j++){							
				for(int t=0; t<inst.t; t++){								
					z[v][j][t].setBounds(round(zValue[v][j][t]), round(zValue[v][j][t]));					
				}
			}
			#endif
			//Fixing X variables
			for(int a=0; a<x[v].getSize(); a++){				
				x[v][a].setBounds(round(xValue[v][a]), round(xValue[v][a]));			   
			}						
			#ifndef NOperateAndGo2
			for(int t=0; t<inst.t; t++){
				w[v][t].setBounds(round(wValue[v][t]), round(wValue[v][t]));
			}
			#endif
		}		
	}else{//If the method return a infeasible solution
		//TODO
		cout << "No  feasile solution" << endl;
		exit(1);
	}
}

void Model::unFixInterval(Instance inst, const int& tS, const int& tF){
	//Unfixing z
	for(int v=0; v<inst.speed.getSize(); v++){
		#ifndef NFixZvar
		for(int j=0; j<inst.numTotalPorts; j++){							
			for(int t=tS; t<tF; t++){								
				z[v][j][t].setBounds(0, 1);					
			}
		}
		#endif
		//Unfixing X variables
		for(int a=0; a<x[v].getSize(); a++){
			int t1,t2, arcType;		
			arcType = inst.travelArcType(inst.arcs[v][a],t1,t2);	
			if( (t2 >= tS && t2 < tF) ||  // Arc times are between relaxed times and before 				
				((arcType == 3 || arcType == 4) && (t1 >= tS && t1 < tF)) ){ //Is a sink arc in the interval 
				x[v][a].setBounds(0, 1);
		   }
		}
		#ifndef NOperateAndGo2
		//~ for(int t=tS; t<tF; t++){								
			//~ w[v][t].setBounds(0, 1);					
		//~ }
		#endif
	}
}
/*Pre-req: Fix All solution*/
void Model::improvementPhase(IloEnv& env, Instance inst, const double& mIntervals, const double& timePerIter, const double& gap, const double& overlap, Timer<chrono::milliseconds>& timer_cplex,float& opt_time,
const double& timeLimit, float& elapsed_time, double& incumbent){	
	double prevObj = 1.0e+10;
	double objValue = incumbent;
	int i;			
	Timer<chrono::milliseconds> timer_LS;
	timer_LS.start();
	cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
	if (gap > 1e-04)
		cplex.setParam(IloCplex::EpGap, gap/100);
	int tS, tF;
	while((prevObj - objValue > 0.0001) && elapsed_time/1000 < timeLimit){				
		for(i=1;i<=ceil(mIntervals*(1+overlap/100));i++){
			if(i==1)
				tS = inst.t/mIntervals*(i-1);							
			else
				tS = inst.t/mIntervals*(i-1)*(1-overlap/100);						
			tF = min(tS+inst.t/mIntervals, (double)inst.t);
			cout << "Unfixing " << tS << "..." << tF << endl;
			unFixInterval(inst, tS, tF);
			
			cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
			timer_cplex.start();
			cplex.solve();
			opt_time += timer_cplex.total();
			elapsed_time += timer_LS.total();			
			incumbent = cplex.getObjValue();
			if(elapsed_time/1000 >= timeLimit){
				cout << "Reached LS time limit " << timeLimit << ": " << elapsed_time/1000 << endl;
				break;
			}
			timer_LS.start();						
			cout << "Objective: " << incumbent << endl;
			//~ if (i < ceil(mIntervals*(1+overlap/100))){
				cout << "Re-fixing " << tS << "..." << tF;
				fixSolution(env, inst, tS, tF,0);	 //Does not fix the last interval either the sink arcs
				cout << ". Done! " << endl;
			//~ }			
		}
		prevObj = objValue;
		objValue = incumbent;		
		if(elapsed_time/1000 >= timeLimit){
			cout << "Reached LS time limit " << timeLimit << ": " << elapsed_time/1000 << endl;
			break;
		}				
	}
}

void Model::getSolVals(IloEnv& env, const Instance& inst){
	if((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){				
		for(int v=0;v<inst.speed.getSize();v++){
			xValue[v] = IloNumArray(env, xValue.getSize());			
			cplex.getValues(x[v], xValue[v]);			
			for(int j=0; j<inst.numTotalPorts; j++){				
				zValue[v][j] = IloNumArray(env, zValue[v][j].getSize());				
				cplex.getValues(z[v][j], zValue[v][j]);
			}
			#ifndef NOperateAndGo2
			//~ wValue[v] = IloNumArray(env, wValue.getSize());			
			//~ cplex.getValues(w[v], wValue[v]);			
			#endif
		}
	}else{
		cout << "Impossible to get feasible solution values!" << endl;
		exit(1);
	}
}
void Model::getSolValsW(IloEnv& env, Instance inst, const int& tS, const int& tF){
	if((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){				
		for(int v=0;v<inst.speed.getSize();v++){
			xValue[v] = IloNumArray(env, xValue[v].getSize());			
			for(int a=0;a<x[v].getSize();a++){
				if(inst.shouldAddArc(tS,tF,v,a)) 
					xValue[v][a] = cplex.getValue(x[v][a]);
			}			
			for(int j=0; j<inst.numTotalPorts; j++){
				zValue[v][j] = IloNumArray(env, zValue[v][j].getSize());				
				for(int t=tS; t<tF; t++){
					zValue[v][j][t] = cplex.getValue(z[v][j][t]);
				}				
				
			}
			#ifndef NOperateAndGo2
			//~ wValue[v] = IloNumArray(env, wValue.getSize());			
			//~ cplex.getValues(w[v], wValue[v]);			
			#endif
		}
	}else{
		cout << "Impossible to get feasible solution values!" << endl;
		exit(1);
	}
}

/*
 * Fix all x variables and allow only z variables to be dicided
 */
void Model::polish(IloEnv& env, Instance inst, const double& timeLimit,const double& gap){
	cplex.setParam(IloCplex::TiLim, timeLimit);
	if (gap > 1e-04)
		cplex.setParam(IloCplex::EpGap, gap/100);
	
	getSolVals(env, inst);	
	
	cout << "Fixiing Z and allowing X..." << endl;	
	//Unfixing X variables
	for(int v=0; v<inst.speed.getSize(); v++){
		for(int a=0; a<x[v].getSize(); a++){		
			x[v][a].setBounds(0, 1);
		}
	}
	
	//Fixing z variables	
	for(int v=0; v<inst.speed.getSize(); v++){
		for(int j=0; j<inst.numTotalPorts; j++){							
			for(int t=0; t<inst.t; t++){											
				z[v][j][t].setBounds(zValue[v][j][t], zValue[v][j][t]);								
			}
		}
	}		
	cplex.solve();	
	
	cout << "Fixiing X and allowing Z..." << endl;
	//Gets the X values
	for(int v=0; v<inst.speed.getSize(); v++){
		xValue[v] = IloNumArray(env, xValue.getSize());			
		cplex.getValues(x[v], xValue[v]);
	}
	//Fixing X variables
	for(int v=0; v<inst.speed.getSize(); v++){
		for(int a=0; a<x[v].getSize(); a++){		
			x[v][a].setBounds(xValue[v][a], xValue[v][a]);
		}
	}
	
	//Unfixing z variables	
	for(int v=0; v<inst.speed.getSize(); v++){
		for(int j=0; j<inst.numTotalPorts; j++){							
			for(int t=0; t<inst.t; t++){											
				z[v][j][t].setBounds(0, 1);								
			}
		}
	}		
	cplex.solve();
}

void Model::regionLocalSearch(IloEnv env, Instance inst, const double& timePerIter, const int& gap, Timer<chrono::milliseconds>& timer_cplex,float& opt_time,
const double& timeLimit, float& elapsed_time, double& incumbent){
	Timer<chrono::milliseconds> timer_LS;
	timer_LS.start();	
	unsigned int v,t,a, idPort;
	int i,r,j;
	double objValue = incumbent;
	double prevObj = 1.0e+10;
	
	cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
	if (gap > 1e-04)
		cplex.setParam(IloCplex::EpGap, gap/100);
	
	//~ //Get solution values
	//~ getSolVals(env, inst);
	
	//~ //FixSolution
	//~ fixAllSolution(env, inst);
	
	while((prevObj - objValue > 0.0001) && elapsed_time/1000 < timeLimit){
		prevObj = objValue;
		for(i=1;i>=0;i--){ //Type region (first allow discharging region)
			cout << "TYPE REGION = " << i << endl;
			for(r=0;r<inst.identifyPort[i].getSize();r++){ //For each region				
				for(j=0;j<inst.identifyPort[i][r].getSize();j++){ //Port					
					idPort = inst.identifyPort[i][r][j];						
					#ifndef NFixZvar
					//Allow Z variables to be solved					
					for(v=0;v<inst.speed.getSize();v++){
						for(t=0;t<inst.t;t++){							
							z[v][idPort][t].setBounds(0,1);
						}
					}				
					#endif						
					//Allow X variables to be solved
					for(v=0;v<inst.speed.getSize();v++){
						for(a=0;a<x[v].getSize();a++){
							int j1,j2,t1,t2,type;
							type = inst.getArcType(inst.arcs[v][a], t1, t2, j1, j2);
							//~ if (type != 0 && ((i==0 && type != 2 && type != 4) || (i==1 && type != 1 && type != 3)) ) //Do not allow incoming regional type arcs
							if (type != 0 && (j1 == idPort || j2 == idPort))
								x[v][a].setBounds(0,1);							
						}
					}
				}					
			}
			//Solve
			cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
			timer_cplex.start();
			cplex.solve();		
			opt_time += timer_cplex.total();
			elapsed_time += timer_LS.total();
			objValue = cplex.getObjValue();
			cout << "Objective: " << objValue << endl;
			if(elapsed_time/1000 >= timeLimit){				
				break;
			}			
			timer_LS.start();
						
			//Get solution values
			getSolVals(env, inst);			
			
			//Re-fix ports
			for(r=0;r<inst.identifyPort[i].getSize();r++){ //For each region
				for(j=0;j<inst.identifyPort[i][r].getSize();j++){ //Port					
					idPort = inst.identifyPort[i][r][j];					
					#ifndef NFixZvar
					for(v=0;v<inst.speed.getSize();v++){
						for(t=0;t<inst.t;t++){
							z[v][idPort][t].setBounds(round(zValue[v][idPort][t]),round(zValue[v][idPort][t]));
						}
					}
					#endif
					//Fix X variables										
					for(v=0;v<inst.speed.getSize();v++){						
						for(a=0;a<inst.arcs[v].size();a++){		
							int j1,j2,t1,t2,type;
							type = inst.getArcType(inst.arcs[v][a], t1, t2, j1, j2);
							if (type != 0 && (j1 == idPort || j2 == idPort)){								
								x[v][a].setBounds(round(xValue[v][a]), round(xValue[v][a]));								
							}
						}
					}					
				}				
			}			
		}		
	}
	incumbent = objValue;
}

void Model::printSolution(IloEnv env, Instance& inst){
	//Verify the solution	
	unsigned int a, j, t, v;
	unsigned int V = inst.speed.getSize();	
	unsigned int J = inst.numTotalPorts;
	unsigned int T = inst.t;
	unordered_map<string,double>::const_iterator it;
	IloArray<IloNumArray> xVals(env, V); 
	IloArray<IloNumArray> sVVals(env, V); 
	#ifndef NSpotMarket
	IloArray<IloNumArray> alphaVals(env, J); 
	#endif
	#ifndef NBetas
	IloArray<IloNumArray> betaVals(env, J); 
	#endif
	IloArray<IloNumArray> sPVals(env, J); 
	IloArray<IloArray<IloNumArray> > zVals(env, V);
	IloArray<IloArray<IloNumArray> > fVals(env, V);
	
	for(v=0;v<V;v++){
		xVals[v] = IloNumArray(env, inst.arcs[v].size());
		sVVals[v] = IloNumArray(env, T);
		cplex.getValues(x[v], xVals[v]);
		cplex.getValues(sV[v], sVVals[v]);
		zVals[v] = IloArray<IloNumArray>(env, J);
		fVals[v] = IloArray<IloNumArray>(env, J);
		for(j=0;j<J;j++){
			zVals[v][j] = IloNumArray(env,T);
			fVals[v][j] = IloNumArray(env,T);
			cplex.getValues(z[v][j],zVals[v][j]);
			cplex.getValues(f[v][j],fVals[v][j]);
		}
	}
	for(j=0;j<J;j++){
		#ifndef NSpotMarket
		alphaVals[j] = IloNumArray(env, T);
		cplex.getValues(alpha[j], alphaVals[j]);
		#endif
		#ifndef NBetas
		betaVals[j] = IloNumArray(env, T);
		cplex.getValues(beta[j], betaVals[j]);
		#endif
		sPVals[j] = IloNumArray(env, T);
		cplex.getValues(sP[j], sPVals[j]);
	}

	//Calculate the objective value	
	//Separated costs
	double arcsCost=0, arcsCost2=0, revenues1=0, revenue2=0, finishEarly=0, finishEarly2=0, spot=0, spot2=0;	
	unsigned int j1,t1,j2,t2;
	cout << "SOLUTION LOG: \n";
	for(v=0;v<V;v++){
		cout << "Vessel " << v << " route: " << endl;
		stringstream ss;
		unsigned int count = 0;
		for(a=0;a<xVals[v].getSize();a++){
			if(xVals[v][a] > 0.01){				
				it = inst.c_va[v].find(inst.arcs[v][a]);								
				arcsCost += it->second*xVals[v][a];				
				arcsCost2 += it->second;				
				if(inst.arcs[v][a].size() == 10){ //Travel Arc
					j1 = stoi(inst.arcs[v][a].substr(0,2));
					t1 = stoi(inst.arcs[v][a].substr(2,3));
					j2 = stoi(inst.arcs[v][a].substr(5,2));
					t2 = stoi(inst.arcs[v][a].substr(7,3));
					cout << " (" << j1 << ","<< t1 << ") Stock: " << sPVals[j1][t1+1] <<
					" -> (" << j2 << "," << t2  << ") Stock " << sPVals[j2][t2+1] <<
					" | Vessel cargo " << sVVals[v][t1+1] << "->" << sVVals[v][t2+1] <<
					" Arc cost " << it->second*xVals[v][a];
					if(zVals[v][j2][t2] > 0.1){
						finishEarly += t2*inst.perPeriodRewardForFinishingEarly;
						revenues1 += fVals[v][j2][t2]*inst.r_jt[j2][t2];
						if (inst.typePort[j2] == 0)
							cout << " Loaded: " << fVals[v][j2][t2];
						else
							cout << " Discharged: " << fVals[v][j2][t2];
						if (fVals[v][j2][t2] < inst.f_min_jt[j2][t2] || inst.f_max_jt[j2][t2] - fVals[v][j2][t2] < zero)
							cout << "Exceded Min/Max capacity of Load(unload) " << inst.f_min_jt[j2][t2] << "(" << inst.f_max_jt[j2][t2] << ")" << endl;						 
						else
							cout << endl;
					}else cout << endl;
					if(inst.idRegion[j1] != inst.idRegion[j2]){
						count++;
						ss << "\n Travel " << count << " visits ports: " << j2 << " ";
					}else 
						ss  << j2 << " ";
					
				}else{ //Source or sink arc
					j1 = stoi(inst.arcs[v][a].substr(1,2));
					t1 = stoi(inst.arcs[v][a].substr(3,3));						
					if (stoi(inst.arcs[v][a].substr(0,1)) == 0){
						cout << " Start at (" << j1 << "," << t1 << ")" << "with Stock " << sPVals[j1][t1] << "/" << sPVals[j1][t1+1] << " | Vessel cargo " << sVVals[v][t1+1];
						ss << " Travel " << count << " visits ports: " << j1 << " ";
						if(zVals[v][j1][t1] > 0.1){
							finishEarly += t1*inst.perPeriodRewardForFinishingEarly;
							revenues1 += fVals[v][j1][t1]*inst.r_jt[j1][t1];
							if (inst.typePort[j1] == 0)
								cout << " Loaded: " << fVals[v][j1][t1];
							else
								cout << " Discharged: " << fVals[v][j1][t1];
						}
					}else
						cout << " End at (" << j1 << "," << t1 << ")" << "with Stock " << sPVals[j1][t1] << "/" << sPVals[j1][t1+1] << " | Vessel cargo " << sVVals[v][t1+1];					
					cout << endl;
									
				}
			}
		}
		cout << ss.str() << endl;
	}
	
	//Ports info
	double totalBeta;
	IloArray<IloIntArray> minOperation(env,J); //Caluclate the minimum number of operations that must be ocur at port j from begin of time horizon until time t
	IloArray<IloIntArray> numOperation(env,J); //Count the operations made to port j from the start of time horizon until time t
	for (j=0;j<J;j++){
		numOperation[j] = IloIntArray(env,T);
		minOperation[j] = IloIntArray(env,T);
		for (t=0;t<T;t++){			
			#ifndef NSpotMarket
			spot2 += alphaVals[j][t]*inst.p_jt[j][t];
			cout << "Port " << j << " ("<<t<<"): " << sPVals[j][t] ;
			//~ cout << " Revenue " << inst.r_jt[j][t];					
			if (alphaVals[j][t] > 0.1){
				 spot += alphaVals[j][t]*inst.p_jt[j][t];
				 cout << " +- spot " << alphaVals[j][t];
				 if (alphaVals[j][t] > inst.alp_max_jt[j][t]) cout << "- Exceeded! Max value is: " << inst.alp_max_jt[j][t];
			 }
			#endif
			if (sPVals[j][t] < inst.sMin_jt[j][0] || (inst.sMax_jt[j][0]-sPVals[j][t]) <= - 0.01 ) cout << " ERROR! Port stock out of bounds";			
			
			
			//Calculating the minimum number of operations
			double sum_d_jt=0, denominator=0, max_dj=0;
			for(int t1=0;t1<=t;t1++){
				sum_d_jt += inst.d_jt[j][t1];
				double dj = inst.sMax_jt[j][t1] + inst.d_jt[j][t1] - inst.sMin_jt[j][t1];
				if (max_dj < dj)
					max_dj = dj;				
			}
			denominator = min(inst.f_max_jt[j][t], max_dj);			
			//Define the numerator according to the type of port
			if (inst.typePort[j] == 0) //Loading port
				sum_d_jt += inst.s_j0[j] - inst.sMax_jt[j][t];
			else //Discharging port
				sum_d_jt += inst.sMin_jt[j][t] - inst.s_j0[j];							
			minOperation[j][t] = ceil(sum_d_jt/denominator);
			
			//Counting number of operations
			if(t>0)
				numOperation[j][t] += numOperation[j][t-1];
			for(v=0;v<V;v++){
				if (zVals[v][j][t] > 0.1)
					numOperation[j][t]++;					
			}
			cout << " Min(executed) operations: " << minOperation[j][t] << "(" << numOperation[j][t] << ")";						
			cout << endl;	
		}
		#ifndef NBetas
		totalBeta += IloSum(betaVals[j]);
		#endif
	}
	//Vessel info
	for (v=0;v<V;v++){
		for (t=0;t<=T;t++){
			cout << "Vessel " << v << " ("<<t<<"): " << sVVals[v][t] ;
			if (sVVals[v][t] < zero || (inst.q_v[v] - sVVals[v][t] < zero)) cout << " ERROR! Vessel stock out of bounds";
			cout << endl;
		}
	}
	
	for(v=0;v<V;v++){
		for(j=0;j<J;j++){
			for(t=0;t<T;t++){				
				if (zVals[v][j][t] > 0.1){
					cout << "Z[" << v << "][" << j << "][" << t << "] = " << inst.r_jt[j][t] << " * " << fVals[v][j][t] << endl;
					revenue2 += inst.r_jt[j][t] * fVals[v][j][t];
					finishEarly2 += inst.perPeriodRewardForFinishingEarly*t;
				}
			}
		}
	}		
	
	cout << "Solution value " << cplex.getObjValue() << endl
	<< "Names1 refers to coefficent*variable values / Names2 refers to coefficent values only annd(or) computing all vector values " << endl
	<< "Arcs1: " <<  arcsCost << endl 
	<< "Arcs2: " <<  arcsCost2 << endl 
	<< "revenues1 " << revenues1 << endl
	<< "revenues2 " << revenue2 << endl
	<< "Discount per 'delayed operations'1 " << finishEarly << endl
	<< "Discount per 'delayed operations'2 " << finishEarly2 << endl
	<< "Spot1 " << spot << endl
	<< "Spot2 " << spot2 << endl
	#ifndef NBetas
	<< "Total betas " << totalBeta << endl
	#endif
	<< "Total1 " << (arcsCost + finishEarly + spot) - revenues1  << endl
	<< "Total2 " << (arcsCost2 + finishEarly2 + spot2) - revenue2  << endl;
}

void mirp::fixAndRelax(string file, string optStr, const double& nIntervals, const double& gapFirst, const int& f, const double& overLap, const int& endBlock,
const int& timePerIterFirst, const double& mIntervals, const int& timePerIterSecond, const double& gapSecond, const double& overlap2, const double& timeLimit){
	///Time parameters
	Timer<chrono::milliseconds> timer_cplex;
	Timer<chrono::milliseconds> timer_global;	
	Timer<chrono::milliseconds> timer_1stPhase;	
	timer_global.start();
	float global_time {0};
	float opt_time {0};
	float elapsed_time {0};
	float elapsed_time_it {0};
	
	///Read input files
	IloEnv env;
	
	try{
		Instance inst(env);
		inst.readInstance(env, file);					
		//Verify if the number of intervals is compatible with the instance		
		if(fmod((double)inst.t ,nIntervals) != 0 || fmod((double)inst.t,mIntervals != 0)){
			cout << "Error: the number of intervals n or m must be a divisor of Time instace without rest" << endl;
			exit(1);
		}
		
		int v, j, t, t1S, t1F, t2S, t2F, t3S, t3F;
		int J = inst.numTotalPorts;
		int T = inst.t;
		int V = inst.speed.getSize(); //# of vessels		
		double obj1stPhase = 0, obj2ndPhase = 0, time1stPhase=0, time2ndPhase=0, intPart;
		
		/// NEW MODEL
		timer_1stPhase.start();
		Model model(env);		
		model.buildFixAndRelaxModel(env,inst, nIntervals, endBlock); 	
		model.setParameters(env, timePerIterFirst, gapFirst);				
		//Rolling horizon
		#ifndef NRollingHorizon
		double p = T/nIntervals*(1-overLap/100); // Units of t that are add at each iteration to the model.
		int s = T-(T/nIntervals*endBlock); 		 // Last t of model when starting relax-and-fix.
		for(v=1; v<= ceil(T-s/p); v++){		
			cout << "ITERATION " << v << "/" << ceil(T-s/p) << endl;
			timer_cplex.start();		
			cout << "Solving...\n";
			if(!model.cplex.solve())
				cout << model.cplex.getStatus() << endl;
				
			opt_time += timer_cplex.total();
			cout << "Solution Status " << model.cplex.getStatus() << " Value: " << model.cplex.getObjValue() << endl;
			
			t3S = ceil(T/nIntervals * (v-1) * (1-overLap/100));
			t3F = min(ceil(T/nIntervals * v * (1-overLap/100)),(double)T);
			cout << "Fixing interval " << t3S << "..." << t3F << endl ;
			model.fixSolution(env, inst, t3S, t3F,0);
			cout << "Fixed! \n";	
			
			t2S = ceil(s+p*(v-1));
			t2F = min(T,(int)ceil(s+p*v)); 
						
			if(t2S < T){
				cout << "Adding to the model " << t2S << "..." << t2F << endl;			
				model.decreaseEndBlock (env, inst, nIntervals, t2S, t2F);
			}
			double newGap = max(0.001, (gapFirst - (gapFirst/ceil(T-s/p))*v) / 100);
			cout << "New GAP " << newGap*100 << " %" << endl;
			model.cplex.setParam(IloCplex::EpGap, newGap);		
		}
		#endif
		#ifdef NRollingHorizon
		//Relax-and-fix
		double p = T/nIntervals*(1-overLap/100); // Units of t that are add at each iteration to the model.
		int s = T-(T/nIntervals*endBlock); 		 // Last t (relaxed) of model when starting relax-and-fix.
		int sizeInterval = T/nIntervals;		 // Last t of integer block (always starting with 1 integer interval)
		for (v=1; v <= ceil(nIntervals); v++){				
		//~ for (v=1; v <= ceil(T-sizeInterval/p); v++){				
			timer_cplex.start();		
			//~ model.cplex.exportModel("Saida.lp");		
			cout << "Solving..." << endl;
			if(!model.cplex.solve())
				cout << model.cplex.getStatus() << endl;
				
			opt_time += timer_cplex.total();
					
			cout << "Solution Status " << model.cplex.getStatus() << " Value: " << model.cplex.getObjValue() << endl;
			//~ model.printSolution(env, inst);
			
			t3S = ceil(T/nIntervals * (v-1) * (1-overLap/100));
			t3F = min(ceil(T/nIntervals * v * (1-overLap/100)),(double)T);
			//~ t3S = ceil(T/nIntervals * (v-1) * (1-overLap/100));
			//~ t3F = min(ceil(T/nIntervals * v * (1-overLap/100)),(double)T);
			cout << "t3 (Fixing interval) " << t3S << "..." << t3F << endl ;
			model.fixSolution(env, inst, t3S, t3F,0);
			cout << "Fixed! \n";	
			
			t2S = T - (T/nIntervals* max(0.0, endBlock - modf(nIntervals, &intPart) - (v-1)*f)); //modf return the fractional part of nIntervals
			t2F = T - (T/nIntervals* max(0.0, endBlock - modf(nIntervals, &intPart) - (v*f)));
			
			//~ t2S = ceil(s+p*(v-1));
			//~ t2F = ceil(s+p*v); 
								
			if(t2S < T){
				cout << "t2 (Adding to the model) " << t2S << "..." << t2F << endl;			
				model.decreaseEndBlock (env, inst, nIntervals, t2S, t2F);
			}			
			t1S = T/nIntervals * v;
			t1F = min(T/nIntervals * (v + 1), (double) T);		
			//~ t1S = ceil(sizeInterval+p*(v-1));
			//~ t1F = ceil(sizeInterval+p*v); 
			if (t1S < T){
				cout << "t1 (Integralizing) " << t1S << "..." << t1F << endl;
				model.reIntegralize(env, inst, t1S, t1F);				
			}
			double newGap = max(0.0, (gapFirst - (gapFirst/ceil(nIntervals))*v) / 100);
			//~ double newGap = max(0.001, (gapFirst - (gapFirst/ceil(T-sizeInterval/p))*v) / 100);
			cout << "New GAP " << newGap*100 << " %" << endl;
			model.cplex.setParam(IloCplex::EpGap, newGap);		
		}
		#endif
		//Last iteration (after fixing the penultimate interval)
		timer_cplex.start();
		model.cplex.solve();
		opt_time += timer_cplex.total();
		time1stPhase = timer_1stPhase.total();
		global_time += timer_global.total();
		timer_global.start();		
		obj1stPhase = model.cplex.getObjValue();		
		double incumbent = obj1stPhase;
		cout << "Solution Status " << model.cplex.getStatus() << " Value: " << obj1stPhase << endl;
		
		cout << "Improving solution..." << endl;	
		model.fixAllSolution(env, inst);
		
		double tLimit=0;

		model.fixAndOptTW(env, inst, mIntervals, timePerIterSecond, gapSecond, overlap2, timer_cplex, opt_time, timeLimit/2, elapsed_time, incumbent);
		
		//~ model.improvementPhase(env, inst, mIntervals, timePerIterSecond, gapSecond, overlap2, timer_cplex, opt_time, timeLimit/3, elapsed_time, incumbent);

		//~ tLimit = (timeLimit - elapsed_time/1000)/2;
		//~ model.fixAndOptmizeH(env, inst, timePerIterSecond, gapSecond, incumbent, timer_cplex, opt_time, tLimit, elapsed_time);
		
		model.regionLocalSearch(env, inst,timePerIterSecond, gapSecond, timer_cplex, opt_time, timeLimit, elapsed_time, incumbent);
		
		//~ cout << "Pool of solutions : " << model.cplex.getSolnPoolNsolns() << endl;
		//~ cout << "Polishing" << endl;
		//~ model.cplex.setParam(IloCplex::TiLim, 0.5);
		//~ model.cplex.setParam(IloCplex::TiLim, timePerIterFirst);
		//~ model.cplex.setParam(IloCplex::PolishAfterTime, timePerIterSecond); 
		///TESTE - warmstart
		//~ model.cplex.setParam(IloCplex::TiLim, 36000);
		//~ model.cplex.setParam(IloCplex::NodeLim, 1);
		//UnFix solution		
		//~ for(int v=0; v<inst.speed.getSize(); v++){
			//~ #ifndef NFixZvar
			//~ for(int j=0; j<inst.numTotalPorts; j++){							
				//~ for(int t=0; t<inst.t; t++){									
					//~ model.z[v][j][t].setBounds(0,1);					
				//~ }
			//~ }
			//~ #endif
			//~ //UnFixing X variables
			//~ for(int a=0; a<model.x[v].getSize(); a++){				
				//~ model.x[v][a].setBounds(0,1);			   
			//~ }						
		//~ }
		//~ timer_cplex.start();
		//~ model.cplex.solve();
		//~ opt_time += timer_cplex.total();
	
		global_time += timer_global.total();		
		time2ndPhase = elapsed_time;//global_time - time1stPhase;
		obj2ndPhase	= incumbent;//model.cplex.getObjValue();
		
		//For getting information about solution
		model.cplex.setParam(IloCplex::TiLim, 1000);
		model.cplex.setParam(IloCplex::NodeLim, 1);
		model.cplex.solve();
		#ifndef NBetas
		double totalBeta=0;
		IloArray<IloNumArray> betaVals(env, J); 
		for(j=0;j<J;j++){
			betaVals[j] = IloNumArray(env, T);
			model.cplex.getValues(model.beta[j], betaVals[j]);
			totalBeta += IloSum(betaVals[j]);
		}		
		#endif
		
		model.printSolution(env, inst);
		cout << endl
		<< nIntervals << "\t"
		<< endBlock << "\t"
		<< overLap << "\t" 
		<< timePerIterFirst << "\t" 
		<< gapFirst << "\t"
		<< mIntervals << "\t"
		<< timePerIterSecond << "\t" 
		<< timeLimit << "\t"
		<< gapSecond << "\t"
		<< opt_time/1000 << "\t"
		<< global_time/1000 << "\t"
		<< time1stPhase/1000 << "\t"
		<< time2ndPhase/1000 << "\t"
		<< obj1stPhase << "\t"
		<< obj2ndPhase << "\t"
		<< abs((obj2ndPhase/obj1stPhase - 1)*100) << "\t"
		//<< model.cplex.getBestObjValue() << "\t"		
		<< endl;
		#ifndef NBetas
		cout << "Total betas = " << totalBeta;
		#endif
		
		
		//~ cout << "First phase -> n : " << nIntervals << " GAP: " << gapFirst << "; f : " << f << "; Overlap: " << overLap << "%; |EndBlock| " << endBlock << "; Time per block " << timePerIterFirst << endl
		//~ << "Second phase -> m: " << mIntervals << "; GAP " << gapSecond << "; Time per local search " << timePerIterSecond << endl
		//~ << "CPLEX time: " << opt_time/1000 << endl << "Other times: " << (global_time-opt_time)/1000 << endl
		//~ << "Total time: " << global_time/1000 << endl
		//~ << "1st phase time: " << time1stPhase/1000 << endl
		//~ << "2nd phase time: " << time2ndPhase/1000 << endl
		//~ << "1st phase obj: " << obj1stPhase << endl
		//~ << "2nd phase obj: " << obj2ndPhase << endl
		//~ << "Improvment from 1st to 2nd phase: " << abs((obj2ndPhase/obj1stPhase - 1)*100) << "%" << endl;
		//~ << "Final GAP " << model.cplex.getMIPRelativeGap() 
		//~ << endl;
	}catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;		
		e.end();
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}
	env.end();
}
