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
#define NEndHorizon
//~ #define NSpotMarket
//~ #define NTravelAtCapacity
//~ #define NTravelEmpty
//~ #define NBerthLimit

ILOSTLBEGIN

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloNumArray> NumNumMatrix;

using namespace std;
using mirp::Model;
using mirp::Instance;

#define zero -0.0001

void Model::buildFixAndRelaxHVModel(IloEnv& env,Instance inst, const double& nIntervals, const int& endBlock, const int& outVessels){
	int j,t,v,a;
	int timePerInterval = inst.t/nIntervals;
	int J = inst.numTotalPorts;
	int T = inst.t;
	int tOEB = inst.t - (timePerInterval*endBlock) ; //Time periods out of End Block (index tOEB is the first in the end block)
	int V = inst.speed.getSize(); //# of vessels
	int vesselsOnModel = V-outVessels;	
	unordered_map<string,double>::const_iterator it;
	
	//Init Variables, converters and storage of variable values
	#ifndef NSpotMarket
	alpha = NumVarMatrix(env, J);
	#endif
	sP = NumVarMatrix(env,J);
	f = IloArray<NumVarMatrix>(env, V);
	sV = NumVarMatrix(env,V);
	x = IloArray<IloBoolVarArray>(env, V);
	z = IloArray<IloArray<IloBoolVarArray> >(env,V);
	
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
		
		for(a=0; a<x[v].getSize(); a++){
			stringstream ss;			
			ss << "x_(" << inst.arcs[v][a] << ")," << v;
			x[v][a].setName(ss.str().c_str());			
			if(inst.shouldRelaxArc(timePerInterval, v, a) || v > 0){ //Kepp integrality only in the first interval of first vessel
				convertX[v][a] = IloConversion(env,x[v][a], ILOFLOAT); 
				model.add(convertX[v][a]);
			}
		}
		for(t=0;t<T+1;t++){
			stringstream ss;			
			ss << "supplyOnVessel_" << v << "," << t;
			if (t > inst.firstTimeAv[v])
				sV[v][t] = IloNumVar(env, 0, inst.q_v[v], ss.str().c_str());
			else
				sV[v][t] = IloNumVar(env, inst.s_v0[v], inst.s_v0[v], ss.str().c_str());
		}		
		
		for(j=0;j<J;j++){
			f[v][j] = IloNumVarArray(env, T);
			z[v][j] = IloBoolVarArray(env, T);
			zValue[v][j] = IloNumArray(env, T);
			convertZ[v][j] = IloArray<IloConversion> (env, T);	
			for(t=0;t<T;t++){
				if (t >= timePerInterval || v > 0){ //Kepp integrality only in the first interval and vessel
					convertZ[v][j][t] = IloConversion (env, z[v][j][t], ILOFLOAT);	
					model.add(convertZ[v][j][t]);
				}
				stringstream ss;
				ss << "f_(" << j << "," << t << ")," << v;
				int maxAmountOperation = min(inst.q_v[v], inst.f_max_jt[j][t]);
				f[v][j][t] = IloNumVar(env, 0, maxAmountOperation, ss.str().c_str());				
				
				ss.str(string());
				ss << "z_(" << j << "," << t << ")," << v;
				z[v][j][t].setName(ss.str().c_str());	
				#ifndef UselessVariables
				if(inst.inArcs[v][j][t].getSize() < 1){
					f[v][j][t].setBounds(0,0);
					z[v][j][t].setBounds(0,0);
				}
				#endif							
			}			
		}
	}
	
	for(j=0;j<J;j++){
		#ifndef NSpotMarket
		//~ alpha[j] = IloNumVarArray(env,0,inst.alp_max_jt[j]);
		alpha[j] = IloNumVarArray(env,T,0,IloInfinity); //Spotmarket is initially unlimited
		#endif
		sP[j] = IloNumVarArray(env, inst.sMin_jt[j], inst.sMax_jt[j]);
		sP[j].add(IloNumVar(env, inst.sMin_jt[j][0], inst.sMax_jt[j][0])); //Add more one variable to represent inventory at source node (note: using index 0 instead ind T does not interfer in the current tested instances)
		for(t=0;t<T+1;t++){
			stringstream ss;
			#ifndef NSpotMarket
			if(t<T){
				ss << "alpha_(" << j << "," << t << ")";
				alpha[j][t].setName(ss.str().c_str());
			}
			#endif
			ss.str(string());
			ss << "sP_(" << j << "," << t << ")";			
			sP[j][t].setName(ss.str().c_str());			
			
		}				
		//Fixing the initial inventory at each port
		sP[j][0].setBounds(inst.s_j0[j],inst.s_j0[j]); 
		stringstream ss;
		ss << "sP_(" << j << "," << 0 << ")";
		sP[j][0].setName(ss.str().c_str());
	}
	
	///Objective function - Just add variables z and x for the n-e interval and vesselsOnModel
	IloExpr expr(env);
	IloExpr expr1(env);
	for(v=0;v<vesselsOnModel;v++){
		for(j=0;j<J;j++){
			for(t=0;t<tOEB;t++){
				#ifndef UselessVariables				
				if (inst.inArcs[v][j][t].getSize() > 0){ 
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
		}		
	}	
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
					expr.clear();
					if(v < vesselsOnModel){
						#ifndef UselessVariables
						if (inst.inArcs[v][j][t].getSize() > 0){ //Only if there is a intering arc in the node
						#endif						
							for (a=0;a<inst.outArcs[v][j][t].getSize();a++){						
								int idA = inst.outArcs[v][j][t][a];
								//Verify if arc and node are in the model
								if ( inst.isInModel(tOEB, v, idA) && t < tOEB)
									expr += x[v][idA];						
							}
							for (a=0;a<inst.inArcs[v][j][t].getSize();a++){
								int idA = inst.inArcs[v][j][t][a];						
								//Verify if arc and node is in model
								if ( inst.isInModel(tOEB, v, idA) && t < tOEB)
									expr += -x[v][idA];		
							}							
						#ifndef UselessVariables
						}
						#endif							
					}
					flowBalance[v][j][t].setExpr(expr);						
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
				if(v < vesselsOnModel){
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
				}else{
					flowBalance[v][j][0].setBounds(0,0);
					flowBalance[v][j][1].setBounds(0,0);
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
				for(v=0;v<vesselsOnModel;v++){
					#ifndef UselessVariables
					if (inst.inArcs[v][j][t].getSize() > 0) //Only if there is a intering arc in the node
					#endif
					expr += -f[v][j][t];
				}
				#ifndef NSpotMarket			
				expr += -alpha[j][t];
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
		model.add(portInvBalance[j]);
	}
	
	//Vessels inventory balance 
	vesselInvBalance = IloArray<IloRangeArray>(env,V);
	for(v=0;v<V;v++){
		vesselInvBalance[v] = IloRangeArray(env,T);			
		for(t=0;t<T;t++){
			expr.clear();
			stringstream ss;
			ss << "invBalanceVessel_" << v << "," << t;
			if (v < vesselsOnModel){
				if (t<tOEB){
					expr += sV[v][t];
					for(j=0;j<J;j++){
						#ifndef UselessVariables
						if (inst.inArcs[v][j][t].getSize() > 0) //Only if there is a intering arc in the node
						#endif
						expr += inst.delta[j]*f[v][j][t];
					}
				}			
				if (t < tOEB)
					vesselInvBalance[v][t] = IloRange(env, 0, -sV[v][t+1] + expr, 0, ss.str().c_str());			
				else
					vesselInvBalance[v][t] = IloRange(env, 0, expr, 0, ss.str().c_str());			
			}else
				vesselInvBalance[v][t] = IloRange(env, 0, expr, 0, ss.str().c_str());			
		}		
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
				for(v=0;v<vesselsOnModel;v++){
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
					if (t < tOEB && v < vesselsOnModel){
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
			if ( inst.isInModel(tOEB, v, a) && v < vesselsOnModel){ 
				travelAtCapacity[v][a] = IloRange(env, 0, 0, ss.str().c_str());		
				travelEmpty[v][a] = IloRange(env, 0, 0, ss1.str().c_str());		
				switch (arcType)
				{
					case 1:
						#ifndef NTravelAtCapacity
						 //timeJ1 is incresead in 1 because in vector sV time 0 is considered as source node										
						travelAtCapacity[v][a].setExpr(-sV[v][timeJ1+1] + inst.q_v[v]*x[v][a]);	
						travelAtCapacity[v][a].setBounds(-IloInfinity, 0);
						//~ cout << "Travel at capacity " << v << " arc " << inst.arcs[v][a] << endl;														
						#endif
						break;
					case 2:
						#ifndef NTravelEmpty
						travelEmpty[v][a].setExpr(inst.q_v[v]*x[v][a] + sV[v][timeJ1+1]);
						travelEmpty[v][a].setBounds(-IloInfinity, inst.q_v[v]);						
						#endif					
						break;
					case 3:
						#ifndef NTravelAtCapacity
						travelAtCapacity[v][a].setExpr(-sV[v][timeJ1+1] + inst.q_v[v]*x[v][a]);
						travelAtCapacity[v][a].setBounds(-IloInfinity, 0);
						#endif					
						break;
					case 4:
						#ifndef NTravelEmpty
						travelEmpty[v][a].setExpr(inst.q_v[v]*x[v][a] + sV[v][timeJ1+1]);						
						travelEmpty[v][a].setBounds(-IloInfinity, inst.q_v[v]);
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
		//~ cumSlack[j] = IloRange(env, expr, inst.alp_max_j[j], ss.str().c_str());
		//~ cumSlack[j] = IloRange(env, inst.alp_max_j[j]*0.2, expr, IloInfinity, ss.str().c_str()); //Minimum of spot
		cumSlack[j] = IloRange(env, expr, IloInfinity, ss.str().c_str());
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
					if (t < tOEB && v < vesselsOnModel){
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
}

void Model::fixVesselInterval(IloEnv env,Instance inst, const int& v,const int& tS, const int& tF){
	if((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){
		//Get the values				
		xValue[v] = IloNumArray(env, xValue.getSize());			
		cplex.getValues(x[v], xValue[v]);		
		//~ cout << "OK" << endl;
		
		for(int j=0; j<inst.numTotalPorts; j++){				
			zValue[v][j] = IloNumArray(env, zValue[v][j].getSize());				
			cplex.getValues(z[v][j], zValue[v][j]);			
		}					
		//Fixing z
		for(int j=0; j<inst.numTotalPorts; j++){							
			for(int t=tS; t<tF; t++){											
				z[v][j][t].setBounds(round(zValue[v][j][t]),round(zValue[v][j][t]));
			}
		}
		//Fixing X variables
		for(int a=0; a<x[v].getSize(); a++){		
			int t1,t2, arcType;		
				arcType = inst.travelArcType(inst.arcs[v][a],t1,t2);					
				if( ( (arcType != 3) && (arcType != 4) && // Is not a sink arc				
				(t2 >= tS && t2 < tF) ) )  // Arc times are between fixing times	(and travessing fixed -> new fixed)					
					x[v][a].setBounds(round(xValue[v][a]), round(xValue[v][a]));				
		}
	}else{
		cout << "Solution infeasible" << endl;
		exit(1);		
	}	
}

void Model::addVesselIntervalToModel(IloEnv& env, Instance inst, const int& vAdd, const double& nIntervals, const int& tS, const int& tF, const int& outVessels, const int& itInt){
	int i,j,t,a;
	int timePerInterval = inst.t/nIntervals;
	int J = inst.numTotalPorts;
	int T = inst.t;	
	int V = inst.speed.getSize(); //# of vessels
	unordered_map<string,double>::const_iterator it;
	
	///Objective function
	IloExpr current = cplex.getObjective().getExpr(); //Gets the current objective function
	IloExpr expr(env);
	IloExpr expr1(env);	
	for(j=0;j<J;j++){
		for(t=tS;t<tF;t++){
			#ifndef UselessVariables
			if (inst.inArcs[vAdd][j][t].getSize() > 0){ //Only if there is a intering arc in the node
			#endif
				expr += inst.r_jt[j][t]*f[vAdd][j][t];				
				expr1 += (t*inst.perPeriodRewardForFinishingEarly)*z[vAdd][j][t];
			#ifndef UselessVariables
			}
			#endif
		}		
	}		
	for (a=0;a<inst.arcs[vAdd].size();a++){
		//Add only the arc in the model that are not in the end block
		if (inst.shouldAddArc(tS, tF, vAdd, a)){
			int t1,t2;
			it = inst.c_va[vAdd].find(inst.arcs[vAdd][a]);
			#ifndef NDEBUG
			int arcType = inst.travelArcType(inst.arcs[vAdd][a], t1,t2);
			if (arcType == 0){
				cout << "Erro - source node was not in model in the first iteration" << endl;
				exit(1);				
			}				
			if (it == inst.c_va[vAdd].end()){
				cout << "Not found key for arc " << inst.arcs[vAdd][a] << endl;
				exit(1);
			}				
			#endif
			expr1 +=  it->second * x[vAdd][a];						
		}			
	}

	for(j=0;j<J;j++){
		for(t=tS;t<tF;t++){			
			#ifndef NSpotMarket
			expr1 += inst.p_jt[j][t]*alpha[j][t]*10;
			#endif
		}		
	}
	obj.setExpr(current + expr1 - expr);
	///Flow balance
	int vesselsOnModel = V-outVessels;
	for(j=0;j<J+1;j++){ //J+1 used for access the information of source and sink node
		if(j<J){ //If is a real port				
			for(t=0;t<tF;t++){ //May have arcs in constraints already in the model	(TODO Improve for not initalize t=0 ever)
				#ifndef UselessVariables
				if (inst.inArcs[vAdd][j][t].getSize() > 0){ //Only if there is a intering arc in the node
				#endif
					//Get the current expr (empty or not) and updates it with new arcs 
					expr.clear();
					expr += flowBalance[vAdd][j][t].getExpr();					
					for (a=0;a<inst.outArcs[vAdd][j][t].getSize();a++){						
						int idA = inst.outArcs[vAdd][j][t][a];
						//Verify if the arc should be added 
						
						if (inst.shouldAddArc(tS, tF, vAdd, idA)) 
							expr += x[vAdd][idA];						
					}
					for (a=0;a<inst.inArcs[vAdd][j][t].getSize();a++){
						int idA = inst.inArcs[vAdd][j][t][a];						
						//Verify if the arc should be added 
						int t1,t2;
						if (inst.shouldAddArc(tS, tF, vAdd, idA) || 
							(inst.travelArcType(inst.arcs[vAdd][idA],t1,t2) == 0 && vAdd >= vesselsOnModel && itInt == 0))  //TODO Add source arc when vessel is not in the model (first iteration) passing itint as argument
							expr += -x[vAdd][idA];		
					}										
					flowBalance[vAdd][j][t].setExpr(expr);
				#ifndef UselessVariables
				}
				#endif					
			}			
		}else{ //Is source (0 - assuming it does not need to be updated) or sink (1) node								
			expr.clear();
			expr1.clear();						
			for (a=0;a<inst.outArcs[vAdd][J][0].getSize();a++){
				//summing for source node
				int idA = inst.outArcs[vAdd][J][0][a];										
				expr += x[vAdd][idA];										
			}					
				
			for (a=0;a<inst.inArcs[vAdd][J][0].getSize();a++){
				//summing for sink nodes
				int idA = inst.inArcs[vAdd][J][0][a];					
					if ( inst.shouldAddArc(tS, tF, vAdd, idA))
					expr1 += x[vAdd][idA];		
			}
			
			//Flow Balance source node			
			flowBalance[vAdd][j][0].setExpr(expr);
			
			//Flow Balance sink node (need to be add?)			
			//~ expr1 += flowBalance[vAdd][j][1].getExpr();
			flowBalance[vAdd][j][1].setExpr(expr1);									
			
			flowBalance[vAdd][j][0].setBounds(1,1);
			flowBalance[vAdd][j][1].setBounds(1,1);
		}
	}
	//Inventory balance ports
	// A - When the constraint was already starded - just add the variable refering to the vessel vAdd
	for(j=0;j<J;j++){		
		for(t=tS;t<tF;t++){
			if (vAdd != 0){
				expr.clear();					
				#ifndef UselessVariables
				if (inst.inArcs[vAdd][j][t].getSize() > 0)
				#endif			
					expr += -f[vAdd][j][t];			
				portInvBalance[j][t].setExpr(portInvBalance[j][t].getExpr() - inst.delta[j] * expr);				
			}else{// B - 	when it is necessary to add the entire constraint, i.e. when vAdd == 0  
				expr.clear();			
				#ifndef UselessVariables
				if (inst.inArcs[vAdd][j][t].getSize() > 0) //Only if there is a intering arc in the node
				#endif			
					expr += -f[vAdd][j][t];								
				#ifndef NSpotMarket			
				expr += -alpha[j][t];
				#endif				
				portInvBalance[j][t].setBounds(inst.delta[j]*inst.d_jt[j][t], inst.delta[j]*inst.d_jt[j][t]);
				portInvBalance[j][t].setExpr(sP[j][t+1] - sP[j][t] - inst.delta[j] * expr);				
			}
		}	
	}
	//Vessels inventory balance 	
	for(t=tS;t<tF;t++){
		expr.clear();	
		expr += sV[vAdd][t];
		for(j=0;j<J;j++){
			#ifndef UselessVariables
			if (inst.inArcs[vAdd][j][t].getSize() > 0) //Only if there is a intering arc in the node
			#endif
				expr += inst.delta[j]*f[vAdd][j][t];
		}			
		vesselInvBalance[vAdd][t].setExpr(-sV[vAdd][t+1] + expr);			
	}					
	//Berth Limit
	#ifndef NBerthLimit	
	for(j=0;j<J;j++){		
		for(t=tS;t<tF;t++){
			if (vAdd != 0){ //Just update the constraint
				expr.clear();			
				expr += berthLimit[j][t].getExpr();
				#ifndef UselessVariables
				if (inst.inArcs[vAdd][j][t].getSize() > 0) 
				#endif
					expr += z[vAdd][j][t];				
				berthLimit[j][t].setExpr(expr);			
			}else{ // Add the entire constraint
				expr.clear();			
				berthLimit[j][t].setBounds(-IloInfinity, inst.b_j[j]);			
				#ifndef UselessVariables
				if (inst.inArcs[vAdd][j][t].getSize() > 0) 
				#endif
					expr += z[vAdd][j][t];				
				berthLimit[j][t].setExpr(expr);
			}
		}
	}
	#endif
	//Vessels only attempt to load/discharge at node only if the vessel is actually at the node	
	for(j=0;j<J;j++){		
		for(t=tS;t<tF;t++){						
			#ifndef UselessVariables
			if (inst.inArcs[vAdd][j][t].getSize() > 0){ //Only if there is a intering arc in the node
			#endif
				expr.clear();									
				for(a=0;a<inst.inArcs[vAdd][j][t].getSize();a++){
					int idA = inst.inArcs[vAdd][j][t][a];
					expr += x[vAdd][idA];
				}					
				atemptToOperate[j][t][vAdd].setBounds(-IloInfinity, 0);
				atemptToOperate[j][t][vAdd].setExpr(z[vAdd][j][t] - expr);			
			#ifndef UselessVariables
			}
			#endif				
		}
	}
	//Vessels must travel at capacity and empty	
	for(a=0;a<inst.arcs[vAdd].size();a++){			
		int timeJ1, timeJ2;				
		int arcType = inst.travelArcType(inst.arcs[vAdd][a], timeJ1, timeJ2);						
		switch (arcType)
		{
			case 1:
				#ifndef NTravelAtCapacity	
				//~ cout << "Travel Arc V " << v << ": " <<  inst.arcs[vAdd][a] << " T1 " << timeJ1 << " T2 " << timeJ2 << endl;							 
				if (timeJ2 >= tS && timeJ2 < tF){ 
					travelAtCapacity[vAdd][a].setBounds(-IloInfinity, 0);
					travelAtCapacity[vAdd][a].setExpr(-sV[vAdd][timeJ1+1] + inst.q_v[vAdd]*x[vAdd][a]);						
					//~ cout << "Travel at capacity " << v << " arc " << inst.arcs[vAdd][a] << endl;			
				}
				#endif
				break;
			case 2:
				#ifndef NTravelEmpty
				if (timeJ2 >= tS && timeJ2 < tF){ 
					travelEmpty[vAdd][a].setBounds(-IloInfinity, inst.q_v[vAdd]);
					travelEmpty[vAdd][a].setExpr(inst.q_v[vAdd]*x[vAdd][a] + sV[vAdd][timeJ1+1]);						
					//~ cout << "Travel empty " << v << " arc " << inst.arcs[vAdd][a] << endl;			
				}
				#endif					
				break;
			case 3:
				#ifndef NTravelAtCapacity
				if (timeJ1 >= tS && timeJ1 < tF){
					travelAtCapacity[vAdd][a].setBounds(-IloInfinity, 0);
					travelAtCapacity[vAdd][a].setExpr(-sV[vAdd][timeJ1+1] + inst.q_v[vAdd]*x[vAdd][a]);
					//~ cout << "Travel at capacity " << v << " arc " << inst.arcs[vAdd][a] << endl;			
				}
				#endif					
				break;
			case 4:
				#ifndef NTravelEmpty
				if (timeJ1 >= tS && timeJ1 < tF){ 
					travelEmpty[vAdd][a].setBounds(-IloInfinity, inst.q_v[vAdd]);
					travelEmpty[vAdd][a].setExpr(inst.q_v[vAdd]*x[vAdd][a] + sV[vAdd][timeJ1+1]);						
					//~ cout << "Travel empty " << v << " arc " << inst.arcs[vAdd][a] << endl;			
				}
				#endif
				break;
			default:
				break;
		}
	}				
	//Cumulative amount of product that can be purchased from or sold to a spot market by each port
	#ifndef NSpotMarket	
	for(j=0;j<J;j++){		
		expr.clear();
		expr = cumSlack[j].getExpr();
		for(t=tS;t<tF;t++){			
			expr += alpha[j][t];							
		}
		cumSlack[j].setExpr(expr);
	}
	model.add(cumSlack);
	#endif
	
	//Amount loaded/discharged must be in the pre-specified interval
	for(j=0;j<J;j++){
		for(t=tS;t<tF;t++){
			#ifndef UselessVariables
			if (inst.inArcs[vAdd][j][t].getSize() > 0){ //Only if there is a intering arc in the node
			#endif								
				operationLowerLimit[vAdd][j][t].setBounds(0, IloInfinity);
				operationLowerLimit[vAdd][j][t].setExpr(f[vAdd][j][t] - inst.f_min_jt[j][t]*z[vAdd][j][t]);								 	
				
				IloNum minfMax = min(inst.f_max_jt[j][t], inst.q_v[vAdd]);
				operationUpperLimit[vAdd][j][t].setBounds(-IloInfinity, 0);
				operationUpperLimit[vAdd][j][t].setExpr( f[vAdd][j][t] - minfMax*z[vAdd][j][t]); 									
			#ifndef UselessVariables
			}
			#endif
		}		
	}			
}

void Model::integralizeVesselInterval(Instance inst, const int& vInt, const int& tS, const int& tF){
	int j,t,a;	
	int J = inst.numTotalPorts;
	int T = inst.t;	
	int V = inst.speed.getSize(); 
	for(a=0; a<x[vInt].getSize(); a++){		
		if (inst.shouldIntegralizeArc(vInt, a, tS, tF)){
			model.remove(convertX[vInt][a]);
		}
	}		
	for(j=0;j<J;j++){
		for(t=tS;t<tF;t++){
			model.remove(convertZ[vInt][j][t]);				
		}
	}
}



void mirp::fixAndRelaxHV(string file, const double& nIntervals, const int& f, const double& overLap, const int& endBlock,
	const double& gapFirst, const int& outVessels, const int& timeLimitFirst, 
	const double& mIntervals, const int& timeLimitSecond, const double& gapSecond){	
	///Time parameters
	Timer<chrono::milliseconds> timer_cplex;
	Timer<chrono::milliseconds> timer_global;
	timer_global.start();
	float global_time {0};
	float opt_time {0};
	
	///Read input files
	IloEnv env;
	Instance inst(env);
	inst.readInstance(env, file);			
	int v, addVessel, vInt, t1S, t1F, t2S, t2F, t3S=0, t3F;
	int J = inst.numTotalPorts;
	int T = inst.t;
	int V = inst.speed.getSize(); //# of vessels
	double obj1stPhase = 0, obj2ndPhase = 0, time1stPhase=0, time2ndPhase=0;
	
	//Verify if the number of outVessels is too big
	if (V - outVessels < 2){
		cout << "Error: number of vessels out of model should not exced V - 2" << endl;
		exit(1);
	}	
	//Verify if the number of intervals is compatible with the instance
	if (fmod(inst.t ,nIntervals) != 0 || fmod(inst.t,mIntervals) != 0){
		cout << "Error: the number of intervals n or m must be a divisor of Time instace without rest" << endl;
		exit(1);
	}
	/// NEW MODEL
	Model model(env);		
	model.buildFixAndRelaxHVModel(env,inst,nIntervals, endBlock, outVessels);	
	//~ model.cplex.exportModel("first.lp");
	model.setParameters(env, timeLimitFirst, gapFirst);	
	unsigned int contInterval = 0;
	unsigned int contV = 0;		
	int itInt1=0,itInt2=0,itInt3=0,itV1=0,itV2=0,itV3=0;
	
	//~ for (int i=0; i <= nIntervals*V; i++){
	while (itV3 <= V && t3S < T){
		//~ model.cplex.exportModel("HorizonVerticall.lp");		
		timer_cplex.start();
		model.cplex.solve();
		opt_time += timer_cplex.total();		
		///Fixing
		if(itV3 > 1 && itV3 % V == 0){			
			itV3 = 0;
			itInt3++;
		}
		//Update GAP
		cout << "Gap " << gapFirst/(itInt3+1) << endl;
		model.cplex.setParam(IloCplex::EpGap, gapFirst/(itInt3+1)/100);
		
		t3S = floor(T/nIntervals * itInt3 * (1.0 - overLap/100));
		t3F = min(double(T), floor(T/nIntervals * (itInt3+1) * (1.0 - overLap/100)));		
		
		if (t3F == T && itV3 == V-1) break;		
		if (t3S < T){			
			cout << "Fix vessel " << itV3 <<" and interval " << t3S << "..." << t3F << endl ;
			model.fixVesselInterval(env, inst, itV3, t3S, t3F);
		}
		itV3++;		
		
		///Adding
		if((itInt2 == 0 && itV2 >= 1 && itV2 % outVessels == 0) ||
			(itV2 > 1 && itV2 % V == 0)){			
			itV2 = 0;
			itInt2++;					
		}
		
		if (itInt2 == 0){
			t2S = 0;
			t2F = floor(T/nIntervals * (nIntervals - endBlock));
		}else{
			t2S = T - floor(T/nIntervals * max(0, endBlock - (itInt2-1)*f));
			t2F = T - floor(T/nIntervals * max(0, endBlock - itInt2*f));		
		}
		
		//~ t2S = T - (T/nIntervals* max(0, endBlock - (v-1)*f));
		//~ t2F = T - (T/nIntervals* max(0, endBlock - (v*f)));		
		
		if (itInt2 == 0 && itV2 < outVessels) //Special case for first iteration
			addVessel = V - outVessels + itV2;
		else
			addVessel = itV2;
				
		if (t2S < T){
			cout << "Add vessel " << addVessel << " and interval " << t2S << "..." << t2F << endl;						
			model.addVesselIntervalToModel(env,inst,addVessel, nIntervals, t2S, t2F, outVessels, itInt2);
		}
		itV2++;
		
		///Integralizing
		if((itInt1 == 0 && itV1 > 1 && itV1 - (V-1) == 0) ||
		(itV1 > 1 && itV1 % V == 0)){			
			itV1 = 0;
			itInt1++;
		}
		t1S = T/nIntervals * itInt1;
		t1F = T/nIntervals * (itInt1 + 1);		
		if (itInt1 == 0 )//&& itV1 < outVessels+1)
			vInt = itV1 + 1;
		else
			vInt = itV1;
		if(vInt < V && t1S < T){
			cout << "Int vessel " << vInt << " and interval " << t1S << "..." << t1F << endl;			
			model.integralizeVesselInterval(inst, vInt, t1S, t1F);		
		}
		//~ model.cplex.exportModel("last.lp");
		itV1++;
		cout << endl;
	}
	
	time1stPhase = opt_time;
	obj1stPhase = model.cplex.getObjValue();
	
	///2nd phase
	cout << "Improving solution..." << endl;		
	//Vertical
	//~ model.improvementPhase(env, inst, mIntervals, timeLimitSecond, gapSecond);
	
	//Horizontal
	model.fixAllSolution(env, inst);
	//~ model.fixAndOptmizeH(env, inst, timeLimitSecond, gapSecond, obj1stPhase, timer_cplex,opt_time);
	
	//~ cout << "Polishing" << endl;
	//~ model.polish(env, inst, timeLimitSecond, gapSecond);		
	
	model.getSolVals(env, inst);		
	
	model.printSolution(env, inst);
	
	model.fixAllSolution(env, inst);
	
	//~ model.cplex.exportModel("last.lp");
	
	//Add constraints of spot market
	for(int j=0;j<J;j++){
		model.cumSlack[j].setBounds(0, inst.alp_max_j[j]);
		for(int t=0;t<T;t++){
			model.alpha[j][t].setBounds(0,inst.alp_max_jt[j][t]);
		}
	}
	
	//Reset objectiveFunction
	cout << "Reset objective function..." << endl;
	model.resetObjFunction(env, inst);	
	cout << "Re optimizing " << endl;
	//~ model.improvementPhase(env, inst, mIntervals, timeLimitSecond, gapSecond);
	//~ model.fixAndOptmizeH(env, inst, timeLimitSecond, gapSecond, obj1stPhase, timer_cplex,opt_time);
	//~ model.regionLocalSearch(env, inst,timeLimitSecond, gapSecond, timer_cplex, opt_time);
	cout << "Polishing again" << endl;
	model.polish(env, inst, timeLimitSecond, gapSecond);	
	
	obj2ndPhase	= model.cplex.getObjValue();
	
	global_time = timer_global.total();		
	time2ndPhase = global_time - time1stPhase;
	
	model.printSolution(env, inst);
	cout << " GAP Fisrt iter: " << gapFirst << " |outVessels| " << outVessels << " Time per vessel " << timeLimitFirst << endl
	<< "GAP Second phase " << gapSecond << " Time Limit second phase " << timeLimitSecond  << endl
	<< "CPLEX time: " << opt_time/1000 << endl << "Other times: " << (global_time-opt_time)/1000 << endl
	<< "Total time: " << global_time/1000 << endl
	<< "1st phase time: " << time1stPhase/1000 << endl
	<< "2nd phase time: " << time2ndPhase/1000 << endl
	<< "Improvment from 1st to 2nd phase: " << abs((obj2ndPhase/obj1stPhase - 1)*100) << "%" << endl;
	
}
