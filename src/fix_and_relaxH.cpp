#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <utility>
#include <algorithm>
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
//~ #define NPairsVessels 
#define NVesselClasses
#define NOperateAndGo
#define NOperateAndGo2
#define NBetas
//~ #define NFixZvar

ILOSTLBEGIN

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloNumArray> NumNumMatrix;

using namespace std;
using mirp::Model;
using mirp::Instance;

#define zero -0.0001
	
/* Build a model for the first iteration of fix and relax horizontally
 * All variables are created (including that are in end block)
 * Just nVessels have integer variables x and z
 * outVessel is the number of vessels that are not used in the first iteration
 * Variables that are in the end block are not included in any constraint and objective function
 * Constraints of arcs/nodes in the and block are created, but are empty 
 */
void Model::buildFixAndRelaxHorizontalModel(IloEnv& env, Instance inst, const int& outVessels){
	int j,t,v,a;	
	int J = inst.numTotalPorts;
	int T = inst.t;
	int V = inst.speed.getSize(); //# of vessels
	int vesselsOnModel = V-outVessels;	
	unordered_map<string,double>::const_iterator it;
	
	if (V - outVessels < 2){
		cout << "Erro! Number of vessels must be at least 2" << endl;
		exit(1);
	}
	
	//Sort the order that the vessels will be considered (ordered by first time avaliable)
	 for (v=0;v<V;v++){		 
		 ordV.push_back(make_pair(inst.firstTimeAv[v],v));		 
	 }
	 sort(ordV.begin(), ordV.end());
	 
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
			if(v != ordV[0].second){ //Kepp integrality only for the first vessel								
			//~ if(v > 0){ //Kepp integrality only for the first vessel
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
				if (v != ordV[0].second){ //Kepp integrality only for the first vessel
				//~ if (v > 0){ //Kepp integrality only for the first vessel					
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
		sP[j].add(IloNumVar(env, inst.sMin_jt[j][0], inst.sMax_jt[j][0])); //Add more one variable to represent inventory at source node (note: using index 0 instead ind T does not interfer in the current tested instances)
		for(t=0;t<T+1;t++){
			stringstream ss, ss1;
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
			
		}				
		//Fixing the initial inventory at each port
		sP[j][0].setBounds(inst.s_j0[j],inst.s_j0[j]); 
		stringstream ss;
		ss << "sP_(" << j << "," << 0 << ")";
		sP[j][0].setName(ss.str().c_str());
	}
	
	///Objective function - Just add variables of vesselsOnModel
	IloExpr expr(env);
	IloExpr expr1(env);
	for(v=0;v<vesselsOnModel;v++){
		for(j=0;j<J;j++){
			for(t=0;t<T;t++){
				#ifndef UselessVariables				
				if (inst.inArcs[ordV[v].second][j][t].getSize() > 0){ 
				#endif
					expr += inst.r_jt[j][t]*f[ordV[v].second][j][t];				
					expr1 += (t*inst.perPeriodRewardForFinishingEarly)*z[ordV[v].second][j][t];
				#ifndef UselessVariables
				}
				#endif
			}		
		}		
		for (a=0;a<inst.arcs[ordV[v].second].size();a++){
			it = inst.c_va[ordV[v].second].find(inst.arcs[ordV[v].second][a]);				
			if (it == inst.c_va[ordV[v].second].end()){
				cout << "Not found key for arc " << inst.arcs[ordV[v].second][a] << endl;
				exit(1);
			}				
			expr1 +=  it->second * x[ordV[v].second][a];									
		}
	}	
	for(j=0;j<J;j++){
		for(t=0;t<T;t++){			
			#ifndef NSpotMarket
			expr1 += inst.p_jt[j][t]*alpha[j][t];		
			#endif
			#ifndef NBetas
			expr1 += 1000*beta[j][t];
			#endif
		}		
	}	
	obj.setExpr(expr1-expr);	
	model.add(obj);
	
	///Constraints
	//Flow balance in port-time nodes, source and sink
	flowBalance = IloArray<IloArray<IloRangeArray> >(env,V);
	for(v=0;v<V;v++){
		flowBalance[ordV[v].second] = IloArray<IloRangeArray>(env,J+1);
		for(j=0;j<J+1;j++){ //J+1 used for access the information of source and sink node
			if(j<J){ //If is a real port
				flowBalance[ordV[v].second][j] = IloRangeArray(env,T, 0, 0);
				for(t=0;t<T;t++){	
					expr.clear();
					if(v < vesselsOnModel){
						#ifndef UselessVariables
						if (inst.inArcs[ordV[v].second][j][t].getSize() > 0){ //Only if there is a entering arc in the node
						#endif											
							for (a=0;a<inst.outArcs[ordV[v].second][j][t].getSize();a++){						
								int idA = inst.outArcs[ordV[v].second][j][t][a];							
								expr += x[ordV[v].second][idA];						
							}
							for (a=0;a<inst.inArcs[ordV[v].second][j][t].getSize();a++){
								int idA = inst.inArcs[ordV[v].second][j][t][a];													
								expr += -x[ordV[v].second][idA];		
							}							
						#ifndef UselessVariables
						}
						#endif		
					}		
					flowBalance[ordV[v].second][j][t].setExpr(expr);			
					stringstream ss;
					ss << "flowBalance_(" << j << "," << t << ")," << v;					
					flowBalance[ordV[v].second][j][t].setName(ss.str().c_str());															
				}
				model.add(flowBalance[ordV[v].second][j]);
			}else{
				flowBalance[ordV[v].second][j] = IloRangeArray(env,2, 1, 1); // Two index, 0 for the source and 1 for sink node
				expr.clear();
				expr1.clear();
				stringstream ss;
				stringstream ss1;
				ss << "flowBalanceSource" << "," << v;
				ss1 << "flowBalanceSink" << "," << v;
				if(v < vesselsOnModel){
					for (a=0;a<inst.outArcs[ordV[v].second][J][0].getSize();a++){
						//summing for source node
						int idA = inst.outArcs[ordV[v].second][J][0][a];										
						expr += x[ordV[v].second][idA];										
					}
					for (a=0;a<inst.inArcs[ordV[v].second][J][0].getSize();a++){
						//summing for sink node
						int idA = inst.inArcs[ordV[v].second][J][0][a];										
						expr1 += x[ordV[v].second][idA];		
					}
				}else{
					flowBalance[ordV[v].second][j][0].setBounds(0,0);
					flowBalance[ordV[v].second][j][1].setBounds(0,0);
				}
				//Flow Balance source node
				flowBalance[ordV[v].second][j][0].setExpr(expr);
				flowBalance[ordV[v].second][j][0].setName(ss.str().c_str());					
				//Flow Balance sink node
				flowBalance[ordV[v].second][j][1].setExpr(expr1);
				flowBalance[ordV[v].second][j][1].setName(ss1.str().c_str());
				
				//~ model.add(flowBalance[ordV[v].second][j]); //Both source and sink node flow				
				model.add(flowBalance[ordV[v].second][j][0]); //Only source node flow TODO Evaluate if it is necessary.
			}			
		}
	}
	//Inventory balance ports 
	portInvBalance = IloArray<IloRangeArray>(env,J);
	for(j=0;j<J;j++){
		portInvBalance[j] = IloRangeArray(env,T);
		for(t=0;t<T;t++){ 
			expr.clear();			
			for(v=0;v<vesselsOnModel;v++){
				#ifndef UselessVariables
				if (inst.inArcs[ordV[v].second][j][t].getSize() > 0) //Only if there is a intering arc in the node
				#endif
				expr += -f[ordV[v].second][j][t];
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
			
		}		
		model.add(portInvBalance[j]);
		//End-of-horzion: Requires that on the last time-period the inventory at ports should be between +-10% of the initial inventory
		#ifndef NEndHorizon
			model.add(sP[j][T-1] <= inst.s_j0[j] + (inst.s_j0[j]*0.3));
			model.add(sP[j][T-1] >= inst.s_j0[j] - (inst.s_j0[j]*0.3));
		#endif
	}
	
	//Vessels inventory balance
	vesselInvBalance = IloArray<IloRangeArray>(env,V);
	for(v=0;v<V;v++){
		vesselInvBalance[ordV[v].second] = IloRangeArray(env,T);			
		for(t=0;t<T;t++){
			expr.clear();			
			stringstream ss;
			ss << "invBalanceVessel_" << v << "," << t;
			if (v < vesselsOnModel){
				expr += sV[ordV[v].second][t];
				for(j=0;j<J;j++){
					#ifndef UselessVariables
					if (inst.inArcs[ordV[v].second][j][t].getSize() > 0) //Only if there is a intering arc in the node
					#endif
					expr += inst.delta[j]*f[ordV[v].second][j][t];
				}							
				vesselInvBalance[ordV[v].second][t] = IloRange(env, 0, -sV[ordV[v].second][t+1] + expr, 0, ss.str().c_str());			
			}else
				vesselInvBalance[ordV[v].second][t] = IloRange(env, 0, expr, 0, ss.str().c_str());			
		}		
		model.add(vesselInvBalance[ordV[v].second]);
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
			for(v=0;v<vesselsOnModel;v++){
				#ifndef UselessVariables
				if (inst.inArcs[ordV[v].second][j][t].getSize() > 0) //Only if there is a intering arc in the node
				#endif
					expr += z[ordV[v].second][j][t];				
			}			
			berthLimit[j][t] = IloRange(env,-IloInfinity, expr, inst.b_j[j], ss.str().c_str());			
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
				if (inst.inArcs[ordV[v].second][j][t].getSize() > 0){ 
				#endif
					if (v < vesselsOnModel){					
						for(a=0;a<inst.inArcs[ordV[v].second][j][t].getSize();a++){
							int idA = inst.inArcs[ordV[v].second][j][t][a];
							expr += x[ordV[v].second][idA];
						}					
						atemptToOperate[j][t][ordV[v].second] = IloRange(env, -IloInfinity, z[ordV[v].second][j][t] - expr, 0, ss.str().c_str());				
					}else //Just descondirer the restriction for the port time and vessel
						atemptToOperate[j][t][ordV[v].second] = IloRange(env, 0, 0, ss.str().c_str());				
				#ifndef UselessVariables
				}else{
					expr.clear();
					atemptToOperate[j][t][ordV[v].second] = IloRange(env, 0, z[ordV[v].second][j][t], 0, ss.str().c_str());				
				}
				#endif
			}
			model.add(atemptToOperate[j][t]);			
		}
	}
	
	//Vessels must travel at capacity and empty 
	travelAtCapacity = IloArray<IloRangeArray> (env, V);
	travelEmpty = IloArray<IloRangeArray> (env, V);
	expr.clear();
	for(v=0;v<V;v++){
		travelAtCapacity[ordV[v].second] = IloRangeArray(env,inst.arcs[ordV[v].second].size()); ///Altough just travel and sink arcs are incluede in this constraint set, the size is on the number of arcs
		travelEmpty[ordV[v].second] = IloRangeArray(env, inst.arcs[ordV[v].second].size());
		for(a=0;a<inst.arcs[ordV[v].second].size();a++){			
			int timeJ1, timeJ2;	
			stringstream ss, ss1;
			ss << "travelAtCap_(" << inst.arcs[ordV[v].second][a] << ")," << v;
			ss1 << "travelEmpty_(" << inst.arcs[ordV[v].second][a] << ")," << v;
			int arcType = inst.travelArcType(inst.arcs[ordV[v].second][a], timeJ1, timeJ2);			
			if (v < vesselsOnModel){
				travelAtCapacity[ordV[v].second][a] = IloRange(env, -IloInfinity, 0, ss.str().c_str());		
				travelEmpty[ordV[v].second][a] = IloRange(env, -IloInfinity, inst.q_v[ordV[v].second], ss1.str().c_str());		
				switch (arcType)
				{
					case 1:
						#ifndef NTravelAtCapacity
						 //timeJ1 is incresead in 1 because in vector sV time 0 is considered as source node										
						travelAtCapacity[ordV[v].second][a].setExpr(-sV[ordV[v].second][timeJ1+1] + inst.q_v[ordV[v].second]*x[ordV[v].second][a]);	
						//~ cout << "Travel at capacity " << v << " arc " << inst.arcs[ordV[v].second][a] << endl;														
						#endif
						break;
					case 2:
						#ifndef NTravelEmpty
						travelEmpty[ordV[v].second][a].setExpr(inst.q_v[ordV[v].second]*x[ordV[v].second][a] + sV[ordV[v].second][timeJ1+1]);						
						#endif					
						break;
					case 3:
						#ifndef NTravelAtCapacity
						travelAtCapacity[ordV[v].second][a].setExpr(-sV[ordV[v].second][timeJ1+1] + inst.q_v[ordV[v].second]*x[ordV[v].second][a]);																	
						#endif					
						break;
					case 4:
						#ifndef NTravelEmpty
						travelEmpty[ordV[v].second][a].setExpr(inst.q_v[ordV[v].second]*x[ordV[v].second][a] + sV[ordV[v].second][timeJ1+1]);						
						#endif
						break;
					default:
						break;
				}
			}else{
				travelAtCapacity[ordV[v].second][a] = IloRange(env, 0, expr, 0, ss.str().c_str());				
				travelEmpty[ordV[v].second][a] = IloRange(env, 0, expr, 0, ss1.str().c_str());				
			}
		}		
		model.add(travelAtCapacity[ordV[v].second]);
		model.add(travelEmpty[ordV[v].second]);
	}
	
	//Cumulative amount of product that can be purchased from or sold to a spot market by each port - Unlimited
	#ifndef NSpotMarket
	cumSlack = IloRangeArray(env, J);
	for(j=0;j<J;j++){
		stringstream ss;
		ss << "cumSlack_" << j;
		expr.clear();
		for(t=0;t<T;t++){			
			expr += alpha[j][t];							
		}
		cumSlack[j] = IloRange(env, expr, inst.alp_max_j[j], ss.str().c_str());
		//~ cumSlack[j] = IloRange(env, expr, IloInfinity, ss.str().c_str());
	}
	model.add(cumSlack);
	#endif
	
	//Amount loaded/discharged must be in the pre-specified interval 
	operationLowerLimit = IloArray<IloArray<IloRangeArray> > (env, V);
	operationUpperLimit = IloArray<IloArray<IloRangeArray> > (env, V);
	for(v=0;v<V;v++){
		operationLowerLimit[ordV[v].second] = IloArray<IloRangeArray>(env, J);
		operationUpperLimit[ordV[v].second] = IloArray<IloRangeArray>(env, J);
		for(j=0;j<J;j++){
			operationLowerLimit[ordV[v].second][j] = IloRangeArray(env, T);
			operationUpperLimit[ordV[v].second][j] = IloRangeArray(env, T);
			for(t=0;t<T;t++){
				expr.clear();
				stringstream ss, ss1;
				ss << "fjmin_(" << j << "," << t << ")," << v;
				ss1 << "fjmax_(" << j << "," << t << ")," << v;
				#ifndef UselessVariables
				if (inst.inArcs[ordV[v].second][j][t].getSize() > 0){ //Only if there is a entering arc in the node
				#endif
					if (v < vesselsOnModel){
						IloNum minfMax = min(inst.f_max_jt[j][t], inst.q_v[ordV[v].second]);
						operationLowerLimit[ordV[v].second][j][t] = IloRange(env, 0, f[ordV[v].second][j][t] - inst.f_min_jt[j][t]*z[ordV[v].second][j][t] , IloInfinity, ss.str().c_str());				
						operationUpperLimit[ordV[v].second][j][t] = IloRange(env, -IloInfinity, f[ordV[v].second][j][t] - minfMax*z[ordV[v].second][j][t],	0, ss1.str().c_str());				
					}else{
						operationLowerLimit[ordV[v].second][j][t] = IloRange(env, 0, expr, 0, ss.str().c_str());				
						operationUpperLimit[ordV[v].second][j][t] = IloRange(env, 0, expr, 0, ss1.str().c_str());			
					}	
				#ifndef UselessVariables
				}else{
					operationLowerLimit[ordV[v].second][j][t] = IloRange(env, 0, f[ordV[v].second][j][t], 0, ss.str().c_str());				
					operationUpperLimit[ordV[v].second][j][t] = IloRange(env, 0, 0, ss1.str().c_str());							
				}
				#endif
				
			}
			model.add(operationLowerLimit[ordV[v].second][j]);
			model.add(operationUpperLimit[ordV[v].second][j]);
		}		
	}	
	
	#ifndef NOperateAndGo
	//Forces a vessel go to another region after operates if a vessel can load(unload) at a port just in 1 time period
	operateAndGo = IloArray<IloArray<IloRangeArray> >(env, V);
	for (v=0;v<V;v++){ 
		operateAndGo[ordV[v].second] = IloArray<IloRangeArray>  (env, J);
		for(j=0;j<J;j++){						
			operateAndGo[ordV[v].second][j] = IloRangeArray(env,T);
			for(t=0;t<T;t++){		
				stringstream ss;
				ss << "operateAndGo_"<<ordV[v].second<<"_"<<j<<","<<t;
				expr.clear();						
				if (v < vesselsOnModel){
					if (inst.q_v[ordV[v].second] < inst.f_max_jt[j][t]){ //Only if a vessel can load(unload) in a port in just 1 time period						
						for(a=0;a<inst.outRegionArcs[ordV[v].second][j][t].getSize();a++){
							int idA = inst.outRegionArcs[ordV[v].second][j][t][a];						
								expr += x[ordV[v].second][idA];							
						}
						expr += -z[ordV[v].second][j][t];						
					}				
				}
				operateAndGo[ordV[v].second][j][t] = IloRange(env, 0, expr, 0,ss.str().c_str());
			}			
			model.add(operateAndGo[ordV[v].second][j]);
		}
	}
	#endif
}

void Model::unFixVessel(Instance inst, const int& v){
	//Unfixing z
	#ifndef NFixZvar
	for(int j=0; j<inst.numTotalPorts; j++){							
		for(int t=0; t<inst.t; t++){											
			z[v][j][t].setBounds(0, 1);								
		}
	}
	#endif
	//Unfixing X variables
	for(int a=0; a<x[v].getSize(); a++){		
		x[v][a].setBounds(0, 1);
   	}
   	#ifndef NOperateAndGo2
	for(int t=0; t<inst.t; t++){
		w[v][t].setBounds(0, 1);
	}
	#endif
}

void Model::fixVessel(IloEnv env, Instance inst, const int& v){
	if((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){		
		//Get the values		
		xValue[v] = IloNumArray(env, xValue.getSize());			
		cplex.getValues(x[v], xValue[v]);
		#ifndef NOperateAndGo2
		wValue[v] = IloNumArray(env,wValue.getSize());
		cplex.getValues(w[v],wValue[v]);
		#endif
		#ifndef NFixZvar
		for(int j=0; j<inst.numTotalPorts; j++){				
			zValue[v][j] = IloNumArray(env, zValue[v][j].getSize());				
			cplex.getValues(z[v][j], zValue[v][j]);			
		}					
		//Fixing z
		for(int j=0; j<inst.numTotalPorts; j++){							
			for(int t=0; t<inst.t; t++){											
				z[v][j][t].setBounds(round(zValue[v][j][t]),round(zValue[v][j][t]));
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
	}else{//Fix to the previous solution (must be already stored in vectos X and Zvalue)
		cout << "Fixing to previous vessel solution (if it is avaliable)" << endl;
		//~ exit(1);
		//Fixing z
		for(int j=0; j<inst.numTotalPorts; j++){							
			for(int t=0; t<inst.t; t++){											
				z[v][j][t].setBounds(round(zValue[v][j][t]),round(zValue[v][j][t]));
			}
		}
		//Fixing X variables
		for(int a=0; a<x[v].getSize(); a++){		
			x[v][a].setBounds(round(xValue[v][a]), round(xValue[v][a]));
		}		
	}	
}

void Model::fixVesselLessInterval(IloEnv env, Instance inst, const int& v, const int& tS, const int& tF){
	if((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){		
		getSolVals(env,inst);	
		//Fixing z
		for(int j=0; j<inst.numTotalPorts; j++){							
			for(int t=0; t<inst.t; t++){											
				if ( (t < tS) || (t > tF) )
					z[v][j][t].setBounds(round(zValue[v][j][t]),round(zValue[v][j][t]));
			}
		}
		//Fixing X variables
		for(int a=0; a<x[v].getSize(); a++){		
			int t1,t2, arcType;		
			arcType = inst.travelArcType(inst.arcs[v][a],t1,t2);							
			if ( ((arcType == 3 || arcType == 4) && (t1 < tS || t1 > tF) ) || //If is a sink ark out of 'free' interval
			(t2 < tS || t1 > tF) ){  // Travel times 'free' of interval unfixed - intersection remain unfixed
				x[v][a].setBounds(round(xValue[v][a]), round(xValue[v][a]));
		   }			
		}
	}else{//Fix to the previous solution (must be already stored in vectos X and Zvalue)
		cout << "No solution available " << endl;
		exit(1);
		//Fixing
	}	
}

void Model::fixAndOptTW(IloEnv& env, Instance inst, const double& mIntervals, const double& timePerIter, const double& gap, const double& overlap, Timer<chrono::milliseconds>& timer_cplex,float& opt_time, 
const double& timeLimit, float& elapsed_time, double& incumbent){
	double objValue = incumbent;
	double prevObj = 1.0e+10;	
	double objValue1 = incumbent;
	double prevObj1 = 1.0e+10;	
	int i,v,V;	
	Timer<chrono::milliseconds> timer_LS;
	timer_LS.start();
	
	V = inst.speed.getSize();
	cplex.setParam(IloCplex::TiLim, timePerIter);
	if (gap > 1e-04)
		cplex.setParam(IloCplex::EpGap, gap/100);
	int count=0, tS, tF;	
	while((prevObj - objValue > 0.1) && count <= mIntervals*V  && elapsed_time/1000 <= timeLimit){		
		for(i=1;i<=ceil(mIntervals);i++){			
			if(i==1)
				tS = inst.t/mIntervals*(i-1);							
			else
				tS = inst.t/mIntervals*(i-1)*(1-overlap/100);						
			tF = min(tS+inst.t/mIntervals, (double)inst.t);
			cout << "Unfixing " << tS << "..." << tF << endl;
			unFixInterval(inst, tS, tF);
			for (v=0;v<V;v++){
				cout << "Unfixing vessel " << v << " of " << V-1 << endl;
				unFixVessel(inst,v);
				cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
				timer_cplex.start();
				cplex.solve();				
				opt_time += timer_cplex.total();
				elapsed_time += timer_LS.total();
				timer_LS.start();
				prevObj1 = objValue1;
				objValue1 = cplex.getObjValue();
				cout << "Objective: " << objValue1 << endl;				
				if(prevObj1 - objValue1 > 0.0001)
					count = 0;
				else
					//~ count++;				
					
				cout << "Fixing vessel " << v << " of " << V-1 << " minus interval " << tS << "-" << tF << endl;				
				fixVesselLessInterval(env,inst, v, tS, tF);
				if(count >= mIntervals*V || (elapsed_time/1000 >= timeLimit))
					break;
			}			
			//Solve with only interval i free			
			//~ cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
			//~ timer_cplex.start();
			//~ cplex.solve();
			//~ opt_time += timer_cplex.total();
			//~ elapsed_time += timer_LS.total();
			//~ timer_LS.start();
			//~ prevObj1 = objValue1;
			//~ objValue1 = cplex.getObjValue();
			//~ if(prevObj1 - objValue1 > 0.0001)
				//~ count = 0;
			//~ else
				//~ count++;
			//~ cout << "Objective: " << objValue1 << endl;			
			
			//~ if (i < mIntervals){
				cout << "Re-fixing " << tS << "..." << tF;
				fixSolution(env, inst, tS, tF,1);	 //Does not fix the last interval either the sink arcs
				cout << ". Done! " << endl;
			//~ }
			if(count >= mIntervals*V || elapsed_time/1000 >= timeLimit){
				break;
				cout << "Stopped by number of no imrpoved solutions or time" << endl;
			}			
		}
		prevObj = objValue;
		objValue = objValue1;//cplex.getObjValue();	
		if (count >= mIntervals*V+mIntervals || elapsed_time/1000 >= timeLimit){
			cout << "STOP BY V or TIME \n";
			break;
		}
	}
	incumbent = objValue1;
}

void Model::fixVesselPair(IloEnv env, Instance inst, const int& v,const int& v1){
	if((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){		
		//Get the values		
		xValue[v] = IloNumArray(env, xValue.getSize());			
		xValue[v1] = IloNumArray(env, xValue.getSize());			
		cplex.getValues(x[v], xValue[v]);
		cplex.getValues(x[v1], xValue[v1]);
		
		for(int j=0; j<inst.numTotalPorts; j++){				
			zValue[v][j] = IloNumArray(env, zValue[v][j].getSize());				
			cplex.getValues(z[v][j], zValue[v][j]);			
			zValue[v1][j] = IloNumArray(env, zValue[v1][j].getSize());				
			cplex.getValues(z[v1][j], zValue[v1][j]);			
		}					
		//Fixing z
		for(int j=0; j<inst.numTotalPorts; j++){							
			for(int t=0; t<inst.t; t++){											
				z[v][j][t].setBounds(round(zValue[v][j][t]),round(zValue[v][j][t]));
				z[v1][j][t].setBounds(round(zValue[v1][j][t]),round(zValue[v1][j][t]));
			}
		}
		//Fixing X variables
		for(int a=0; a<x[v].getSize(); a++){		
			x[v][a].setBounds(round(xValue[v][a]), round(xValue[v][a]));
		}
		for(int a=0; a<x[v1].getSize(); a++){		
			x[v1][a].setBounds(round(xValue[v1][a]), round(xValue[v1][a]));
		}
	}else{//Fix to the previous solution (must be already stored in vectos X and Zvalue)
		cout << "Erro! Model is not feasible or was modified after solved" << endl;
		exit(1);		
	}	
}

bool pairCompare(const std::pair<int,double> &left, const std::pair<int,double> &right) {
    return left.second < right.second;
}

/*
 * Fix the overall solution, and for each vessel v \in V, unfix it variables and solve the model
 * Pre-req: call of method fixAllSolution() 
 */ 
void Model::fixAndOptmizeH(IloEnv& env, Instance inst, const double& timePerIter, const double& gap, double& incumbent, Timer<chrono::milliseconds>& timer_cplex,float& opt_time,
const double& timeLimit, float& elapsed_time){
	Timer<chrono::milliseconds> timer_LS;
	timer_LS.start();	
	int v, v1, V = inst.speed.getSize();				
	if (gap > 1e-04)
		cplex.setParam(IloCplex::EpGap, gap/100);
	cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
	double currentObj = incumbent;
	double previousObj = 1e+20;
	float elapsed_local_time{0};
	#ifdef NPairsVessels
	//Unfixing one vesseal per iteration
	while ((previousObj - currentObj > 0.0001) && elapsed_time < timeLimit){
		for (v=0;v<V;v++){
			cout << "Unfixing vessel " << v << " of " << V << endl;
			unFixVessel(inst,v);
			cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
			timer_cplex.start();
			cplex.solve();
			opt_time += timer_cplex.total();
			elapsed_time += timer_LS.total();
			if(elapsed_time/1000 >= timeLimit){				
				break;
			}			
			timer_LS.start();
			
			cout << "Objective: " << cplex.getObjValue() << endl;
			fixVessel(env,inst, v);
			cout << "Fixing vessel " << v << " of " << V << endl;
		}
		if(elapsed_time/1000 >= timeLimit){
			cout << "Reached LS time limit " << timeLimit << ": " << elapsed_time/1000 << endl;
			break;
		}
		timer_cplex.start();
		cplex.solve();
		opt_time += timer_cplex.total();
		previousObj = currentObj;
		currentObj = cplex.getObjValue();
	}	
	#endif
	
	#ifndef NPairsVessels
	//Get random pairs of vessels
	int count, rId, v2, sumComb = 0;  
	double prevObj;
	vector <pair<unsigned int,unsigned int> > vesselsComb;
	srand(V);

	//Local search
	while ((previousObj - currentObj > 0.0001) && elapsed_local_time/1000 < 
	timeLimit  && count <= V){
		//Make pairs
		sumComb = 0;
		for(v=0;v<V-1;v++){
			for(v1=v+1;v1<V;v1++){			
				vesselsComb.push_back(make_pair(v,v1));
				sumComb++;
			}
		}		
		for (int i=0;i<sumComb;i++){
			cout << "Size: " << vesselsComb.size() << endl;
			rId = iRand(0,vesselsComb.size()-1);			
			v1 = vesselsComb[rId].first;
			v2 = vesselsComb[rId].second;
			vesselsComb.erase(vesselsComb.begin()+rId);
			cout << "Unfixing vessel " << v1 << " and " << v2 << endl;
			unFixVessel(inst,v1);
			unFixVessel(inst,v2);			
			cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_local_time/1000,0.0)));
			timer_cplex.start();
			cplex.solve();
			opt_time += timer_cplex.total();
			prevObj = incumbent;
			incumbent = cplex.getObjValue();
			if(prevObj - incumbent > 0.0001)
					count = 0;
				else
					count++;				
			cout << "Objective " << incumbent << endl;			
			fixVesselPair(env, inst, v1,v2); //Fix the last 2 vessel
			elapsed_local_time += timer_LS.total();
			if(elapsed_local_time/1000 >= timeLimit || count > V){
				break;
			}
			timer_LS.start();			
		}						
		previousObj = currentObj;
		currentObj = incumbent;				
	}	
	elapsed_time += elapsed_local_time;
	#endif	
}

void Model::addVesselToModel(IloEnv& env, Instance inst, const int& vAdd){
	int j,t,a;	
	int J = inst.numTotalPorts;
	int T = inst.t;
	int V = inst.speed.getSize(); //# of vessels	
	unordered_map<string,double>::const_iterator it;

	///Objective function - Just add variables of index vAdd
	IloExpr current = cplex.getObjective().getExpr(); //Gets the current objective function
	IloExpr expr(env);
	IloExpr expr1(env);	
	for(j=0;j<J;j++){
		for(t=0;t<T;t++){
			#ifndef UselessVariables				
			if (inst.inArcs[vAdd][j][t].getSize() > 0){ 
			#endif
				expr += inst.r_jt[j][t]*f[vAdd][j][t];				
				expr1 += (t*inst.perPeriodRewardForFinishingEarly)*z[vAdd][j][t];
			#ifndef UselessVariables
			}
			#endif
		}		
	}		
	for (a=0;a<inst.arcs[vAdd].size();a++){
		it = inst.c_va[vAdd].find(inst.arcs[vAdd][a]);				
		if (it == inst.c_va[vAdd].end()){
			cout << "Not found key for arc " << inst.arcs[vAdd][a] << endl;
			exit(1);
		}				
		expr1 +=  it->second * x[vAdd][a];
	}	
	obj.setExpr(current + expr1-expr);		
	
	///Constraints
	//Flow balance in port-time nodes, source and sink		
	for(j=0;j<J+1;j++){ //J+1 used for access the information of source and sink node
		if(j<J){ //If is a real port			
			for(t=0;t<T;t++){	
				expr.clear();
				#ifndef UselessVariables
				if (inst.inArcs[vAdd][j][t].getSize() > 0){ 
				#endif											
					for (a=0;a<inst.outArcs[vAdd][j][t].getSize();a++){						
						int idA = inst.outArcs[vAdd][j][t][a];							
						expr += x[vAdd][idA];						
					}
					for (a=0;a<inst.inArcs[vAdd][j][t].getSize();a++){
						int idA = inst.inArcs[vAdd][j][t][a];													
						expr += -x[vAdd][idA];		
					}							
				#ifndef UselessVariables
				}
				#endif								
				flowBalance[vAdd][j][t].setExpr(expr);				
			}			
		}else{			
			expr.clear();
			expr1.clear();						
			for (a=0;a<inst.outArcs[vAdd][J][0].getSize();a++){
				//summing for source node
				int idA = inst.outArcs[vAdd][J][0][a];										
				expr += x[vAdd][idA];				
			}
			//~ if (vAdd == 5) cout << "Expr sourceArc5 " << expr << endl;
			for (a=0;a<inst.inArcs[vAdd][J][0].getSize();a++){
				//summing for sink node
				int idA = inst.inArcs[vAdd][J][0][a];										
				expr1 += x[vAdd][idA];		
			}

			flowBalance[vAdd][j][0].setBounds(1,1);
			flowBalance[vAdd][j][1].setBounds(1,1);

			//Flow Balance source node
			flowBalance[vAdd][j][0].setExpr(expr);			
			//Flow Balance sink node
			flowBalance[vAdd][j][1].setExpr(expr1);
		}
	}
	
	//Inventory balance ports 
	for(j=0;j<J;j++){
		for(t=0;t<T;t++){ 
			expr.clear();
			expr = portInvBalance[j][t].getExpr();						
			#ifndef UselessVariables
			if (inst.inArcs[vAdd][j][t].getSize() > 0) //Only if there is a intering arc in the node
			#endif
				expr += (-inst.delta[j] * -f[vAdd][j][t]);			
			portInvBalance[j][t].setExpr(expr);			
		}		
	}
	
	//Vessels inventory balance	
	for(t=0;t<T;t++){		
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
		for(t=0;t<T;t++){
			expr.clear();
			expr = berthLimit[j][t].getExpr();
			#ifndef UselessVariables
			if (inst.inArcs[vAdd][j][t].getSize() > 0) //Only if there is a intering arc in the node
			#endif
				expr += z[vAdd][j][t];				
			
			berthLimit[j][t].setExpr(expr);
		}		
	}
	#endif
	
	//Vessels only attempt to load/discharge at node only if the vessel is actually at the node
	for(j=0;j<J;j++){		
		for(t=0;t<T;t++){			
			expr.clear();			
			#ifndef UselessVariables
			if (inst.inArcs[vAdd][j][t].getSize() > 0){ 
			#endif				
				for(a=0;a<inst.inArcs[vAdd][j][t].getSize();a++){
					int idA = inst.inArcs[vAdd][j][t][a];
					expr += x[vAdd][idA];
				}					
				atemptToOperate[j][t][vAdd].setBounds(-IloInfinity, 0);
				atemptToOperate[j][t][vAdd].setExpr(z[vAdd][j][t] - expr);				
			#ifndef UselessVariables
			}else{
				expr.clear();
				//~ atemptToOperate[j][t][vAdd] = IloRange(env, 0, z[vAdd][j][t], 0, ss.str().c_str());				
			}
			#endif			
		}
	}
	
	//Vessels must travel at capacity and empty 	
	expr.clear();	
	for(a=0;a<inst.arcs[vAdd].size();a++){			
		int timeJ1, timeJ2;			
		int arcType = inst.travelArcType(inst.arcs[vAdd][a], timeJ1, timeJ2);					
		travelAtCapacity[vAdd][a].setBounds(-IloInfinity, 0);
		travelEmpty[vAdd][a].setBounds(-IloInfinity, inst.q_v[vAdd]);
		switch (arcType)
		{
			case 1:
				#ifndef NTravelAtCapacity
				 //timeJ1 is incresead in 1 because in vector sV time 0 is considered as source node										
				travelAtCapacity[vAdd][a].setExpr(-sV[vAdd][timeJ1+1] + inst.q_v[vAdd]*x[vAdd][a]);	
				//~ cout << "Travel at capacity " << v << " arc " << inst.arcs[v][a] << endl;														
				#endif
				break;
			case 2:
				#ifndef NTravelEmpty
				travelEmpty[vAdd][a].setExpr(inst.q_v[vAdd]*x[vAdd][a] + sV[vAdd][timeJ1+1]);						
				#endif					
				break;
			case 3:
				#ifndef NTravelAtCapacity
				travelAtCapacity[vAdd][a].setExpr(-sV[vAdd][timeJ1+1] + inst.q_v[vAdd]*x[vAdd][a]);																	
				#endif					
				break;
			case 4:
				#ifndef NTravelEmpty
				travelEmpty[vAdd][a].setExpr(inst.q_v[vAdd]*x[vAdd][a] + sV[vAdd][timeJ1+1]);						
				#endif
				break;
			default:
				break;
		}		
	}
	
	//Cumulative amount of product that can be purchased from or sold to a spot market by each port - Reduce UB only in the last vessel
	#ifndef NSpotMarket
	if (vAdd == V-1){
		//~ for(j=0;j<J;j++){
			//~ cumSlack[j].setUB(inst.alp_max_j[j]);
			//~ for(t=0;t<T;t++){
				//~ alpha[j][t].setBounds(0,inst.alp_max_jt[j][t]);
			//~ }
		//~ }
	}	
	#endif
	
	//Amount loaded/discharged must be in the pre-specified interval 		
	for(j=0;j<J;j++){		
		for(t=0;t<T;t++){
			expr.clear();			
			#ifndef UselessVariables
			if (inst.inArcs[vAdd][j][t].getSize() > 0){ //Only if there is a entering arc in the node
			#endif				
				IloNum minfMax = min(inst.f_max_jt[j][t], inst.q_v[vAdd]);
				operationLowerLimit[vAdd][j][t].setBounds(0,IloInfinity);
				operationLowerLimit[vAdd][j][t].setExpr(f[vAdd][j][t] - inst.f_min_jt[j][t]*z[vAdd][j][t]);
				
				operationUpperLimit[vAdd][j][t].setBounds(-IloInfinity, 0);
				operationUpperLimit[vAdd][j][t].setExpr(f[vAdd][j][t] - minfMax*z[vAdd][j][t]);					
			#ifndef UselessVariables
			}
			#endif			
		}		
	}		
	
	#ifndef NOperateAndGo
	//Forces a vessel go for another region after operates (and if have no enough capacity to operate)		
	for(j=0;j<J;j++){									
		for(t=0;t<T;t++){
			expr.clear();				
			if (inst.q_v[vAdd] < inst.f_max_jt[j][t]){ //Only if a vessel can load(unload) in a port in just 1 time period										
				for(a=0;a<inst.outRegionArcs[vAdd][j][t].getSize();a++){
					int idA = inst.outRegionArcs[vAdd][j][t][a];
					expr += x[vAdd][idA];
				}				
				expr += -z[vAdd][j][t];
			}
			operateAndGo[vAdd][j][t].setExpr(expr);
		}		
	}	
	#endif	
}

void Model::integralizeVessel(const Instance& inst, const int& vInt){
	int J = inst.numTotalPorts;
	int T = inst.t;
	int a,j,t;
	///Remove conversion for vessel vInt	
	for(a=0; a<x[vInt].getSize(); a++){		
		model.remove(convertX[vInt][a]);		
	}	
	for(j=0;j<J;j++){	
		for(t=0;t<T;t++){			
			model.remove(convertZ[vInt][j][t]);
		}	
	}

}

void Model::resetObjFunction(IloEnv& env, Instance inst){
	int v,j,t,a;
	int V = inst.speed.getSize();
	int J = inst.numTotalPorts;
	int T = inst.t;
	unordered_map<string,double>::const_iterator it;
	IloExpr expr(env);
	IloExpr expr1(env);
	for(v=0;v<V;v++){
		for(j=0;j<J;j++){
			for(t=0;t<T;t++){
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
			it = inst.c_va[v].find(inst.arcs[v][a]);				
			if (it == inst.c_va[v].end()){
				cout << "Not found key for arc " << inst.arcs[v][a] << endl;
				exit(1);
			}				
			expr1 +=  it->second * x[v][a];									
		}
	}	
	for(j=0;j<J;j++){
		for(t=0;t<T;t++){			
			#ifndef NSpotMarket
			expr1 += inst.p_jt[j][t]*alpha[j][t];		
			#endif
		}		
	}	
	obj.setExpr(expr1-expr);	
}

void mirp::fixAndRelaxH(string file, string optStr, const double& gapFirst, const int& outVessels, const int& timeLimitFirst, 
	const double& mIntervals, const int& timeLimitSecond, const double& gapSecond, const double& overlap2){
	///Time parameters
	Timer<chrono::milliseconds> timer_cplex;
	Timer<chrono::milliseconds> timer_global;
	timer_global.start();
	float global_time {0};
	float opt_time {0};
		
	IloEnv env;	
	try{
		///Read input files
		Instance inst(env);
		inst.readInstance(env, file);	
			
		int v, addVesselS, addVesselF, vInt;
		int J = inst.numTotalPorts;
		int T = inst.t;
		int V = inst.speed.getSize(); //# of vessels
		double obj1stPhase = 0, obj2ndPhase = 0, time1stPhase=0, time2ndPhase=0;
		
		//Verify if the number of outVessels is too big
		if (V - outVessels < 2){
			cout << "Error: number of vessels out of model should not exced V - 2" << endl;
			exit(1);
		}	
		/// NEW MODEL
		Model model(env);
		model.buildFixAndRelaxHorizontalModel(env,inst,outVessels);
		model.setParameters(env, timeLimitFirst, gapFirst);	
		
		/// 1st phase
		for (v=0; v<V; v++){
			timer_cplex.start();		
			model.cplex.solve();			
			opt_time += timer_cplex.total();		
			cout << "Solution Status: " << model.cplex.getStatus() << endl;		
			cout << "Objective: " << model.cplex.getObjValue() << endl;
				
			cout << "Fixing vessel " << model.ordV[v].second << "...";
			model.fixVessel(env, inst, model.ordV[v].second);
			cout << "Fixed! \n";
		
			addVesselS = V - outVessels + v;							
			if(addVesselS < V){
				cout << "Adding vessel " << model.ordV[addVesselS].second << " to the model." << endl;						
				model.addVesselToModel(env,inst,model.ordV[addVesselS].second);
			}
			vInt = v + 1;
			if(vInt < V){		
				cout << "Integralizing vessel " << model.ordV[vInt].second <<  endl;			
				model.integralizeVessel(inst, model.ordV[vInt].second);		
			}				
			double newGap = (gapFirst - (gapFirst/V)*(v+1)) / 100;
			cout << "New GAP " << newGap << endl;
			model.cplex.setParam(IloCplex::EpGap, newGap);	
		}
		//Last iteration (after fixing the penultimate vessel)
		timer_cplex.start();
		model.cplex.solve();
		opt_time += timer_cplex.total();		
		time1stPhase = timer_global.total();
		global_time += time1stPhase;
		timer_global.start();
		obj1stPhase = model.cplex.getObjValue();
			
		///2nd phase
		cout << "Improving solution..." << endl;				
		
		//~ model.regionLocalSearch(env, inst,timeLimitSecond, gapSecond*0.5, timer_cplex, opt_time);
		
		//~ model.fixAllSolution(env, inst); 	
		//~ model.fixAndOptmizeH(env, inst, timeLimitSecond, gapSecond, obj1stPhase,timer_cplex,opt_time);	

		//~ model.improvementPhase(env, inst, mIntervals, timeLimitSecond, gapSecond, overlap2, timer_cplex, opt_time);
							
		global_time += timer_global.total();			
		time2ndPhase = global_time - time1stPhase;
		obj2ndPhase	= model.cplex.getObjValue();
		
		model.printSolution(env, inst);
		cout << "GAP Fisrt iter: " << gapFirst << " |outVessels| " << outVessels << " Time per vessel " << timeLimitFirst << endl
		<< "GAP Second phase " << gapSecond << " Time Limit second phase " << timeLimitSecond  << endl
		<< "CPLEX time: " << opt_time/1000 << endl << "Other times: " << (global_time-opt_time)/1000 << endl
		<< "Total time: " << global_time/1000 << endl
		<< "1st phase time: " << time1stPhase/1000 << endl
		<< "2nd phase time: " << time2ndPhase/1000 << endl
		<< "1st phase obj: " << obj1stPhase << endl
		<< "2nt phase obj: " << obj2ndPhase << endl
		<< "Improvment from 1st to 2nd phase: " << abs((obj2ndPhase/obj1stPhase - 1)*100) << "%" << endl;
	}catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;		
		e.end();
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}
	env.end();
}
