/*Fixed charge network flow for a single-product MIRP
 * Index i and j of ports J starts in 1
 * 0 is the source node, and J+1 is the sink node
 * Index t for time period start in 1
 * t=0 is the time in which vessesl starts its voyage from source node
 * Index v for vessels starts in 0
 *  
 */ 
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <ilcplex/ilocplex.h>
#include "../include/util.h"
#include "../include/model.h"
//~ #define NDEBUG
#include <assert.h>
#define UselessVariables
//~ #define NSpotMarket
//~ #define NTravelAtCapacity
//~ #define NTravelEmpty
//~ #define NBerthLimit
#define NBranching
#define NValidInequalities
#define NOperateAndGo
#define N2PortNorevisit
#define NRelaxation

ILOSTLBEGIN

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloNumArray> NumNumMatrix;

using namespace std;
using mirp::Model;

void Model::buildModel(IloEnv& env, Instance inst){
	int i,j,t,v,a;
	int J = inst.numTotalPorts;
	int T = inst.t + 1;			//Init time-periods in 1
	int V = inst.speed.getSize(); //# of vessels
	int N = J + 2; // Number of nodes including source and sink. (0,1...j,j+1)
		
	//Init Variables and info arrays
	#ifndef NSpotMarket
	alpha = NumVarMatrix(env, N);
	#endif
	sP = NumVarMatrix(env,N);
	f = IloArray<NumVarMatrix> (env, V);
	fX = IloArray<IloArray<NumVarMatrix> >(env, V);
	fOA = IloArray<NumVarMatrix> (env, V);
	fOB = IloArray<NumVarMatrix> (env, V);
	fW = IloArray<NumVarMatrix> (env, V);
	
	x = IloArray<IloArray<IloArray<IloBoolVarArray> > >(env, V);
	hasArc = IloArray<IloArray<IloArray<IloIntArray> > >(env, V);
	arcCost = IloArray<IloArray<IloArray<IloNumArray> > >(env, V);
	hasEnteringArc1st = IloArray<IloArray<IloIntArray> >(env, V);
	z = IloArray<IloArray<IloBoolVarArray> >(env,V);
	oA = IloArray<IloArray<IloBoolVarArray> >(env,V);
	oB = IloArray<IloArray<IloBoolVarArray> >(env,V);
	w = IloArray<IloArray<IloBoolVarArray> >(env,V);
	
	//Additional variables for branching
	#ifndef NBranching		
	#endif
	
	for(v=0;v<V;v++){
		f[v] = NumVarMatrix(env,N);
		fX[v] = IloArray<NumVarMatrix>(env,N);
		fOA[v] = NumVarMatrix(env,N);
		fOB[v] = NumVarMatrix(env,N);
		fW[v] = NumVarMatrix(env,N);
		
		z[v] = IloArray<IloBoolVarArray>(env,N);
		oA[v] = IloArray<IloBoolVarArray>(env,N);
		oB[v] = IloArray<IloBoolVarArray>(env,N);
		w[v] = IloArray<IloBoolVarArray>(env,N);
		x[v] = IloArray<IloArray<IloBoolVarArray> >(env, N);				
		hasArc[v] = IloArray<IloArray<IloIntArray> >(env, N);				
		arcCost[v] = IloArray<IloArray<IloNumArray> >(env, N);
		hasEnteringArc1st[v] = IloArray<IloIntArray>(env, N);
			
		for(j=0;j<N;j++){			
			x[v][j] = IloArray<IloBoolVarArray>(env,N);
			hasArc[v][j] = IloArray<IloIntArray>(env,N);
			arcCost[v][j] = IloArray<IloNumArray>(env,N);			
			fX[v][j] = IloArray<IloNumVarArray>(env, N);
			for(i=0;i<N;i++){
				x[v][j][i] = IloBoolVarArray(env,T);
				hasArc[v][j][i] = IloIntArray(env,T);
				arcCost[v][j][i] = IloNumArray(env,T);
				fX[v][j][i] = IloNumVarArray(env, T);
				for(t=0;t<T;t++){				
					stringstream ss;			
					ss << "x_"<<v<<","<<j<<","<<i<<","<<t;
					x[v][j][i][t].setName(ss.str().c_str());
					ss.str(string());
					ss << "fX_"<<v<<","<<j<<","<<i<<","<<t;
					fX[v][j][i][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
				}
			}
			if(j>0 && j<=J){ //Only considering Ports
				f[v][j] = IloNumVarArray(env, T);
				fOA[v][j] = IloNumVarArray(env, T);
				fOB[v][j] = IloNumVarArray(env, T);
				fW[v][j] = IloNumVarArray(env, T);
				z[v][j] = IloBoolVarArray(env, T);
				oA[v][j] = IloBoolVarArray(env, T);
				oB[v][j] = IloBoolVarArray(env, T);
				w[v][j] = IloBoolVarArray(env, T);
				hasEnteringArc1st[v][j] = IloIntArray(env,T);
				for(t=1;t<T;t++){				
					stringstream ss;
					ss << "f_(" << j << "," << t << ")," << v;
					int maxAmountOperation = min(inst.q_v[v], inst.f_max_jt[j-1][t-1]);
					f[v][j][t] = IloNumVar(env, 0, maxAmountOperation, ss.str().c_str());								
					ss.str(string());
					ss << "fOA_(" << j << "," << t << ")," << v;
					fOA[v][j][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
					ss.str(string());
					ss << "z_(" << j << "," << t << ")," << v;
					z[v][j][t].setName(ss.str().c_str());
					ss.str(string());
					ss << "oA_(" << j << "," << t << ")," << v;
					oA[v][j][t].setName(ss.str().c_str());					
					if(t<T-1){ //Variables that haven't last time index
						ss.str(string());
						ss << "fOB_(" << j << "," << t << ")," << v;
						fOB[v][j][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
						ss.str(string());
						ss << "fW_(" << j << "," << t << ")," << v;
						fW[v][j][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
						ss.str(string());
						ss << "oB_(" << j << "," << t << ")," << v;
						oB[v][j][t].setName(ss.str().c_str());
						ss.str(string());
						ss << "w_(" << j << "," << t << ")," << v;
						w[v][j][t].setName(ss.str().c_str());					
					}									
				}		
			}	
		}
	}

	for(i=1;i<N-1;i++){		
		alpha[i] = IloNumVarArray(env,T);		
		sP[i] = IloNumVarArray(env, T);		
		for(t=0;t<T;t++){
			stringstream ss;						
			if (t == 0){								
				ss << "sP_(" << i << ",0)";
				sP[i][t] = IloNumVar(env,inst.s_j0[i-1], inst.s_j0[i-1], ss.str().c_str()); //Initial inventory		
			}else{
				ss.str(string());
				ss << "sP_(" << i << "," << t << ")";				
				sP[i][t] = IloNumVar(env,inst.sMin_jt[i-1][0], inst.sMax_jt[i-1][0], ss.str().c_str()); //As the port capacity is fixed, always used the data from index 0
				
				ss.str(string());
				ss << "alpha_(" << i << "," << t << ")";
				alpha[i][t] = IloNumVar(env, 0, inst.alp_max_jt[i-1][t-1], ss.str().c_str());
			}			
		}		
	}
	
	///Objective function
	IloExpr expr(env);
	IloExpr expr1(env);	
	for(v=0;v<V;v++){
		///Builds the arcs between nodes 												//2nd term
		//Source arc	
		j = inst.initialPort[v]+1;
		hasArc[v][0][j][0] = 1;			//Time between source node and initial ports is equal to first time available
		arcCost[v][0][j][0] = inst.portFee[j-1];
		expr1 += arcCost[v][0][j][0]*x[v][0][j][0];
		hasEnteringArc1st[v][j][inst.firstTimeAv[v]+1] = 1;
		x[v][0][j][0].setBounds(1,1); // Fixing initial port position
				
		//Travel and sink arcs		
		i = inst.initialPort[v]+1;					//Arcs from initial port i
		for (t=inst.firstTimeAv[v]+1;t<T;t++){		//and initial time available t			
			for(j=1;j<N-1;j++){						//Not necessary to account sink node
				if (i != j){
					int t2 = t + inst.travelTime[v][i-1][j-1]; 
					if (t2<T){ 	//If exists time to reach port j 
						double arc_cost;
						if (inst.typePort[i-1]==1 && inst.typePort[j-1]==0){ //If port i is consuming and j is a producer port
							arc_cost = inst.costKm[v]*inst.distanceMatrix[i-1][j-1]*(1-inst.trav_empt[v]) + inst.portFee[j-1];		
						}else{
							arc_cost = inst.costKm[v]*inst.distanceMatrix[i-1][j-1] + inst.portFee[j-1];
						}
						hasArc[v][i][j][t] = 1;	
						//~ cout << v << " " << i << " " << j << " " << t << "=" << hasArc[v][i][j][t] << endl;
						arcCost[v][i][j][t] = arc_cost;
						hasEnteringArc1st[v][j][t2] = 1;
						if (t2+1<T-1) 			//Waiting arc 
							hasEnteringArc1st[v][j][t2+1] = 1;
							
						expr1 += arc_cost*x[v][i][j][t];
												
						//Sink arc from port j 
						hasArc[v][j][N-1][t2] = 1;
						arcCost[v][j][N-1][t2] = -(T-t2-1)*inst.perPeriodRewardForFinishingEarly;
						expr1 += arcCost[v][j][N-1][t2]*x[v][j][N-1][t2];						
						//when the time t reach a time that can be built a mirror between j and i
						if (t >= inst.firstTimeAv[v]+1 + inst.travelTime[v][i-1][j-1]){ //+1 for the data convert 
							if (inst.typePort[j-1]==1 && inst.typePort[i-1]==0){ 		  //If port j is consuming and i is a producer port
								arc_cost = inst.costKm[v]*inst.distanceMatrix[j-1][i-1]*(1-inst.trav_empt[v]) + inst.portFee[i-1];		
							}else{
								arc_cost = inst.costKm[v]*inst.distanceMatrix[j-1][i-1] + inst.portFee[i-1];
							}
							hasArc[v][j][i][t] = 1;
							arcCost[v][j][i][t] = arc_cost;
							hasEnteringArc1st[v][i][t2] = 1;
							expr1 += arc_cost*x[v][j][i][t];							
						}
						//~ //Create arc from j,t2 to others ports (j2) in time t3
						for(int j2=1;j2<=J;j2++){
							if(j2 != i && j2 != j){
								int t3 = t2+inst.travelTime[v][j-1][j2-1];  
								if(t3<T){
									if (inst.typePort[j-1]==1 && inst.typePort[j2-1]==0){
										arc_cost = inst.costKm[v]*inst.distanceMatrix[j-1][j2-1]*(1-inst.trav_empt[v]) + inst.portFee[j2-1];		
									}else{
										arc_cost = inst.costKm[v]*inst.distanceMatrix[j-1][j2-1] + inst.portFee[j2-1];
									}
									hasArc[v][j][j2][t2] = 1;
									arcCost[v][j][j2][t2] = arc_cost;
									expr1 += arc_cost*x[v][j][j2][t2];
									hasEnteringArc1st[v][j2][t3] = 1;
								}
							}
						}						
					}
				}
			}
			if(t+1<T)
				hasEnteringArc1st[v][i][t+1] = 1;
			//Sink arc from port i
			hasArc[v][i][N-1][t] = 1;
			arcCost[v][i][N-1][t] = -(T-t-1)*inst.perPeriodRewardForFinishingEarly;
			expr1 += arcCost[v][i][N-1][t]*x[v][i][N-1][t];			
		}		
		
		
		for(i=1;i<N-1;i++){		//Only considering ports
			for(t=1;t<T;t++){					
				if(hasEnteringArc1st[v][i][t]){
					expr += inst.r_jt[i-1][t-1]*f[v][i][t];									//1st term
					expr1 += (t-1)*inst.attemptCost*z[v][i][t];								//3rd term									   
				}
			}			
		}		
		
	}	
	for(j=1;j<N-1;j++){
		for(t=1;t<T;t++){			
			expr1 += inst.p_jt[j-1][t-1]*alpha[j][t];									//4rd term
		}
	}
	obj.setExpr(expr1-expr);	
	model.add(obj);
	//~ cout << hasEnteringArc1st[0][4][7] << endl;
	///Constraints
	//Flow balance source and sink nodes
	IloExpr expr_sinkFlow(env);
	//~ IloExpr expr_sourceFlow;
	sinkNodeBalance = IloRangeArray(env,V,1,1); 
	sourceNodeBalance = IloRangeArray(env,V,1,1);
	//First level balance
	IloExpr expr_1stLevel(env);
	IloExpr expr_1stFlow(env);
	firstLevelBalance = IloArray<IloArray<IloRangeArray> > (env, V);
	firstLevelFlow = IloArray<IloArray<IloRangeArray> > (env, V);
	
	//Second level balance
	IloExpr expr_2ndLevel(env);
	IloExpr expr_2ndFlow(env);
	secondLevelBalance = IloArray<IloArray<IloRangeArray> > (env, V);
	secondLevelFlow = IloArray<IloArray<IloRangeArray> > (env, V);
	
	//Link betwenn 1st and 2nd level
	linkBalance = IloArray<IloArray<IloRangeArray> > (env, V);
	
	//Berth Limit
	berthLimit = IloArray<IloRangeArray>(env, N-1);		//Index 0 ignored
	IloExpr expr_berth(env);
	
	//Travel at capacity and empty
	travelAtCapacity = IloArray<IloArray<IloArray<IloRangeArray> > > (env, V);
	travelEmpty = IloArray<IloArray<IloArray<IloRangeArray> > > (env, V);
	
	//Inventory balance at ports
	portInventory = IloArray<IloRangeArray> (env, N-1);
	IloExpr expr_invBalancePort(env);
	
	//Cumulative spot market
	cumSlack = IloRangeArray(env,N-1);
	IloExpr expr_cumSlack(env);
	
	//Limits on operation values
	operationLowerLimit = IloArray<IloArray<IloRangeArray> >(env,V);
	operationUpperLimit = IloArray<IloArray<IloRangeArray> >(env,V);
	
	//Flow limits
	flowCapacityX = IloArray<IloArray<IloArray<IloRangeArray> > > (env,V);
	flowCapacityOA = IloArray<IloArray<IloRangeArray> >(env,V);
	flowCapacityOB = IloArray<IloArray<IloRangeArray> >(env,V);
	flowCapacityW = IloArray<IloArray<IloRangeArray> >(env,V);
	
	for(i=1;i<N-1;i++){			
		stringstream ss1;
		berthLimit[i] = IloRangeArray(env,T);
		expr_cumSlack.clear();
		portInventory[i] = IloRangeArray(env,T);				
		for(t=1;t<T;t++){
			expr_berth.clear();
			expr_invBalancePort.clear();
			stringstream ss, ss2;
			ss << "berthLimit_(" << i << "," << t << ")";
			bool emptyExpr = true;
			for(v=0;v<V;v++){				
				if (hasEnteringArc1st[v][i][t]==1){ //Only if exists an entering arc in the node				
					expr_berth += z[v][i][t];
					emptyExpr = false;
				}
				expr_invBalancePort += -f[v][i][t];
			}			
			if(!emptyExpr){ //Only if there are some Z variable in the expr
				berthLimit[i][t] = IloRange(env,-IloInfinity, expr_berth, inst.b_j[i-1], ss.str().c_str());									
				model.add(berthLimit[i][t]);
			}
			
			expr_cumSlack += alpha[i][t];
			expr_invBalancePort += -alpha[i][t];
			ss2 << "invBalancePort_(" << i << "," << t << ")";
			portInventory[i][t] = IloRange(env, inst.delta[i-1]*inst.d_jt[i-1][t-1],
				sP[i][t]-sP[i][t-1]-inst.delta[i-1]*expr_invBalancePort,
				inst.delta[i-1]*inst.d_jt[i-1][t-1], ss2.str().c_str());
			model.add(portInventory[i][t]);
			
		}		
		ss1 << "cum_slack("<<i<<")";
		cumSlack[i] = IloRange(env, expr_cumSlack, inst.alp_max_j[i-1], ss1.str().c_str());
		model.add(cumSlack[i]);
		
	}	
	expr_cumSlack.end();
	expr_berth.end();
	expr_invBalancePort.end();
	
	for(v=0;v<V;v++){
		expr_sinkFlow.clear();		
						
		firstLevelBalance[v] = IloArray<IloRangeArray>(env,J+1); //Only for ports - id 0 not used		
		firstLevelFlow[v] = IloArray<IloRangeArray>(env,J+1); 
		secondLevelBalance[v] = IloArray<IloRangeArray>(env,J+1); 
		secondLevelFlow[v] = IloArray<IloRangeArray>(env,J+1); 
		linkBalance[v] = IloArray<IloRangeArray>(env,J+1); 
		
		travelAtCapacity[v] = IloArray<IloArray<IloRangeArray> > (env, N-1);
		travelEmpty[v] = IloArray<IloArray<IloRangeArray> > (env, N-1);
		
		operationLowerLimit[v] = IloArray<IloRangeArray> (env,N-1);
		operationUpperLimit[v] = IloArray<IloRangeArray> (env,N-1);
		
		flowCapacityX[v] = IloArray<IloArray<IloRangeArray> >(env,N);
		flowCapacityOA[v] = IloArray<IloRangeArray> (env,N-1);
		flowCapacityOB[v] = IloArray<IloRangeArray> (env,N-1);
		flowCapacityW[v] = IloArray<IloRangeArray> (env,N-1);
		
		for(i=0;i<N;i++){
			if(i>0 && i <= J){ //Only considering ports				//ERROR index superior inside 
				firstLevelBalance[v][i] = IloRangeArray(env,T,0,0); //Id 0 is not used				
				secondLevelBalance[v][i] = IloRangeArray(env,T,0,0); 
				firstLevelFlow[v][i] = IloRangeArray(env,T,0,0); 
				secondLevelFlow[v][i] = IloRangeArray(env,T,0,0); 
				linkBalance[v][i] = IloRangeArray(env,T,0,0); 
				
				travelAtCapacity[v][i] = IloArray<IloRangeArray> (env, N);
				travelEmpty[v][i] = IloArray<IloRangeArray> (env, N);
				
				operationLowerLimit[v][i] = IloRangeArray (env,T);
				operationUpperLimit[v][i] = IloRangeArray (env,T);
				
				for(j=0;j<N;j++){
					if(j>0){ //Considering ports and sink node
						travelAtCapacity[v][i][j] = IloRangeArray(env, T, -IloInfinity, 0);
						travelEmpty[v][i][j] = IloRangeArray (env, T, -IloInfinity, inst.q_v[v]);
						for(t=1;t<T;t++){
							stringstream ss, ss1;
							if (j<=J){				//When j is a port
								if(inst.typePort[i-1] == 0 && inst.typePort[j-1] == 1){
									ss << "travelAtCap_(" << i << "," << j << ")_" << t << "," << v;
									travelAtCapacity[v][i][j][t].setExpr(-fX[v][i][j][t] + inst.q_v[v]*x[v][i][j][t]);	
									travelAtCapacity[v][i][j][t].setName(ss.str().c_str());
								}else if (inst.typePort[i-1] == 1 && inst.typePort[j-1] == 0){
									ss1 << "travelEmpty_(" << i << "," << j << ")_" << t << "," << v;
									travelEmpty[v][i][j][t].setExpr(inst.q_v[v]*x[v][i][j][t] + fX[v][i][j][t]);
									travelEmpty[v][i][j][t].setName(ss1.str().c_str());
								}
							}else{ // when j is the sink node
								if(inst.typePort[i-1] == 0){
									ss << "travelAtCap_(" << i << ",snk)_" << t << "," << v;
									travelAtCapacity[v][i][j][t].setExpr(-fX[v][i][j][t] + inst.q_v[v]*x[v][i][j][t]);	
									travelAtCapacity[v][i][j][t].setName(ss.str().c_str());
								}else if (inst.typePort[i-1] == 1){ 
									ss1 << "travelEmpty_(" << i << ",snk)_" << t << "," << v;
									travelEmpty[v][i][j][t].setExpr(inst.q_v[v]*x[v][i][j][t] + fX[v][i][j][t]);
									travelEmpty[v][i][j][t].setName(ss1.str().c_str());
								}
							}
							model.add(travelAtCapacity[v][i][j][t]);
							model.add(travelEmpty[v][i][j][t]);
						}						
					}
				}
				
				flowCapacityOA[v][i] = IloRangeArray(env,T);
				flowCapacityOB[v][i] = IloRangeArray(env,T);
				flowCapacityW[v][i] = IloRangeArray(env,T);
				
				for(t=1;t<T;t++){
					stringstream ss, ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8, ss9;
					ss << "First_level_balance_" << v << "(" << i << "," << t << ")";
					ss1 << "Second_level_balance_" << v << "(" << i << "," << t << ")";
					ss2 << "link_balance_" << v << "(" << i << "," << t << ")";
					ss3 << "First_level_flow_" << v << "(" << i << "," << t << ")";
					ss4 << "Second_level_flow_" << v << "(" << i << "," << t << ")";
					
					if(hasArc[v][i][N-1][t]==1)
						expr_sinkFlow += x[v][i][N-1][t];
					
					expr_1stLevel.clear();
					expr_1stFlow.clear();
					expr_2ndLevel.clear();
					expr_2ndFlow.clear();
										
					//~ if(hasEnteringArc1st[v][i][t]==1){
						for(j=0;j<N;j++){
							if(j<N-1){ //No consider sink arc (first level balance)
								//If j is the source node and reach i in time t
								if(j==0 && inst.initialPort[v]+1 == i && t==inst.firstTimeAv[v]+1){
									expr_1stLevel += x[v][j][i][0];							
									expr_1stFlow += fX[v][j][i][0];
									fX[v][j][i][0].setBounds(inst.s_v0[v], inst.s_v0[v]); //Fixing initial inventory 
									//~ break; //If is the entering arc of vessel, there is no other arc entering port i for vessel v
								}
								else if (j>0){ //When j is a port
									if (t - inst.travelTime[v][j-1][i-1] >= 0){ //If is possible to exist an arc from j to i										
										if(hasArc[v][j][i][t-inst.travelTime[v][j-1][i-1]] == 1){ //If the arc exists
											expr_1stLevel += x[v][j][i][t-inst.travelTime[v][j-1][i-1]];									
											expr_1stFlow += fX[v][j][i][t-inst.travelTime[v][j-1][i-1]];
										}
									}
								}
							}
							if(j>0){ //No consider source arc (second level balance)
								//~ cout << v << " " << i << " " << j << " " << t << " " << hasArc[v][i][j][t] << endl;
								if (hasArc[v][i][j][t] == 1){
									//~ cout << v << " " << i << " " << j << " " << t << endl;
									expr_2ndLevel += x[v][i][j][t];
									expr_2ndFlow += fX[v][i][j][t];
								}
							}
						}
						if (t==1 || hasEnteringArc1st[v][i][t-1]==0){ //First time period or not entering arc in the previous time period
							firstLevelBalance[v][i][t].setExpr(expr_1stLevel - w[v][i][t] - oA[v][i][t]);
							secondLevelBalance[v][i][t].setExpr(oA[v][i][t]-oB[v][i][t] - expr_2ndLevel);
							linkBalance[v][i][t].setExpr(oA[v][i][t] - z[v][i][t]);
							firstLevelFlow[v][i][t].setExpr(expr_1stFlow - fW[v][i][t] - fOA[v][i][t]);
							secondLevelFlow[v][i][t].setExpr(fOA[v][i][t] + inst.delta[i-1]*f[v][i][t] - fOB[v][i][t] - expr_2ndFlow);
						}else if (t==T-1){ //Last time period
							firstLevelBalance[v][i][t].setExpr(expr_1stLevel + w[v][i][t-1] - oA[v][i][t]);
							secondLevelBalance[v][i][t].setExpr(oA[v][i][t] + oB[v][i][t-1] - expr_2ndLevel);
							linkBalance[v][i][t].setExpr(oA[v][i][t] + oB[v][i][t-1] - z[v][i][t]);
							firstLevelFlow[v][i][t].setExpr(expr_1stFlow + fW[v][i][t-1] - fOA[v][i][t]);
							secondLevelFlow[v][i][t].setExpr(fOA[v][i][t] + fOB[v][i][t-1] + inst.delta[i-1]*f[v][i][t] - expr_2ndFlow);
						}else{ //Other times
							firstLevelBalance[v][i][t].setExpr(expr_1stLevel + w[v][i][t-1] - w[v][i][t] - oA[v][i][t]);
							secondLevelBalance[v][i][t].setExpr(oA[v][i][t] + oB[v][i][t-1] - oB[v][i][t] - expr_2ndLevel);					
							linkBalance[v][i][t].setExpr(oA[v][i][t] + oB[v][i][t-1] - z[v][i][t]);
							firstLevelFlow[v][i][t].setExpr(expr_1stFlow + fW[v][i][t-1] - fW[v][i][t] - fOA[v][i][t]);
							secondLevelFlow[v][i][t].setExpr(fOA[v][i][t] + fOB[v][i][t-1] + inst.delta[i-1]*f[v][i][t] -fOB[v][i][t] - expr_2ndFlow);
						}
						firstLevelBalance[v][i][t].setName(ss.str().c_str());
						model.add(firstLevelBalance[v][i][t]);
						
						secondLevelBalance[v][i][t].setName(ss1.str().c_str());
						model.add(secondLevelBalance[v][i][t]);
						
						linkBalance[v][i][t].setName(ss2.str().c_str());
						model.add(linkBalance[v][i][t]);
						
						firstLevelFlow[v][i][t].setName(ss3.str().c_str());
						model.add(firstLevelFlow[v][i][t]);
						
						secondLevelFlow[v][i][t].setName(ss4.str().c_str());						
						model.add(secondLevelFlow[v][i][t]);
					
					
						ss5 << "fjmin_(" << i << "," << t << ")," << v;
						ss6 << "fjmax_(" << i << "," << t << ")," << v;
						IloNum minfMax = min(inst.f_max_jt[i-1][t-1], inst.q_v[v]);
						operationLowerLimit[v][i][t] = IloRange(env, 0, f[v][i][t] - inst.f_min_jt[i-1][t-1]*z[v][i][t] , IloInfinity, ss5.str().c_str());
						model.add(operationLowerLimit[v][i][t]);
						
						operationUpperLimit[v][i][t] = IloRange(env, -IloInfinity, f[v][i][t] - minfMax*z[v][i][t],	0, ss6.str().c_str());				
						model.add(operationUpperLimit[v][i][t]);
						
						ss7 << "flowLimitOA_"<<v<<","<<i<<","<<t;					
						flowCapacityOA[v][i][t] = IloRange(env, -IloInfinity, fOA[v][i][t]-inst.q_v[v]*oA[v][i][t], 0, ss7.str().c_str());					
						model.add(flowCapacityOA[v][i][t]);					
						
						if(t<T-1){ //Constraints with no last time index
							ss8 << "flowLimitOB_"<<v<<","<<i<<","<<t;
							ss9 << "flowLimitW_"<<v<<","<<i<<","<<t;
							flowCapacityOB[v][i][t] = IloRange(env, -IloInfinity, fOB[v][i][t]-inst.q_v[v]*oB[v][i][t], 0, ss8.str().c_str());
							flowCapacityW[v][i][t] = IloRange(env, -IloInfinity, fW[v][i][t]-inst.q_v[v]*w[v][i][t], 0, ss9.str().c_str());
							model.add(flowCapacityOB[v][i][t]);
							model.add(flowCapacityW[v][i][t]);
						}
					//~ }
				}				
			
			}			
			flowCapacityX[v][i] = IloArray<IloRangeArray>(env,N);
			for(j=0;j<N;j++){
				if(i != j){
					flowCapacityX[v][i][j] = IloRangeArray(env,T);
					for(t=0;t<T;t++){
						stringstream ss;
						ss << "flowLimitX_" << v << "," << i << "," << j << "," << t;
						flowCapacityX[v][i][j][t] = IloRange(env, -IloInfinity, fX[v][i][j][t] - inst.q_v[v]*x[v][i][j][t], 0, ss.str().c_str());
					}
					model.add(flowCapacityX[v][i][j]);
				}
			}
		
		}		
		stringstream ss;
		ss << "flowBalanceSink_" << v;
		sinkNodeBalance[v].setExpr(expr_sinkFlow);
		sinkNodeBalance[v].setName(ss.str().c_str());
		model.add(sinkNodeBalance);
	}
	expr_1stLevel.end();
	expr_2ndLevel.end();
	expr_2ndLevel.end();
	expr_2ndFlow.end();
	expr_sinkFlow.end();
	
	
	#ifndef NBranching
		
	#endif	
	
	#ifndef NValidInequalities
	
	#endif	
	
	#ifndef NOperateAndGo
	
	#endif	
	
	#ifndef N2PortNorevisit

	#endif
}

void Model::printSolution(IloEnv env, Instance& inst){
	int i,j,t,v,a;
	int J = inst.numTotalPorts;
	int T = inst.t + 1;			//Init time-periods in 1
	int V = inst.speed.getSize(); //# of vessels
	int N = J + 2; // Number of nodes including source and sink. (0,1...j,j+1)
	
	double costArc, revenue, costSpot, costAtempt;
	
	//Store binary variables
	xValue = IloArray<IloArray<IloArray<IloNumArray> > >(env,V) ;
	zValue = IloArray<IloArray<IloNumArray> >(env,V);
	wValue = IloArray<IloArray<IloNumArray> >(env,V);
	oAValue = IloArray<IloArray<IloNumArray> >(env,V);
	oBValue = IloArray<IloArray<IloNumArray> > (env,V);
	//Store fractional variables
	fValue = IloArray<IloArray<IloNumArray> >(env,V);
	fXValue = IloArray<IloArray<IloArray<IloNumArray> > >(env,V);
	foAValue = IloArray<IloArray<IloNumArray> > (env,V);
	foBValue = IloArray<IloArray<IloNumArray> >(env,V);
	fWValue = IloArray<IloArray<IloNumArray> >(env,V);
	sPValue = IloArray<IloNumArray>(env,N-1);
	alphaValue = IloArray<IloNumArray>(env,N-1);
	
	for (v=0;v<V;v++){
		xValue[v] = IloArray<IloArray<IloNumArray>>(env, N);
		zValue[v] = IloArray<IloNumArray>(env,N-1);
		wValue[v] = IloArray<IloNumArray>(env,N-1);
		oAValue[v] = IloArray<IloNumArray>(env,N-1);
		oBValue[v] = IloArray<IloNumArray>(env,N-1);
		
		fXValue[v] = IloArray<IloArray<IloNumArray>>(env, N);
		fValue[v] = IloArray<IloNumArray>(env,N-1);
		foAValue[v] = IloArray<IloNumArray>(env,N-1);
		foBValue[v] = IloArray<IloNumArray>(env,N-1);
		fWValue[v] = IloArray<IloNumArray>(env,N-1);		
		
		for (i=0;i<N;i++){
			if(i > 0 && i < N-1){ //Not consider sink node
				zValue[v][i] = IloNumArray(env,T);
				wValue[v][i] = IloNumArray(env,T);
				oAValue[v][i] = IloNumArray(env,T);
				oBValue[v][i] = IloNumArray(env,T);
				
				fValue[v][i] = IloNumArray(env,T);
				foAValue[v][i] = IloNumArray(env,T);
				foBValue[v][i] = IloNumArray(env,T);
				fWValue[v][i] = IloNumArray(env,T);
				
				for(t=1;t<T;t++){
					if(hasEnteringArc1st[v][i][t]){
						zValue[v][i][t] = cplex.getValue(z[v][i][t]);					
						cout << "Z value " << v << " " << i << " " << t << " " << zValue[v][i][t] << endl;
						if(zValue[v][i][t]>=0.1)
							costAtempt += (t-1)*inst.attemptCost;
						
						oAValue[v][i][t] = cplex.getValue(oA[v][i][t]);				
						
						fValue[v][i][t] = cplex.getValue(f[v][i][t]);
						if(zValue[v][i][t]>=0.1)
							revenue += fValue[v][i][t]*inst.r_jt[i-1][t-1];
						
						foAValue[v][i][t] = cplex.getValue(fOA[v][i][t]);						
						if(t<T-1){ //Variables not associated with last time index
							wValue[v][i][t] = cplex.getValue(w[v][i][t]);
							oBValue[v][i][t] = cplex.getValue(oB[v][i][t]);							
							foBValue[v][i][t] = cplex.getValue(fOB[v][i][t]);	
							fWValue[v][i][t] = cplex.getValue(fW[v][i][t]);						
						}
					}
				}
			}
			xValue[v][i] = IloArray<IloNumArray>(env,N);
			fXValue[v][i] = IloArray<IloNumArray>(env,N);
			for(j=0;j<N;j++){
				xValue[v][i][j] = IloNumArray(env,T);
				fXValue[v][i][j] = IloNumArray(env,T);
				for(t=0;t<T;t++){					
					if(hasArc[v][i][j][t] == 1){ // Only if the arc exists
						xValue[v][i][j][t] = cplex.getValue(x[v][i][j][t]);
						fXValue[v][i][j][t] = cplex.getValue(fX[v][i][j][t]);
						if(xValue[v][i][j][t] >= 0.1){
							costArc += arcCost[v][i][j][t];
						}
					}					
				}
			}
		}	
	}
	for(i=1;i<=J;i++){
		sPValue[i] = IloNumArray(env,T);
		alphaValue[i] = IloNumArray(env,T);
		for(t=0;t<T;t++){
			sPValue[i][t] = cplex.getValue(sP[i][t]);
			if(t>0){
				alphaValue[i][t] = cplex.getValue(alpha[i][t]);
				if (alphaValue[i][t] >= 0.1)
					costSpot += inst.p_jt[i-1][t-1];
				cout << "Port stock (" << i << "," << t << "): " << sPValue[i][t] << 
			" Alpha: " << alphaValue[i][t] << endl;
			}else
				cout << "Port stock (" << i << "," << t << "): " << sPValue[i][t] << endl;
		}
	}
	
	//Printing solution
	cout << "SOLUTION LOG: \n";
	for(v=0;v<V;v++){
		cout << "Vessel " << v << " route: " << endl;
		for(i=0;i<N-1;i++){			
			for(j=0;j<N;j++){	
				for(t=0;t<T;t++){
					if (i==0 && xValue[v][i][j][t]>=0.1) //Starting
						cout << "Depart at (" << i << "," << t << ") -> (" << j << "," << t+inst.firstTimeAv[v]+1 << ") with load: " << fXValue[v][i][j][t] << endl;
					else if(xValue[v][i][j][t] >= 0.1 && j < N-1) //Travelling
						cout << "Travel from (" << i << "," << t << ") -> (" << j << "," << t+inst.travelTime[v][i-1][j-1] << ") with load: " << fXValue[v][i][j][t] << " Cost: " << arcCost[v][i][j][t] << endl;
					else if (j == N-1 && xValue[v][i][j][t] >= 0.1) //Sink arcs
						cout << "Exit system at (" << i << "," << t << ") with load: " << fXValue[v][i][j][t] << " Cost: " << arcCost[v][i][j][t] << endl;
				}
			}
			if(i>0){
				for(t=1;t<T;t++){
					if(hasEnteringArc1st[v][i][t]){
						if(t<T-1 && wValue[v][i][t] >= 0.1)						
							cout << "Waiting in (" << i << "," << t << ") with load: " << fWValue[v][i][t] << endl;
						if(zValue[v][i][t] >= 0.1){
							if (t==1)
								cout << "Operate at (" << i << "," << t << ") - oA= " << oAValue[v][i][t] << "(" << foAValue[v][i][t] << ") - f= " << fValue[v][i][t] << endl;
							else{
								cout << "Operate at (" << i << "," << t << ") - oA= " << oAValue[v][i][t] << "(" << foAValue[v][i][t] << ") oB= " << oBValue[v][i][t-1] << "(" << foBValue[v][i][t-1] << 
								")- f= " << fValue[v][i][t] << endl;
							}
						}else{
							if (fValue[v][i][t] > 0.01)
								cout << "Operated but not in node (" << i << "," << t << ")_" << v << ": "  << fValue[v][i][t] << endl;
							assert(fValue[v][i][t] <= 0.01);
						}
					}
				}
			}
		}
	}
	cout << "Final value: " << costArc + costAtempt + costSpot - revenue << endl;
	cout << "Arcs: " << costArc << endl <<
			"Attempts: " << costAtempt << endl <<
			"Spots: " << costSpot << endl <<
			"Revenue: " << revenue << endl;

}
void Model::setParameters(IloEnv& env, const double& timeLimit, const double& gap=1e-04){
	cplex.exportModel("mip.lp");
	//~ env.setNormalizer(IloFalse);
	//Pressolve
	//~ cplex.setParam(IloCplex::PreInd,0);
	//~ cplex.setParam(IloCplex::RelaxPreInd,0);
	//~ cplex.setParam(IloCplex::PreslvNd,-1);
	
	//~ cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	//~ cplex.setParam(IloCplex::ConflictDisplay, 2); 
	//~ cplex.setParam(IloCplex::MIPDisplay, 1); 
	//~ cplex.setParam(IloCplex::MIPEmphasis, 4); // RINS
	//~ cplex.setParam(IloCplex::DiveType, 3); // Guied dive
	cplex.setParam(IloCplex::WorkMem, 8192);
	cplex.setParam(IloCplex::WorkDir, "workDir/");
	cplex.setParam(IloCplex::ClockType, 2);
	cplex.setParam(IloCplex::TiLim, timeLimit);
	//~ cplex.setParam(IloCplex::MIPEmphasis, 1); // 1 - Feasibility
	
	//~ cplex.setParam(IloCplex::NodeLim, 1);
	
	if (gap > 1e-04)
		cplex.setParam(IloCplex::EpGap, gap/100);
	//~ cplex.setParam(IloCplex::IntSolLim, 1); //Set number of solutions which must be found before stopping
	
	//Define high priority on the branching vairables
	#ifndef NBranching	
	for (int j=0;j<y.getSize();j++){
		cplex.setPriority(y[j],1);
	}	
	#endif
}
void mirp::milp(string file, const double& timeLimit, string optStr){
	///Time parameters
	Timer<chrono::milliseconds> timer_cplex;
	Timer<chrono::milliseconds> timer_global;
	timer_global.start();
	float global_time {0};
	float opt_time {0};
	
	///Read input files
	IloEnv env;
	try{
		Instance inst(env);
		inst.readInstance(env, file);
		
		int i,j,t,v,a;
		int J = inst.numTotalPorts;
		int T = inst.t;
		int V = inst.speed.getSize(); //# of vessels
		
		/// NEW MODEL
		Model model(env);
		model.buildModel(env,inst); 		
		model.setParameters(env, timeLimit);
		
		#ifndef NRelaxation
		//Remove the integrality constraints
		
		#endif

		//Start optiization	
		timer_cplex.start();
		if(model.cplex.solve()){
			global_time = timer_global.total();
			opt_time = timer_cplex.total();
			cout << model.cplex.getNbinVars() << " \t & " << model.cplex.getNintVars() << " \t & " << model.cplex.getNcols() << " \t & " << model.cplex.getNrows() << " \t & "
			<< opt_time/1000 << " \t & " << model.cplex.getObjValue() << " \t & " << model.cplex.getMIPRelativeGap()*100 
			<<"\% \t & " << model.cplex.getNiterations() << " \t & " << model.cplex.getNnodes() << " \t & " << endl;
			//~ 
			//~ model.cplex.setParam(IloCplex::IntSolLim, 2100000000); //Set number of solutions which must be found before stopping
			//~ model.cplex.setParam(IloCplex::TiLim, timeLimit-opt_time/1000);
			//~ 
			//~ timer_cplex.start();
			//~ timer_global.start();
			//~ model.cplex.solve();
			//~ global_time += timer_global.total();
			//~ opt_time += timer_cplex.total();
			//log << opt_time/1000 << " \t & " << model.cplex.getObjValue() << " \t & " << model.cplex.getMIPRelativeGap()*100 
			//<<"\% \t & " << model.cplex.getNiterations() << " \t & " << model.cplex.getNnodes() << " \t & \\\\";
					
			//Verify the solution
			
		}else{		
			cout << model.cplex.getStatus() << endl;		
			global_time = timer_global.total();
			opt_time = timer_cplex.total();
			cout << model.cplex.getNbinVars() << " \t & " << model.cplex.getNintVars() << " \t & " << model.cplex.getNcols() << " \t & " << model.cplex.getNrows() << " \t & "
			<< "\t & \t & \t & \t & \t & \t" <<  opt_time << " & \t & \t & " << model.cplex.getNiterations() << " \t & " << model.cplex.getNnodes() << " \\\\";		
		}
		model.printSolution(env, inst);
	}catch (IloException& e) {
	cerr << "Concert exception caught: " << e << endl;		
	e.end();
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}
	env.end();
	
	
	
}
