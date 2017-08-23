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
//~ #define NTravelAtCapacity
//~ #define NTravelEmpty
//~ #define NBerthLimit
#define NBranching
#define NBetas
#define WaitAfterOperate 				//If defined, allows a vessel to wait after operates at a port.
#define NKnapsackInequalities
//~ #define NSimplifyModel				//Remove arcs between port i and j for vessel v if min_f_i + min_f_j > Q_v

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
 * Just 1st interval have integer variables. 
 * Transition variables (x, w, wB, oB) are relaxed if the next time period belongs to a relaxed interval
 * Variables that are in the end block are not included in any constraint and objective function
 * Constraints of arcs/nodes in the and block are created, but are empty 
 */
void Model::buildFixAndRelaxModel(IloEnv& env, Instance inst, const double& nIntervals, const int& endBlock){
	int i, j,t,v,a;
	int timePerInterval = inst.t/nIntervals; //Number of time periods in each interval
	int J = inst.numTotalPorts;
	int T = inst.t+1;
	double intPart;
	int tOEB = inst.t - (timePerInterval*max(0.0,endBlock-modf(nIntervals, &intPart))) ; //Time periods out of End Block (index tOEB+1 is the first in the end block)
	int V = inst.speed.getSize(); //# of vessels
    int N = J + 2;
	cout << "Building model...\n Integer Block: [0," << timePerInterval << 
    "]\n Relaxed block: [" << timePerInterval+1 << "," << tOEB <<
    "]\n End block: [" << tOEB+1 << "," << T-1 << "]\n";
	
    ///Variables, converters and arrays to storage the values
    alpha = NumVarMatrix(env, N);    
	#ifndef NBetas
	beta = NumVarMatrix(env, N);
	#endif
    
    sP = NumVarMatrix(env,N);
	f = IloArray<NumVarMatrix> (env, V);
	fX = IloArray<IloArray<NumVarMatrix> >(env, V);	
	
	fOA = IloArray<NumVarMatrix> (env, V);
	#ifndef WaitAfterOperate
	fOB = IloArray<NumVarMatrix> (env, V);
	#endif
	fW = IloArray<NumVarMatrix> (env, V);

	hasArc = IloArray<IloArray<IloArray<IloIntArray> > >(env, V);
	arcCost = IloArray<IloArray<IloArray<IloNumArray> > >(env, V);
	hasEnteringArc1st = IloArray<IloArray<IloIntArray> >(env, V);
	
	x = IloArray<IloArray<IloArray<IloBoolVarArray> > >(env, V);    
	z = IloArray<IloArray<IloBoolVarArray> >(env,V);
    convertX = IloArray<IloArray<IloArray<IloArray<IloConversion> > > >(env,V);
    convertZ = IloArray<IloArray<IloArray<IloConversion> > >(env,V);
    xValue = IloArray<IloArray<IloArray<IloNumArray> > >(env,V);
    zValue = IloArray<IloArray<IloNumArray> >(env,V);
    
	#ifndef WaitAfterOperate
	oA = IloArray<IloArray<IloBoolVarArray> >(env,V);
	oB = IloArray<IloArray<IloBoolVarArray> >(env,V);
    convertOA = IloArray<IloArray<IloArray<IloConversion> > >(env,V);
    convertOB = IloArray<IloArray<IloArray<IloConversion> > >(env,V);
    oAValue = IloArray<IloArray<IloNumArray> >(env,V);
    oBValue = IloArray<IloArray<IloNumArray> >(env,V);
	#endif
	w = IloArray<IloArray<IloBoolVarArray> >(env,V);
    convertW = IloArray<IloArray<IloArray<IloConversion> > >(env,V);
    wValue = IloArray<IloArray<IloNumArray> >(env,V);
	
	#ifdef WaitAfterOperate
	fWB = IloArray<NumVarMatrix> (env, V);
	wB = IloArray<IloArray<IloBoolVarArray> >(env,V);
    convertWB = IloArray<IloArray<IloArray<IloConversion> > >(env,V);
    wBValue = IloArray<IloArray<IloNumArray> >(env,V);
	#endif
	
	//Additional variables for branching
	#ifndef NBranching
	#endif
	
    for(v=0;v<V;v++){
		f[v] = NumVarMatrix(env,N);
		fX[v] = IloArray<NumVarMatrix>(env,N);
		fOA[v] = NumVarMatrix(env,N);
		#ifndef WaitAfterOperate
		fOB[v] = NumVarMatrix(env,N);
		#endif
		fW[v] = NumVarMatrix(env,N);
		
		z[v] = IloArray<IloBoolVarArray>(env,N);
        convertZ[v] = IloArray<IloArray<IloConversion> >(env,N);
        zValue[v] = IloArray<IloNumArray>(env,N);
		#ifndef WaitAfterOperate
		oA[v] = IloArray<IloBoolVarArray>(env,N);
		oB[v] = IloArray<IloBoolVarArray>(env,N);
        convertOA[v] = IloArray<IloArray<IloConversion> >(env,N);
        convertOB[v] = IloArray<IloArray<IloConversion> >(env,N);
        oAValue[v] = IloArray<IloNumArray>(env,N);
        oBValue[v] = IloArray<IloNumArray>(env,N);
		#endif
		w[v] = IloArray<IloBoolVarArray>(env,N);
		x[v] = IloArray<IloArray<IloBoolVarArray> >(env, N);        
        convertW[v] = IloArray<IloArray<IloConversion> >(env,N);
        convertX[v] = IloArray<IloArray<IloArray<IloConversion> > >(env,N);
        wValue[v] = IloArray<IloNumArray>(env,N);
        xValue[v] = IloArray<IloArray<IloNumArray> >(env,N);
					
		hasArc[v] = IloArray<IloArray<IloIntArray> >(env, N);
		arcCost[v] = IloArray<IloArray<IloNumArray> >(env, N);
		hasEnteringArc1st[v] = IloArray<IloIntArray>(env, N);
		
		#ifdef WaitAfterOperate
		fWB[v] = NumVarMatrix(env,N);
		wB[v] = IloArray<IloBoolVarArray>(env,N);
        convertWB[v] = IloArray<IloArray<IloConversion> >(env,N);
        wBValue[v] = IloArray<IloNumArray>(env,N);
		#endif

		for(j=0;j<N;j++){
			x[v][j] = IloArray<IloBoolVarArray>(env,N);
            convertX[v][j] = IloArray<IloArray<IloConversion> >(env,N);
            xValue[v][j] = IloArray<IloNumArray>(env,N);
			hasArc[v][j] = IloArray<IloIntArray>(env,N);
			arcCost[v][j] = IloArray<IloNumArray>(env,N);
			fX[v][j] = IloArray<IloNumVarArray>(env, N);			
            for(i=0;i<N;i++){ //TODO Optimize: use i=j+1 (there are consequences on build arcs)
				x[v][j][i] = IloBoolVarArray(env,T);
                convertX[v][j][i] = IloArray<IloConversion>(env,T);
                xValue[v][j][i] = IloNumArray(env, T);
				hasArc[v][j][i] = IloIntArray(env,T);
				arcCost[v][j][i] = IloNumArray(env,T);
				fX[v][j][i] = IloNumVarArray(env, T);                
				for(t=0;t<T;t++){
					stringstream ss;
					ss << "x_"<< v <<","<< j <<","<< i <<","<<t;
					x[v][j][i][t].setName(ss.str().c_str());
					ss.str(string());
                    if(j > 0 && j <= J){ //If j is a port
                        if (i > 0 && i <= J){ //If i is a port
                            if(t + inst.travelTime[v][j-1][i-1] > timePerInterval){ 	//If the arrive time is out of first interval, relax it                                
                                convertX[v][j][i][t] = IloConversion(env, x[v][j][i][t], ILOFLOAT);
                                model.add(convertX[v][j][i][t]);
                            }
                        }else if (i == N-1){ // If i is the sink node
                            if(t > timePerInterval){ 
                                convertX[v][j][i][t] = IloConversion(env, x[v][j][i][t], ILOFLOAT);
                                model.add(convertX[v][j][i][t]);
                            }
                        }
                    }
					ss << "fX_"<<v<<","<<j<<","<<i<<","<<t;
					fX[v][j][i][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
				}
			}

			if(j>0 && j<=J){ //Only considering Ports
				f[v][j] = IloNumVarArray(env, T);
				fOA[v][j] = IloNumVarArray(env, T);
				#ifndef WaitAfterOperate
				fOB[v][j] = IloNumVarArray(env, T);
				#endif
				fW[v][j] = IloNumVarArray(env, T);
				
				z[v][j] = IloBoolVarArray(env, T);
                convertZ[v][j] = IloArray<IloConversion>(env,T);
                zValue[v][j] = IloNumArray(env,T);
				#ifndef WaitAfterOperate
				oA[v][j] = IloBoolVarArray(env, T);
				oB[v][j] = IloBoolVarArray(env, T);
                convertOA[v][j] = IloArray<IloConversion>(env,T);
                convertOB[v][j] = IloArray<IloConversion>(env,T);
                oAValue[v][j] = IloNumArray(env,T);
                oBValue[v][j] = IloNumArray(env,T);
				#endif
				w[v][j] = IloBoolVarArray(env, T);
				convertW[v][j] = IloArray<IloConversion>(env,T);
				wValue[v][j] = IloNumArray(env,T);
				hasEnteringArc1st[v][j] = IloIntArray(env,T);
				#ifdef WaitAfterOperate
				wB[v][j] = IloBoolVarArray(env, T);
                convertWB[v][j] = IloArray<IloConversion>(env,T);
                wBValue[v][j] = IloNumArray(env,T);
				fWB[v][j] = IloNumVarArray(env, T);
				#endif
				
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
                    if(t > timePerInterval){
                        convertZ[v][j][t] = IloConversion(env, z[v][j][t], ILOFLOAT);
                        model.add(convertZ[v][j][t]);
                    }
					#ifndef WaitAfterOperate
					ss.str(string());
					ss << "oA_(" << j << "," << t << ")," << v;
					oA[v][j][t].setName(ss.str().c_str());
                    if(t > timePerInterval){
                        convertOA[v][j][t] = IloConversion(env, oA[v][j][t], ILOFLOAT);
                        model.add(convertOA[v][j][t]);
                    }
					#endif
					if(t<T-1){ //Variables that haven't last time index
						ss.str(string());
						ss << "fW_(" << j << "," << t << ")," << v;
						fW[v][j][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
                        
                        ss.str(string());
						ss << "w_(" << j << "," << t << ")," << v;
						w[v][j][t].setName(ss.str().c_str());
                        if(t >= timePerInterval){       //Note the use of >= when considering 'horizontal' transition arcs
                            convertW[v][j][t] = IloConversion(env, w[v][j][t], ILOFLOAT);
                            model.add(convertW[v][j][t]);
                        }

						#ifndef WaitAfterOperate
                        ss.str(string());
						ss << "fOB_(" << j << "," << t << ")," << v;
						fOB[v][j][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
                        
						ss.str(string());
						ss << "oB_(" << j << "," << t << ")," << v;
						oB[v][j][t].setName(ss.str().c_str());
                        if(t >= timePerInterval){
                            convertOB[v][j][t] = IloConversion(env, oB[v][j][t], ILOFLOAT);
                            model.add(convertOB[v][j][t]);
                        }
						#endif
						
						#ifdef WaitAfterOperate
                        ss.str(string());
                        ss << "wB_(" << j << "," << t << ")," << v;
                        wB[v][j][t].setName(ss.str().c_str());
                        if(t >= timePerInterval){
                            convertWB[v][j][t] = IloConversion(env, wB[v][j][t], ILOFLOAT);
                            model.add(convertWB[v][j][t]);
                        }
                        ss.str(string());
                        ss << "fWB_(" << j << "," << t << ")," << v;
                        fWB[v][j][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
                        #endif
					}
				}
			}
		}
	}	
	for(i=1;i<N-1;i++){
		alpha[i] = IloNumVarArray(env,T);
		sP[i] = IloNumVarArray(env, T);	
        #ifndef NBetas
        beta[i] = IloNumVarArray(env,T);
        #endif
		for(t=0;t<T;t++){
			stringstream ss;
			if (t == 0){
				ss << "sP_(" << i << ",0)";
				sP[i][t] = IloNumVar(env,inst.s_j0[i-1], inst.s_j0[i-1], ss.str().c_str()); //Initial inventory
			}else{
				ss.str(string());
				ss << "sP_(" << i << "," << t << ")";
				sP[i][t] = IloNumVar(env, inst.sMin_jt[i-1][0], inst.sMax_jt[i-1][0], ss.str().c_str()); //As the port capacity is fixed, always used the data from index 0
				ss.str(string());
				ss << "alpha_(" << i << "," << t << ")";
				alpha[i][t] = IloNumVar(env, 0, inst.alp_max_jt[i-1][t-1], ss.str().c_str());
                #ifndef NBetas
                ss.str(string());
				ss << "beta_(" << i << "," << t << ")";
				beta[i][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
                #endif
			}
		}
	}
	
	///Objective function
	IloExpr expr(env);
	IloExpr expr1(env);	
	for(v=0;v<V;v++){
		///Builds the arcs between nodes. Only arcs in the integer and relaxed intervals are added							//2nd term
		//Source arcs -- all are included, assuming that the first time available is in the integer or relaxed block
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
				#ifdef NSimplifyModel
				if (i != j){
				#endif
				#ifndef NSimplifyModel
				if (i != j && (inst.typePort[i-1] != inst.typePort[j-1] || inst.f_min_jt[i-1][t-1] + inst.f_min_jt[j-1][t-1] <= inst.q_v[v]) ){
				#endif
					int t2 = t + inst.travelTime[v][i-1][j-1]; 
					if (t2<T){ 	//If exists time to reach port j 
						double arc_cost;
						if (inst.typePort[i-1]==1 && inst.typePort[j-1]==0){ //If port i is consuming and j is a producer port
							arc_cost = inst.costKm[v]*inst.distanceMatrix[i-1][j-1]*(1-inst.trav_empt[v]) + inst.portFee[j-1];
						}else{
							arc_cost = inst.costKm[v]*inst.distanceMatrix[i-1][j-1] + inst.portFee[j-1];
						}
						hasArc[v][i][j][t] = 1;	
						arcCost[v][i][j][t] = arc_cost;
						if(t2 <= tOEB)      //If arriving node is in the model
                            expr1 += arc_cost*x[v][i][j][t];
                        hasEnteringArc1st[v][j][t2] = 1;
						if (t2+1<T-1) 			//Waiting arc 
							hasEnteringArc1st[v][j][t2+1] = 1;

						//Sink arc from port j 
						hasArc[v][j][N-1][t2] = 1;
						arcCost[v][j][N-1][t2] = -(T-t2-1)*inst.perPeriodRewardForFinishingEarly;
						if(t2 <= tOEB)      //If arriving node is in the model
                            expr1 += arcCost[v][j][N-1][t2]*x[v][j][N-1][t2];
						//when time t reach a time that can be built a mirror between j and i
						if (t >= inst.firstTimeAv[v]+1 + inst.travelTime[v][i-1][j-1]){ 
							if (inst.typePort[j-1]==1 && inst.typePort[i-1]==0){ 		  //If port j is consuming and i is a producer port
								arc_cost = inst.costKm[v]*inst.distanceMatrix[j-1][i-1]*(1-inst.trav_empt[v]) + inst.portFee[i-1];		
							}else{
								arc_cost = inst.costKm[v]*inst.distanceMatrix[j-1][i-1] + inst.portFee[i-1];
							}
							hasArc[v][j][i][t] = 1;
							arcCost[v][j][i][t] = arc_cost;
							hasEnteringArc1st[v][i][t2] = 1;
							if(t2 <= tOEB)      //If arriving node is in the model
                                expr1 += arc_cost*x[v][j][i][t];
						}
						//Create arc from j,t2 to others ports (j2) in time t3
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
									if(t3 <= tOEB)      //Ir arriving node is in the model
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
			if(t <= tOEB)      //If arriving node is in the model
                expr1 += arcCost[v][i][N-1][t]*x[v][i][N-1][t];
		}

		for(i=1;i<N-1;i++){		//Only considering ports
			for(t=1;t<=tOEB;t++){
				if(hasEnteringArc1st[v][i][t]){
					expr += inst.r_jt[i-1][t-1]*f[v][i][t];									//1st term
					expr1 += (t-1)*inst.attemptCost*z[v][i][t];								//3rd term
				}
			}
		}
	}
	
	for(j=1;j<N-1;j++){
		for(t=1;t<=tOEB;t++){			
			expr1 += inst.p_jt[j-1][t-1]*alpha[j][t];									//4rd term
            #ifndef NBetas
            expr1 += 1000*beta[j][t];
            #endif
		}
	}
	obj.setExpr(expr1-expr);
	model.add(obj);
	
	///Constraints
	//Flow balance source and sink nodes
	IloExpr expr_sinkFlow(env);
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
	#ifndef WaitAfterOperate
	linkBalance = IloArray<IloArray<IloRangeArray> > (env, V);
	#endif
	
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
	flowCapacityW = IloArray<IloArray<IloRangeArray> >(env,V);
	#ifndef WaitAfterOperate
	flowCapacityOB = IloArray<IloArray<IloRangeArray> >(env,V);
	#endif
	#ifdef WaitAfterOperate
	flowCapacityWB = IloArray<IloArray<IloRangeArray> >(env,V);
	#endif
	
	#ifndef NKnapsackInequalities
	knapsack_P_1 = IloArray<IloRangeArray> (env, J+1);			//Altough created for all J ports, it is only considered for loading(production) or consumption(discharging) ports
	knapsack_P_2 = IloArray<IloRangeArray> (env, J+1);	
	knapsack_D_1 = IloArray<IloRangeArray> (env, J+1);	
	knapsack_D_2 = IloArray<IloRangeArray> (env, J+1);	
	#ifndef WaitAfterOperate
	knapsack_D_3 = IloArray<IloRangeArray> (env, J+1);	
	#endif
	#endif
	
	for(i=1;i<N-1;i++){
		stringstream ss1;
		berthLimit[i] = IloRangeArray(env,T);
		expr_cumSlack.clear();
		portInventory[i] = IloRangeArray(env,T);				
		#ifndef NKnapsackInequalities
		if (inst.typePort[i-1] == 0){ 	//Loading port
			knapsack_P_1[i] = IloRangeArray(env, (T-3)*2);		//(T-2)*2 = Number of combinations 1,...t + t,...,T for all t \in T. Using T-3 because de increment
			knapsack_P_2[i] = IloRangeArray(env, (T-3)*2);
		}else{							//Discharging port
			knapsack_D_1[i] = IloRangeArray(env, (T-3)*2);
			knapsack_D_2[i] = IloRangeArray(env, (T-3)*2);
			#ifndef WaitAfterOperate
			knapsack_D_3[i] = IloRangeArray(env, (T-3)*2);
			#endif
		}        
		int it,it2=0,k,l;
		IloExpr expr_kP1_LHS(env), expr_kP2_LHS(env), expr_kD1_LHS(env), expr_kD2_LHS(env);		
		for(it=0;it<((T-3)*2)-(T-1-tOEB);it++){	//For each valid inequality - limited to the constraint that uses the interval [tOEB-1, tOEB]
			double kP1_RHS=0, kP2_RHS=0, kD1_RHS=0, kD2_RHS=0, sum_alphaMax=0, alphaUB=0;
			expr_kP1_LHS.clear();
			expr_kP2_LHS.clear();
			expr_kD1_LHS.clear();
			expr_kD2_LHS.clear();
			//Definining the size of set T =[k,l] - CAUTION: It is inverse of the paper of Agra et al (2013)
			k=1;
			l=tOEB;
			if(it<tOEB-2){
				l = it+2;
				it2++;		//It2 gets the same value of it until it reach the value tOEB-2
			}else if(it == tOEB-2){
				it = T-3;	//it jumps to half of array of the IloRangeArrey (starting the constraints of range [t,...|T|], t \in T
				k = it2+4-l;
				it2++;
			}
			else{				
				k = it2+4-l;
				it2++;
			}			
            for(v=0;v<V;v++){
				#ifdef WaitAfterOperate
				expr_kP1_LHS += wB[v][i][k];
				if(k>1 || hasEnteringArc1st[v][i][k-1]==1)
					expr_kD1_LHS += w[v][i][k-1] + wB[v][i][k-1];
				#endif
				#ifndef WaitAfterOperate
				expr_kP1_LHS += oB[v][i][k];
				if(k>1 || hasEnteringArc1st[v][i][k-1]==1)
					expr_kD1_LHS += w[v][i][k-1] + oB[v][i][k-1];
				#endif
				#ifndef WaitAfterOperate
				if(k>1 || hasEnteringArc1st[v][i][k-1]==1)
					expr_kD2_LHS += oB[v][i][k-1];
				#endif
				for(t=k;t<=l;t++){
					expr_kP2_LHS += z[v][i][t];	//It is used for both loading and discharging ports
					#ifndef WaitAfterOperate
					expr_kD2_LHS += oA[v][i][t];
					#endif
					#ifdef WaitAfterOperate
					expr_kD2_LHS += z[v][i][t];
					#endif
					for(j=0;j<N;j++){						
						if(j==0 && inst.initialPort[v]+1 == i && t==inst.firstTimeAv[v]+1) 	//Source arc
							expr_kD1_LHS += x[v][j][i][0];								
						else if (j == N-1 && hasArc[v][i][j][t] == 1)						//Sink arc
							expr_kP1_LHS += x[v][i][j][t];
						else if (j > 0 && j < N-1){											//"Normal" arcs	
							if(i != j){
								if(hasArc[v][i][j][t] == 1 && t+inst.travelTime[v][i-1][j-1] <= tOEB)  //If arc exists and arrives at port j in the integer or relaxed block
									expr_kP1_LHS += x[v][i][j][t];
								if(t - inst.travelTime[v][j-1][i-1] >= 0){					//Only if it is possible an entering arc due the travel time
									if(hasArc[v][j][i][t-inst.travelTime[v][j-1][i-1]] == 1)
										expr_kD1_LHS += x[v][j][i][t-inst.travelTime[v][j-1][i-1]];
								}
							}
						}
					}
				}
			}
			for(t=k;t<=l;t++){
				//~ sum_alphaMax += inst.alp_max_jt[i-1][t-1];
				kP1_RHS += inst.d_jt[i-1][t-1];
				kD1_RHS += inst.d_jt[i-1][t-1];
			}
			//~ alphaUB = min(sum_alphaMax, inst.alp_max_j[i-1]);
			kP1_RHS += -inst.sMax_jt[i-1][0] - alphaUB;
			kP2_RHS = ceil(kP1_RHS/inst.f_max_jt[i-1][0]);
			kP1_RHS = ceil(kP1_RHS/inst.maxCapacity);
			
			if(k==1)
				kD1_RHS += -inst.s_j0[i-1] + inst.sMin_jt[i-1][0] - alphaUB;
			else
				kD1_RHS += -inst.sMax_jt[i-1][k-1] + inst.sMin_jt[i-1][0] - alphaUB;
			kD2_RHS = ceil(kD1_RHS/inst.f_max_jt[i-1][0]);
			kD1_RHS = ceil(kD1_RHS/inst.maxCapacity);
			
			stringstream ss, ss1, ss2;
			if(inst.typePort[i-1] == 0){
				ss << "knpasackP1_" << i << "," << it;
				knapsack_P_1[i][it] = IloRange(env, kP1_RHS, expr_kP1_LHS, IloInfinity, ss.str().c_str());
				ss1 << "knapsackP2_" << i << "," << it;
				knapsack_P_2[i][it] = IloRange(env, kP2_RHS, expr_kP2_LHS, IloInfinity, ss1.str().c_str());
                model.add(knapsack_P_1[i][it]);
                model.add(knapsack_P_2[i][it]);
			}else{
				ss << "knpasackD1_" << i << "," << it;
				knapsack_D_1[i][it] = IloRange(env, kD1_RHS, expr_kD1_LHS, IloInfinity, ss.str().c_str());
				ss1 << "knpasackD2_" << i << "," << it;
				knapsack_D_2[i][it] = IloRange(env, kD1_RHS, expr_kD2_LHS, IloInfinity, ss1.str().c_str());
                model.add(knapsack_D_1[i][it]);
                model.add(knapsack_D_2[i][it]);
				#ifndef WaitAfterOperate
				ss2 << "knpasackD3_" << i << "," << it;
				knapsack_D_3[i][it] = IloRange(env, kD2_RHS, expr_kP2_LHS, IloInfinity, ss2.str().c_str());
                model.add(knapsack_D_3[i][it]);
				#endif
			}
		}
		#endif
		
		for(t=1;t<=tOEB;t++){
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
			#ifndef NBetas
            expr_invBalancePort += -beta[i][t];
            #endif
            
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
		#ifndef WaitAfterOperate
		linkBalance[v] = IloArray<IloRangeArray>(env,J+1); 
		#endif
		
		travelAtCapacity[v] = IloArray<IloArray<IloRangeArray> > (env, N-1);
		travelEmpty[v] = IloArray<IloArray<IloRangeArray> > (env, N-1);
		
		operationLowerLimit[v] = IloArray<IloRangeArray> (env,N-1);
		operationUpperLimit[v] = IloArray<IloRangeArray> (env,N-1);
		
		flowCapacityX[v] = IloArray<IloArray<IloRangeArray> >(env,N);
		flowCapacityOA[v] = IloArray<IloRangeArray> (env,N-1);
		#ifndef WaitAfterOperate
		flowCapacityOB[v] = IloArray<IloRangeArray> (env,N-1);
		#endif
		flowCapacityW[v] = IloArray<IloRangeArray> (env,N-1);
		#ifdef WaitAfterOperate
		flowCapacityWB[v] = IloArray<IloRangeArray> (env,N-1);
		#endif
		
		for(i=0;i<N;i++){
			if(i>0 && i <= J){ //Only considering ports				
				firstLevelBalance[v][i] = IloRangeArray(env,T,0,0); //Id 0 is not used
				secondLevelBalance[v][i] = IloRangeArray(env,T,0,0); 
				firstLevelFlow[v][i] = IloRangeArray(env,T,0,0); 
				secondLevelFlow[v][i] = IloRangeArray(env,T,0,0); 
				#ifndef WaitAfterOperate
				linkBalance[v][i] = IloRangeArray(env,T,0,0); 
				#endif
				travelAtCapacity[v][i] = IloArray<IloRangeArray> (env, N);
				travelEmpty[v][i] = IloArray<IloRangeArray> (env, N);
				
				operationLowerLimit[v][i] = IloRangeArray (env,T);
				operationUpperLimit[v][i] = IloRangeArray (env,T);
				
				for(j=0;j<N;j++){
					if(j>0){ //Considering ports and sink node
						travelAtCapacity[v][i][j] = IloRangeArray(env, T, -IloInfinity, 0);
						travelEmpty[v][i][j] = IloRangeArray (env, T, -IloInfinity, inst.q_v[v]);
						for(t=1;t<=tOEB;t++){
							stringstream ss, ss1;
                            if (j<=J){				//When j is a port
                                if(t + inst.travelTime[v][i-1][j-1] <= tOEB){ //Only if arc departs and arrives in the integer or relaxed block
                                    if(inst.typePort[i-1] == 0 && inst.typePort[j-1] == 1){
                                        ss << "travelAtCap_(" << i << "," << j << ")_" << t << "," << v;
                                        travelAtCapacity[v][i][j][t].setExpr(-fX[v][i][j][t] + inst.q_v[v]*x[v][i][j][t]);	
                                        travelAtCapacity[v][i][j][t].setName(ss.str().c_str());
                                        model.add(travelAtCapacity[v][i][j][t]);
                                    }else if (inst.typePort[i-1] == 1 && inst.typePort[j-1] == 0){
                                        ss1 << "travelEmpty_(" << i << "," << j << ")_" << t << "," << v;
                                        travelEmpty[v][i][j][t].setExpr(inst.q_v[v]*x[v][i][j][t] + fX[v][i][j][t]);
                                        travelEmpty[v][i][j][t].setName(ss1.str().c_str());
                                        model.add(travelEmpty[v][i][j][t]);
                                    }
                                }
							}else{ // when j is the sink node
								if(inst.typePort[i-1] == 0){
									ss << "travelAtCap_(" << i << ",snk)_" << t << "," << v;
									travelAtCapacity[v][i][j][t].setExpr(-fX[v][i][j][t] + inst.q_v[v]*x[v][i][j][t]);	
									travelAtCapacity[v][i][j][t].setName(ss.str().c_str());
                                    model.add(travelAtCapacity[v][i][j][t]);
								}else if (inst.typePort[i-1] == 1){ 
									ss1 << "travelEmpty_(" << i << ",snk)_" << t << "," << v;
									travelEmpty[v][i][j][t].setExpr(inst.q_v[v]*x[v][i][j][t] + fX[v][i][j][t]);
									travelEmpty[v][i][j][t].setName(ss1.str().c_str());
                                    model.add(travelEmpty[v][i][j][t]);
								}
							}
						}
					}
				}
				flowCapacityOA[v][i] = IloRangeArray(env,T);
				#ifndef WaitAfterOperate
				flowCapacityOB[v][i] = IloRangeArray(env,T);
				#endif
				flowCapacityW[v][i] = IloRangeArray(env,T);
				#ifdef WaitAfterOperate
				flowCapacityWB[v][i] = IloRangeArray(env,T);
				#endif
				
				for(t=1;t<=tOEB;t++){
					stringstream ss, ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8, ss9, ss10;
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

					for(j=0;j<N;j++){
						if(j<N-1){ //No consider sink arc (first level balance)
							//If j is the source node and reach i in time t
							if(j==0 && inst.initialPort[v]+1 == i && t==inst.firstTimeAv[v]+1){
								expr_1stLevel += x[v][j][i][0];
								expr_1stFlow += fX[v][j][i][0];
								fX[v][j][i][0].setBounds(inst.s_v0[v], inst.s_v0[v]); //Fixing initial inventory 
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
							if (hasArc[v][i][j][t] == 1){
                                if(j == N-1){ //If j is the sink node, add to expr
                                    expr_2ndLevel += x[v][i][j][t];
                                    expr_2ndFlow += fX[v][i][j][t];
                                }else if (t + inst.travelTime[v][i-1][j-1] <= tOEB){ //If j is a port, it is then necessary that the arrival is in the model
                                    expr_2ndLevel += x[v][i][j][t];
                                    expr_2ndFlow += fX[v][i][j][t];
                                }
							}
						}
					}
					IloExpr expr_link;
					#ifndef WaitAfterOperate
					if (t==1 || (t < tOEB && hasEnteringArc1st[v][i][t-1]==0)){ //First time period or not last time period nor entering arc in the previous time period
						expr_1stLevel += - w[v][i][t] - oA[v][i][t];
                        expr_2ndLevel += -oA[v][i][t] + oB[v][i][t];
                        expr_link = oA[v][i][t] - z[v][i][t];						
                        expr_1stFlow+= -fW[v][i][t] - fOA[v][i][t];						
						expr_2ndFlow += -fOA[v][i][t] - inst.delta[i-1]*f[v][i][t] + fOB[v][i][t];
					}else if (t==tOEB){ //Last time period
                        if (hasEnteringArc1st[v][i][t-1]==1){
                            expr_1stLevel += w[v][i][t-1] - oA[v][i][t];
                            expr_2ndLevel += -oA[v][i][t] -oB[v][i][t-1];
                            expr_link = oA[v][i][t] + oB[v][i][t-1] - z[v][i][t];   
                            expr_1stFlow += fW[v][i][t-1] - fOA[v][i][t];
                            expr_2ndFlow += -fOA[v][i][t] - fOB[v][i][t-1] - inst.delta[i-1]*f[v][i][t];
                        }else{
                            expr_1stLevel += - oA[v][i][t];
                            expr_2ndLevel += -oA[v][i][t];
                            expr_link = oA[v][i][t] - z[v][i][t];   
                            expr_1stFlow += - fOA[v][i][t];
                            expr_2ndFlow += -fOA[v][i][t] - inst.delta[i-1]*f[v][i][t];
                        }
					}else{ //Other times
						expr_1stLevel += w[v][i][t-1] - w[v][i][t] - oA[v][i][t];
						expr_2ndLevel += - oA[v][i][t] - oB[v][i][t-1] + oB[v][i][t];
						expr_link = oA[v][i][t] + oB[v][i][t-1] - z[v][i][t];
						expr_1stFlow += fW[v][i][t-1] - fW[v][i][t] - fOA[v][i][t];
						expr_2ndFlow += -fOA[v][i][t] - fOB[v][i][t-1] - inst.delta[i-1]*f[v][i][t] + fOB[v][i][t];	
					}
					#endif
					#ifdef WaitAfterOperate
					if (t==1 || (t < tOEB && hasEnteringArc1st[v][i][t-1]==0)){ //First time period or not last time period nor entering arc in the previous time period
						expr_1stLevel += - w[v][i][t] - z[v][i][t];
						expr_2ndLevel += -z[v][i][t] + wB[v][i][t];	
						expr_1stFlow+= -fW[v][i][t] - fOA[v][i][t];	
						expr_2ndFlow += -fOA[v][i][t] - inst.delta[i-1]*f[v][i][t] + fWB[v][i][t];
					}else if (t==tOEB){ //Last time period
						if (hasEnteringArc1st[v][i][t-1]==1){
                            expr_1stLevel += w[v][i][t-1] - z[v][i][t] + wB[v][i][t-1];                                
                            expr_1stFlow += fW[v][i][t-1] - fOA[v][i][t] + fWB[v][i][t-1];                                
                        }else{
                            expr_1stLevel += - z[v][i][t];                                
                            expr_1stFlow += - fOA[v][i][t];                                
                        }
                        expr_2ndLevel += -z[v][i][t];
                        expr_2ndFlow += -fOA[v][i][t] - inst.delta[i-1]*f[v][i][t];
					}else{ //Other times
						expr_1stLevel += w[v][i][t-1] - w[v][i][t] - z[v][i][t] + wB[v][i][t-1];
						expr_2ndLevel += - z[v][i][t] + wB[v][i][t];	
						expr_1stFlow += fW[v][i][t-1] - fW[v][i][t] - fOA[v][i][t] + fWB[v][i][t-1];
						expr_2ndFlow += -fOA[v][i][t] - inst.delta[i-1]*f[v][i][t] + fWB[v][i][t];
					}
					#endif
					
					firstLevelBalance[v][i][t].setExpr(expr_1stLevel);
					secondLevelBalance[v][i][t].setExpr(expr_2ndLevel);
					#ifndef WaitAfterOperate
					linkBalance[v][i][t].setExpr(expr_link);
					#endif
					firstLevelFlow[v][i][t].setExpr(expr_1stFlow);
					secondLevelFlow[v][i][t].setExpr(expr_2ndFlow);
					
					firstLevelBalance[v][i][t].setName(ss.str().c_str());
					model.add(firstLevelBalance[v][i][t]);
					
					secondLevelBalance[v][i][t].setName(ss1.str().c_str());
					model.add(secondLevelBalance[v][i][t]);
					
					#ifndef WaitAfterOperate
					linkBalance[v][i][t].setName(ss2.str().c_str());
					model.add(linkBalance[v][i][t]);
					#endif
					
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
					#ifndef WaitAfterOperate
					flowCapacityOA[v][i][t] = IloRange(env, -IloInfinity, fOA[v][i][t]-inst.q_v[v]*oA[v][i][t], 0, ss7.str().c_str());					
					#endif
					#ifdef WaitAfterOperate
					flowCapacityOA[v][i][t] = IloRange(env, -IloInfinity, fOA[v][i][t]-inst.q_v[v]*z[v][i][t], 0, ss7.str().c_str());					
					#endif
					model.add(flowCapacityOA[v][i][t]);					
					
					if(t<tOEB){ //Constraints with no last time index
						#ifndef WaitAfterOperate
						ss8 << "flowLimitOB_"<<v<<","<<i<<","<<t;
						flowCapacityOB[v][i][t] = IloRange(env, -IloInfinity, fOB[v][i][t]-inst.q_v[v]*oB[v][i][t], 0, ss8.str().c_str());
						model.add(flowCapacityOB[v][i][t]);
						#endif
						ss9 << "flowLimitW_"<<v<<","<<i<<","<<t;
						flowCapacityW[v][i][t] = IloRange(env, -IloInfinity, fW[v][i][t]-inst.q_v[v]*w[v][i][t], 0, ss9.str().c_str());							
						model.add(flowCapacityW[v][i][t]);
						#ifdef WaitAfterOperate
						ss10 << "flowLimitWB_"<<v<<","<<i<<","<<t;
						flowCapacityWB[v][i][t] = IloRange(env, -IloInfinity, fWB[v][i][t]-inst.q_v[v]*wB[v][i][t], 0, ss10.str().c_str());
						model.add(flowCapacityWB[v][i][t]);
						#endif							
					}
				}													
				flowCapacityX[v][i] = IloArray<IloRangeArray>(env,N);
				for(j=1;j<N;j++){ //No arriving arc in source arc j=0
					if(i != j && i < N-1){ //There is no departing arc from sink node
						flowCapacityX[v][i][j] = IloRangeArray(env,T);
						for(t=0;t<=tOEB;t++){
							stringstream ss;
							ss << "flowLimitX_" << v << "," << i << "," << j << "," << t;
							if(j == N-1){ //If j is sink node                            
								flowCapacityX[v][i][j][t] = IloRange(env, -IloInfinity, fX[v][i][j][t] - inst.q_v[v]*x[v][i][j][t], 0, ss.str().c_str());
								model.add(flowCapacityX[v][i][j][t]);
							}else if (i>0 && t + inst.travelTime[v][i-1][j-1] <= tOEB){ //If x variable reach j in the time period
								flowCapacityX[v][i][j][t] = IloRange(env, -IloInfinity, fX[v][i][j][t] - inst.q_v[v]*x[v][i][j][t], 0, ss.str().c_str());
								model.add(flowCapacityX[v][i][j][t]);
							}
						}
					}
				}		
			}
		}		
		stringstream ss;
		ss << "flowBalanceSink_" << v;
		sinkNodeBalance[v].setExpr(expr_sinkFlow);
		sinkNodeBalance[v].setName(ss.str().c_str());
		model.add(sinkNodeBalance[v]);
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
/* Param p = 0 If algorithm can extract var solution values; 1 otherwise*/
void Model::fixSolution(IloEnv& env, Instance inst, const int& tS, const int& tF,const int& p){
	int J = inst.numTotalPorts;
	int N = J+1;
	
	if(p == 0){
		getSolValsW(env, inst, tS, tF);	//Get only values of the interval
	}
	
	//Fix solution - If p == 1, consider the values from the previous p==0
	for(int v=0; v<inst.speed.getSize(); v++){
		for(int i=1; i<=J; i++){
            //Fixing values of last index of previous fixed block (if is not first iteration)
            if(tS>1){
                if(hasEnteringArc1st[v][i][tS-1]){
                    w[v][i][tS-1].setBounds(round(wValue[v][i][tS-1]), round(wValue[v][i][tS-1]));
                    #ifdef WaitAfterOperate
                    wB[v][i][tS-1].setBounds(round(wBValue[v][i][tS-1]), round(wBValue[v][i][tS-1]));
                    #endif
                    #ifndef WaitAfterOperate
                    oB[v][i][tS-1].setBounds(round(oBValue[v][i][tS-1]), round(oBValue[v][i][tS-1]));
                    #endif
                }
            }
			for(int t=1; t<=tF; t++){
				if(hasEnteringArc1st[v][i][t]){
					if(t>=tS){ //Fixing values only of the current block                        
                        z[v][i][t].setBounds(round(zValue[v][i][t]), round(zValue[v][i][t]));
                        if(t<tF){
                            w[v][i][t].setBounds(round(wValue[v][i][t]), round(wValue[v][i][t]));
                            #ifdef WaitAfterOperate
                            wB[v][i][t].setBounds(round(wBValue[v][i][t]), round(wBValue[v][i][t]));
                            #endif
                            #ifndef WaitAfterOperate
                            oB[v][i][t].setBounds(round(oBValue[v][i][t]), round(oBValue[v][i][t]));
                            #endif
                        }
                        #ifndef WaitAfterOperate
                        oA[v][i][t].setBounds(round(oAValue[v][i][t]), round(oAValue[v][i][t]));
                        #endif
					}
					for(int j=i+1;j<=J;j++){ //Sink arc are not fixed
						if(hasArc[v][i][j][t]){
                            int t2 = t + inst.travelTime[v][i-1][j-1];
                            if (t2>= tS && t2 <= tF) 
                                x[v][i][j][t].setBounds(round(xValue[v][i][j][t]), round(xValue[v][i][j][t]));
						}
					}
				}
			}
		}
	}
}

void Model::decreaseEndBlock (IloEnv& env, Instance inst, const double& nIntervals, const int& tS, const int& tF){	
	int i, j,t,v,a;
	int timePerInterval = inst.t/nIntervals; //Number of time periods in each interval
	int J = inst.numTotalPorts;
	int T = inst.t+1;
	double intPart;	
	int V = inst.speed.getSize(); //# of vessels
    int N = J + 2;
	
    ///Objective function
	IloExpr expr_obj_current = cplex.getObjective().getExpr(); //Gets the current objective function
	IloExpr expr_obj_cost(env);
	IloExpr expr_obj_revenue(env);
	for(v=0;v<V;v++){
		for(i=1;i<=J;i++){
			for(t=0;t<=tF;t++){
				for(j=i+1;j<N;j++){
					if(hasArc[v][i][j][t]==1){
						if(j<=J){ 								//If j is a port
							int t2 = t + inst.travelTime[v][i-1][j-1]; 
							if(t2>=tS && t2<=tF)     //If was not present in the model
								expr_obj_cost += arcCost[v][i][j][t]*x[v][i][j][t];
						}else{									//If j is sink node
							if(t>=tS)
								expr_obj_cost += arcCost[v][i][j][t]*x[v][i][j][t];
						}
					}
				}
				if(t>=tS && hasEnteringArc1st[v][i][t]){
					expr_obj_revenue += inst.r_jt[i-1][t-1]*f[v][i][t];						//1st term
					expr_obj_cost += (t-1)*inst.attemptCost*z[v][i][t];						//3rd term				
				}
			}
		}
	}
	//Port-time iterator
	for(i=1;i<=J;i++){
		for(t=tS;t<=tF;t++){
			expr_obj_cost += inst.p_jt[i-1][t-1]*alpha[i][t];								//4rd term
			#ifndef NBetas
            expr_obj_cost += 1000*beta[i][t];										//Auxiliary variables
            #endif
		}
		#ifndef NKnapsackInequalities
		int it,it2=0,k,l;
		IloExpr expr_kP1_LHS(env), expr_kP2_LHS(env), expr_kD1_LHS(env), expr_kD2_LHS(env);
		for(it=0;it <((T-3)*2 - (T-1-tF));it++){	//'Revise' existing inequalities (x variables) and add the new until [tF-1, tF]
			double kP1_RHS=0, kP2_RHS=0, kD1_RHS=0, kD2_RHS=0, sum_alphaMax=0, alphaUB=0;
			expr_kP1_LHS.clear();
			expr_kP2_LHS.clear();
			expr_kD1_LHS.clear();
			expr_kD2_LHS.clear();
			//Definining the size of set T =[k,l] - CAUTION: It is inverse of the paper of Agra et al (2013)
			k=1;
			l=tF;
			if(it<tF-2){
				l = it+2;
				it2++;
			}else if(it==tF-2){
				it = T-3;
				k = it2+4-l;
				it2++;
			}else{
				k = it2+4-l;
				it2++;
			}            
			if( (it >= tS-2 && it < tF-2) || (it > T-3 + tS-2) ){ //valid inequality is created from scratch - If it is in interval [1...[tS-2..tF-1]], or If it is in interval [[2..tS-1]...tF]
				for(v=0;v<V;v++){
					#ifdef WaitAfterOperate
						expr_kP1_LHS += wB[v][i][k];
					if(k>1 || hasEnteringArc1st[v][i][k-1]==1)      //TODO verificar se não deve ser &&
						expr_kD1_LHS += w[v][i][k-1] + wB[v][i][k-1];
					#endif
					#ifndef WaitAfterOperate
						expr_kP1_LHS += oB[v][i][k];

					if(k>1 || hasEnteringArc1st[v][i][k-1]==1)
						expr_kD1_LHS += w[v][i][k-1] + oB[v][i][k-1];
					#endif
					#ifndef WaitAfterOperate
					if(k>1 || hasEnteringArc1st[v][i][k-1]==1)
						expr_kD2_LHS += oB[v][i][k-1];
					#endif
					for(t=k;t<=l;t++){                       
						expr_kP2_LHS += z[v][i][t];	//It is used for both loading and discharging ports
						#ifndef WaitAfterOperate
						expr_kD2_LHS += oA[v][i][t];
						#endif
						#ifdef WaitAfterOperate
						expr_kD2_LHS += z[v][i][t];
						#endif
						for(j=0;j<N;j++){						
							if(j==0 && inst.initialPort[v]+1 == i && t==inst.firstTimeAv[v]+1) 	//Source arc
								expr_kD1_LHS += x[v][j][i][0];								
							else if (j == N-1 && hasArc[v][i][j][t] == 1)						//Sink arc
								expr_kP1_LHS += x[v][i][j][t];
							else if (j > 0 && j < N-1){											//"Normal" arcs
								if(i != j){
									if(hasArc[v][i][j][t] == 1 && t+inst.travelTime[v][i-1][j-1] <= tF)  //If arc exists and arrives at port j in the integer or relaxed block
										expr_kP1_LHS += x[v][i][j][t];
									if(t - inst.travelTime[v][j-1][i-1] >= 0){					//Only if it is possible an entering arc due the travel time
										if(hasArc[v][j][i][t-inst.travelTime[v][j-1][i-1]] == 1)
											expr_kD1_LHS += x[v][j][i][t-inst.travelTime[v][j-1][i-1]];
									}
								}
							}
						}
					}
				}
			
				for(t=k;t<=l;t++){
					kP1_RHS += inst.d_jt[i-1][t-1];
					kD1_RHS += inst.d_jt[i-1][t-1];
				}
				kP1_RHS += -inst.sMax_jt[i-1][0];
				kP2_RHS = ceil(kP1_RHS/inst.f_max_jt[i-1][0]);
				kP1_RHS = ceil(kP1_RHS/inst.maxCapacity);
				
				if(k==1)
					kD1_RHS += -inst.s_j0[i-1] + inst.sMin_jt[i-1][0];
				else
					kD1_RHS += -inst.sMax_jt[i-1][k-1] + inst.sMin_jt[i-1][0];
				kD2_RHS = ceil(kD1_RHS/inst.f_max_jt[i-1][0]);
				kD1_RHS = ceil(kD1_RHS/inst.maxCapacity);
				
				stringstream ss, ss1, ss2;
				if(inst.typePort[i-1] == 0){
					ss << "knpasackP1_" << i << "," << it;
					knapsack_P_1[i][it] = IloRange(env, kP1_RHS, expr_kP1_LHS, IloInfinity, ss.str().c_str());
					ss1 << "knapsackP2_" << i << "," << it;
					knapsack_P_2[i][it] = IloRange(env, kP2_RHS, expr_kP2_LHS, IloInfinity, ss1.str().c_str());
					model.add(knapsack_P_1[i][it]);
					model.add(knapsack_P_2[i][it]);
				}else{
					ss << "knpasackD1_" << i << "," << it;
					knapsack_D_1[i][it] = IloRange(env, kD1_RHS, expr_kD1_LHS, IloInfinity, ss.str().c_str());
					ss1 << "knpasackD2_" << i << "," << it;
					knapsack_D_2[i][it] = IloRange(env, kD1_RHS, expr_kD2_LHS, IloInfinity, ss1.str().c_str());
					model.add(knapsack_D_1[i][it]);
					model.add(knapsack_D_2[i][it]);
					#ifndef WaitAfterOperate
					ss2 << "knpasackD3_" << i << "," << it;
					knapsack_D_3[i][it] = IloRange(env, kD2_RHS, expr_kP2_LHS, IloInfinity, ss2.str().c_str());
					model.add(knapsack_D_3[i]);
					#endif
				}
			}else if (it >= 0 && it < tS-2 ){		//Only modify the previously added valid inequalities of type [1..j] (only needed add x variables in the case of loading ports)
				//Get the current expr values                 
				if (inst.typePort[i-1] == 0){
                    expr_kP1_LHS = knapsack_P_1[i][it].getExpr();                    
                    #ifndef WaitAfterOperate 
                    expr_kP2_LHS = knapsack_D_3[i][it].getExpr();
                    #endif 
                }else{     //Discharging ports               
                    expr_kD1_LHS = knapsack_D_1[i][it].getExpr();                     
                    expr_kD2_LHS = knapsack_D_2[i][it].getExpr();                    
                }

				for(v=0;v<V;v++){										
					for(t=k;t<=l;t++){						
						#ifndef WaitAfterOperate
						if(t>=tS && t<= tF)			//Variables of 'extended' interval
							expr_kD2_LHS += oA[v][i][t];
							expr_kP2_LHS += z[v][i][t];
						#endif                        
						#ifdef WaitAfterOperate
						if(t>=tS && t<= tF)			//Variables of 'extended' interval
							expr_kD2_LHS += z[v][i][t];
						#endif                        
						for(j=0;j<N;j++){													
							if (j == N-1 && hasArc[v][i][j][t] == 1 && t>=tS)					//Sink arc departing after the last 'l'
								expr_kP1_LHS += x[v][i][j][t];
							else if (j > 0 && j < N-1){											//"Normal" arcs	
								if(i != j){
									int t2 = t+inst.travelTime[v][i-1][j-1];
									if(hasArc[v][i][j][t] == 1 && (t2>=tS && t2<=tF) )  //If arc exists, was not added in the previous iteration and arrives at port j in the integer or relaxed block
										expr_kP1_LHS += x[v][i][j][t];
									if (t >= tS && t <= tF){									//If entering arc was not added in the previous iteration
										if(t - inst.travelTime[v][j-1][i-1] >= 0){					//Only if it is possible an entering arc due the travel time
											if(hasArc[v][j][i][t-inst.travelTime[v][j-1][i-1]] == 1)
												expr_kD1_LHS += x[v][j][i][t-inst.travelTime[v][j-1][i-1]];
										}
									}
								}
							}
						}
					}
				}
			
				for(t=k;t<=l;t++){ //Calculate kD1_RHS from scratch
					kD1_RHS += inst.d_jt[i-1][t-1];
				}				
				kP2_RHS = ceil(kP1_RHS/inst.f_max_jt[i-1][0]);
							
				if(k==1)
					kD1_RHS += -inst.s_j0[i-1] + inst.sMin_jt[i-1][0] - alphaUB;
				else
					kD1_RHS += -inst.sMax_jt[i-1][k-1] + inst.sMin_jt[i-1][0] - alphaUB;
				kD2_RHS = ceil(kD1_RHS/inst.f_max_jt[i-1][0]);
				kD1_RHS = ceil(kD1_RHS/inst.maxCapacity);
				
				if(inst.typePort[i-1] == 0){					
					knapsack_P_1[i][it].setExpr(expr_kP1_LHS);					
				}else{					
					knapsack_D_1[i][it].setExpr(expr_kD1_LHS);
					knapsack_D_1[i][it].setLB(kD1_RHS);					
					knapsack_D_2[i][it].setExpr(expr_kD2_LHS);
					knapsack_D_2[i][it].setLB(kD1_RHS);
					#ifndef WaitAfterOperate					
					knapsack_D_3[i][it].setExpr(expr_kP2_LHS);
					knapsack_D_3[i][it].setLB(kD2_RHS);					
					#endif
				}				
			}
		}
		#endif
	}
    
	//New objective
	obj.setExpr(expr_obj_current + expr_obj_cost - expr_obj_revenue);
    IloExpr expr_cumSlack(env), expr_berth(env), expr_invBalancePort(env);
	for(i=1;i<N-1;i++){        
		expr_cumSlack.clear();
        expr_cumSlack = cumSlack[i].getExpr();
		for(t=tS;t<=tF;t++){
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
			#ifndef NBetas
            expr_invBalancePort += -beta[i][t];
            #endif
			ss2 << "invBalancePort_(" << i << "," << t << ")";
			portInventory[i][t] = IloRange(env, inst.delta[i-1]*inst.d_jt[i-1][t-1],
				sP[i][t]-sP[i][t-1]-inst.delta[i-1]*expr_invBalancePort,
				inst.delta[i-1]*inst.d_jt[i-1][t-1], ss2.str().c_str());
			model.add(portInventory[i][t]);
		}		
		cumSlack[i].setExpr(expr_cumSlack);
	}
	expr_cumSlack.end();
	expr_berth.end();
	expr_invBalancePort.end();
    
    IloExpr expr_sinkFlow(env), expr_1stLevel(env), expr_2ndLevel(env), expr_1stFlow(env), expr_2ndFlow(env);
	for(v=0;v<V;v++){
        expr_sinkFlow.clear();
        expr_sinkFlow = sinkNodeBalance[v].getExpr();
		for(i=0;i<N;i++){            
			if(i>0 && i <= J){ //Only considering ports
				//Updating flows and limits of last time period of previous interval                 
                expr_1stLevel.clear();
                expr_2ndLevel.clear();
                expr_1stFlow.clear();
                expr_2ndFlow.clear();
                
                if(hasEnteringArc1st[v][i][tS-1]){ //Only necessary if there is entering arc in the node					
					expr_1stLevel = firstLevelBalance[v][i][tS-1].getExpr();
					expr_2ndLevel = secondLevelBalance[v][i][tS-1].getExpr();
					expr_1stFlow = firstLevelFlow[v][i][tS-1].getExpr();
					expr_2ndFlow = secondLevelFlow[v][i][tS-1].getExpr();
					
					expr_1stLevel += - w[v][i][tS-1];
					expr_1stFlow += -fW[v][i][tS-1];
					
					#ifndef WaitAfterOperate                
					expr_2ndLevel += oB[v][i][tS-1];
					expr_2ndFlow += fOB[v][i][tS-1];
					#endif
					
					#ifdef WaitAfterOperate
					expr_2ndLevel += wB[v][i][tS-1];
					expr_2ndFlow += fWB[v][i][tS-1];
					#endif
					
					firstLevelBalance[v][i][tS-1].setExpr(expr_1stLevel);
					firstLevelFlow[v][i][tS-1].setExpr(expr_1stFlow);
					secondLevelBalance[v][i][tS-1].setExpr(expr_2ndLevel);
					secondLevelFlow[v][i][tS-1].setExpr(expr_2ndFlow);
				}
                
                //Add constraints on the flow variables of last time period of previous block.
                stringstream ss8,ss9,ss10;
                #ifndef WaitAfterOperate
                ss8 << "flowLimitOB_"<<v<<","<<i<<","<<tS-1;
                flowCapacityOB[v][i][tS-1] = IloRange(env, -IloInfinity, fOB[v][i][tS-1]-inst.q_v[v]*oB[v][i][tS-1], 0, ss8.str().c_str());
                model.add(flowCapacityOB[v][i][tS-1]);                            
                #endif
                ss9 << "flowLimitW_"<<v<<","<<i<<","<<tS-1;
                flowCapacityW[v][i][tS-1] = IloRange(env, -IloInfinity, fW[v][i][tS-1]-inst.q_v[v]*w[v][i][tS-1], 0, ss9.str().c_str());
                model.add(flowCapacityW[v][i][tS-1]);
                #ifdef WaitAfterOperate
                ss10 << "flowLimitWB_"<<v<<","<<i<<","<<tS-1;
                flowCapacityWB[v][i][tS-1] = IloRange(env, -IloInfinity, fWB[v][i][tS-1]-inst.q_v[v]*wB[v][i][tS-1], 0, ss10.str().c_str());
                model.add(flowCapacityWB[v][i][tS-1]);
                #endif
                //end updating                
                for(j=0;j<N;j++){
					if(j>0){ //Considering ports and sink node
						for(t=1;t<=tF;t++){
                            stringstream ss, ss1;
                            if (j<=J){				//When j is a port
                                int t2 = t + inst.travelTime[v][i-1][j-1];
                                if(t2>=tS && t2<=tF){ // Arc in the new added block or in the intersection with new block and previous interval
                                    if(inst.typePort[i-1] == 0 && inst.typePort[j-1] == 1){
                                        ss << "travelAtCap_(" << i << "," << j << ")_" << t << "," << v;
                                        travelAtCapacity[v][i][j][t].setExpr(-fX[v][i][j][t] + inst.q_v[v]*x[v][i][j][t]);
                                        travelAtCapacity[v][i][j][t].setName(ss.str().c_str());
                                        model.add(travelAtCapacity[v][i][j][t]);
                                    }else if (inst.typePort[i-1] == 1 && inst.typePort[j-1] == 0){
                                        ss1 << "travelEmpty_(" << i << "," << j << ")_" << t << "," << v;
                                        travelEmpty[v][i][j][t].setExpr(inst.q_v[v]*x[v][i][j][t] + fX[v][i][j][t]);
                                        travelEmpty[v][i][j][t].setName(ss1.str().c_str());
                                        model.add(travelEmpty[v][i][j][t]);
                                    }
                                }
							}else{ // when j is the sink node
								if(t >= tS){ //Only for nodes in the current added block
                                    if(inst.typePort[i-1] == 0){
                                        ss << "travelAtCap_(" << i << ",snk)_" << t << "," << v;
                                        travelAtCapacity[v][i][j][t].setExpr(-fX[v][i][j][t] + inst.q_v[v]*x[v][i][j][t]);	
                                        travelAtCapacity[v][i][j][t].setName(ss.str().c_str());
                                        model.add(travelAtCapacity[v][i][j][t]);
                                    }else if (inst.typePort[i-1] == 1){ 
                                        ss1 << "travelEmpty_(" << i << ",snk)_" << t << "," << v;
                                        travelEmpty[v][i][j][t].setExpr(inst.q_v[v]*x[v][i][j][t] + fX[v][i][j][t]);
                                        travelEmpty[v][i][j][t].setName(ss1.str().c_str());
                                        model.add(travelEmpty[v][i][j][t]);
                                    }
                                }
							}
						}
					}
				}
				for(t=1;t<=tF;t++){
                    if(t<tS){   ///Needed for update outgoing flow balance arcs (just 2nd level)
                        expr_2ndLevel = secondLevelBalance[v][i][t].getExpr();
                        expr_2ndFlow = secondLevelFlow[v][i][t].getExpr();
                        for(j=1;j<=J;j++){
                            if (hasArc[v][i][j][t] == 1){
                                if (t + inst.travelTime[v][i-1][j-1] >= tS && t + inst.travelTime[v][i-1][j-1] <= tF){ //If arc departs from previos interval (t<tS) and arrive in the current interval
                                    expr_2ndLevel += x[v][i][j][t];
                                    expr_2ndFlow += fX[v][i][j][t];
                                }
                            }
                        }
                        //Updating previous flow 2nd level
                        secondLevelBalance[v][i][t].setExpr(expr_2ndLevel);
                        secondLevelFlow[v][i][t].setExpr(expr_2ndFlow);
                    }else{ /// if t>= tS : Only considering the time-periods of the added interval
                        stringstream ss, ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8, ss9, ss10;
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
                        for(j=0;j<N;j++){
                            if(j>0 && j<N-1){ //No consider sink arc (first level balance)
                                if (t - inst.travelTime[v][j-1][i-1] >= 0){ //If it is possible to exist an arc from j to i
                                    if(hasArc[v][j][i][t-inst.travelTime[v][j-1][i-1]] == 1){ //If the arc exists
                                        expr_1stLevel += x[v][j][i][t-inst.travelTime[v][j-1][i-1]];
                                        expr_1stFlow += fX[v][j][i][t-inst.travelTime[v][j-1][i-1]];
                                    }
                                }
                            }
                            if(j>0){ //No consider source arc (second level balance)
                                if (hasArc[v][i][j][t] == 1){
                                    if(j == N-1){ //If j is the sink node, add to expr
                                        expr_2ndLevel += x[v][i][j][t];
                                        expr_2ndFlow += fX[v][i][j][t];
                                    }else if (t + inst.travelTime[v][i-1][j-1] <= tF){ //If j is a port, it is necessary that the arrival is in the model
                                        expr_2ndLevel += x[v][i][j][t];
                                        expr_2ndFlow += fX[v][i][j][t];
                                    }
                                }
                            }
                        }
                        IloExpr expr_link;
                        #ifndef WaitAfterOperate 
                        if (t < tF && hasEnteringArc1st[v][i][t-1]==0){ //Not last time period nor entering arc in the previous time period
                            expr_1stLevel += - w[v][i][t] - oA[v][i][t];
                            expr_2ndLevel += -oA[v][i][t] + oB[v][i][t];
                            expr_1stFlow+= -fW[v][i][t] - fOA[v][i][t];
                            expr_link = oA[v][i][t] - z[v][i][t];
                            expr_2ndFlow += -fOA[v][i][t] - inst.delta[i-1]*f[v][i][t] + fOB[v][i][t];
                        }else if (t==tF){ //Last time period
                            if (hasEnteringArc1st[v][i][t-1]==1){
                                expr_1stLevel += w[v][i][t-1] - oA[v][i][t];
                                expr_2ndLevel += -oA[v][i][t] -oB[v][i][t-1];
                                expr_link = oA[v][i][t] + oB[v][i][t-1] - z[v][i][t];   
                                expr_1stFlow += fW[v][i][t-1] - fOA[v][i][t];
                                expr_2ndFlow += -fOA[v][i][t] - fOB[v][i][t-1] - inst.delta[i-1]*f[v][i][t];
                            }else{
                                expr_1stLevel += - oA[v][i][t];
                                expr_2ndLevel += -oA[v][i][t];
                                expr_link = oA[v][i][t] - z[v][i][t];   
                                expr_1stFlow += - fOA[v][i][t];
                                expr_2ndFlow += -fOA[v][i][t] - inst.delta[i-1]*f[v][i][t];
                            }
                        }else{ //t<tF and hasEnteringArc1st = 1
                            expr_1stLevel += w[v][i][t-1] - w[v][i][t] - oA[v][i][t];
                            expr_2ndLevel += - oA[v][i][t] - oB[v][i][t-1] + oB[v][i][t];
                            expr_link = oA[v][i][t] + oB[v][i][t-1] - z[v][i][t];
                            expr_1stFlow += fW[v][i][t-1] - fW[v][i][t] - fOA[v][i][t];
                            expr_2ndFlow += -fOA[v][i][t] - fOB[v][i][t-1] - inst.delta[i-1]*f[v][i][t] + fOB[v][i][t];	
                        }
                        #endif
                        #ifdef WaitAfterOperate
                        if (t < tF && hasEnteringArc1st[v][i][t-1]==0){ //Not last time period nor entering arc in the previous time period
                            expr_1stLevel += - w[v][i][t] - z[v][i][t];
                            expr_2ndLevel += -z[v][i][t] + wB[v][i][t];	
                            expr_1stFlow+= -fW[v][i][t] - fOA[v][i][t];	
                            expr_2ndFlow += -fOA[v][i][t] - inst.delta[i-1]*f[v][i][t] + fWB[v][i][t];
                        }else if (t==tF){ //Last time period
                            if (hasEnteringArc1st[v][i][t-1]==1){
                                expr_1stLevel += w[v][i][t-1] - z[v][i][t] + wB[v][i][t-1];                                
                                expr_1stFlow += fW[v][i][t-1] - fOA[v][i][t] + fWB[v][i][t-1];                                
                            }else{
                                expr_1stLevel += - z[v][i][t];                                
                                expr_1stFlow += - fOA[v][i][t];                                
                            }
                            expr_2ndLevel += -z[v][i][t];
                            expr_2ndFlow += -fOA[v][i][t] - inst.delta[i-1]*f[v][i][t];
                        }else{ //t<tF and hasEnteringArc1st = 1
                            expr_1stLevel += w[v][i][t-1] - w[v][i][t] - z[v][i][t] + wB[v][i][t-1];
                            expr_2ndLevel += - z[v][i][t] + wB[v][i][t];	
                            expr_1stFlow += fW[v][i][t-1] - fW[v][i][t] - fOA[v][i][t] + fWB[v][i][t-1];
                            expr_2ndFlow += -fOA[v][i][t] - inst.delta[i-1]*f[v][i][t] + fWB[v][i][t];
                        }
                        #endif
                        
                        firstLevelBalance[v][i][t].setExpr(expr_1stLevel);
                        secondLevelBalance[v][i][t].setExpr(expr_2ndLevel);
                        #ifndef WaitAfterOperate
                        linkBalance[v][i][t].setExpr(expr_link);
                        #endif
                        firstLevelFlow[v][i][t].setExpr(expr_1stFlow);
                        secondLevelFlow[v][i][t].setExpr(expr_2ndFlow);
                        
                        firstLevelBalance[v][i][t].setName(ss.str().c_str());
                        model.add(firstLevelBalance[v][i][t]);
                        
                        secondLevelBalance[v][i][t].setName(ss1.str().c_str());
                        model.add(secondLevelBalance[v][i][t]);
                        
                        #ifndef WaitAfterOperate
                        linkBalance[v][i][t].setName(ss2.str().c_str());
                        model.add(linkBalance[v][i][t]);
                        #endif
                        
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
                        #ifndef WaitAfterOperate
                        flowCapacityOA[v][i][t] = IloRange(env, -IloInfinity, fOA[v][i][t]-inst.q_v[v]*oA[v][i][t], 0, ss7.str().c_str());
                        #endif
                        #ifdef WaitAfterOperate
                        flowCapacityOA[v][i][t] = IloRange(env, -IloInfinity, fOA[v][i][t]-inst.q_v[v]*z[v][i][t], 0, ss7.str().c_str());
                        #endif
                        model.add(flowCapacityOA[v][i][t]);
                        
                        if(t<tF){ //Constraints with no last time index
                            #ifndef WaitAfterOperate
                            ss8 << "flowLimitOB_"<<v<<","<<i<<","<<t;
                            flowCapacityOB[v][i][t] = IloRange(env, -IloInfinity, fOB[v][i][t]-inst.q_v[v]*oB[v][i][t], 0, ss8.str().c_str());
                            model.add(flowCapacityOB[v][i][t]);                            
                            #endif
                            ss9 << "flowLimitW_"<<v<<","<<i<<","<<t;
                            flowCapacityW[v][i][t] = IloRange(env, -IloInfinity, fW[v][i][t]-inst.q_v[v]*w[v][i][t], 0, ss9.str().c_str());							
                            model.add(flowCapacityW[v][i][t]);
                            #ifdef WaitAfterOperate
                            ss10 << "flowLimitWB_"<<v<<","<<i<<","<<t;
                            flowCapacityWB[v][i][t] = IloRange(env, -IloInfinity, fWB[v][i][t]-inst.q_v[v]*wB[v][i][t], 0, ss10.str().c_str());
                            model.add(flowCapacityWB[v][i][t]);
                            #endif
                        }
                    }
				}
                for(j=1;j<N;j++){
                    if(i != j && i < N-1){ //There is no departing arc from sink node
                        for(t=1;t<=tF;t++){
                            stringstream ss;
                            ss << "flowLimitX_" << v << "," << i << "," << j << "," << t;                            
                            if(j == N-1 && t>=tS){ //If j is sink node and t belong to the current added block
                                flowCapacityX[v][i][j][t] = IloRange(env, -IloInfinity, fX[v][i][j][t] - inst.q_v[v]*x[v][i][j][t], 0, ss.str().c_str());
                                model.add(flowCapacityX[v][i][j][t]);
                            }else if (j != N-1){                                 
                                int t2 = t + inst.travelTime[v][i-1][j-1];                                
                                if (t2>=tS && t2<=tF){ //If fx variable cross the intersection or is in the added block
                                    flowCapacityX[v][i][j][t] = IloRange(env, -IloInfinity, fX[v][i][j][t] - inst.q_v[v]*x[v][i][j][t], 0, ss.str().c_str());
                                    model.add(flowCapacityX[v][i][j][t]);
                                }
                            }
                        }
                    }
                }
            }
		}
		sinkNodeBalance[v].setExpr(expr_sinkFlow);
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
void Model::reIntegralize(IloEnv& env, Instance inst, const int& tS, const int& tF){
	int i,j,t,v;
	int J = inst.numTotalPorts;	
	int V = inst.speed.getSize(); //# of vessels
    int N = J + 2;
	
	for(v=0;v<V;v++){
		for(i=1; i<=J; i++){
            //Integralizing the last index of previous block
            model.remove(convertW[v][i][tS-1]);
            #ifndef WaitAfterOperate
            model.remove(convertOB[v][i][tS-1]);
            #endif
            #ifdef WaitAfterOperate
            model.remove(convertWB[v][i][tS-1]);
            #endif
            
			for(t=1; t<=tF; t++){
				if(t>=tS){
					model.remove(convertZ[v][i][t]);
					#ifndef WaitAfterOperate
					model.remove(convertOA[v][i][t]);
                    #endif
                    if(t<tF){
                        model.remove(convertW[v][i][t]);
                        #ifndef WaitAfterOperate
                        model.remove(convertOB[v][i][t]);
                        #endif
                        #ifdef WaitAfterOperate
                        model.remove(convertWB[v][i][t]);
                        #endif
                    }
				}
				for(j=i+1;j<=N-1;j++){					
					if (j < N-1){   		//If arrives at a port AND
						int t2 = t + inst.travelTime[v][i-1][j-1];
						if(t2 >= tS && t2 <= tF){ //Traveling arc is in the integer block OR crosses previos integer block and the integer block
							model.remove(convertX[v][i][j][t]);							
						}
					}else if (t>= tS && j == N-1){ //Sink arc						
						model.remove(convertX[v][i][j][t]);
					}
				}
			}
		}
	}
}

//~ void Model::fixAllSolution(IloEnv& env,const Instance& inst){
	//~ if((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){				
		//~ getSolVals(env, inst);			
		//~ //Fix solution		
		//~ for(int v=0; v<inst.speed.getSize(); v++){
			//~ #ifndef NFixZvar
			//~ for(int j=0; j<inst.numTotalPorts; j++){							
				//~ for(int t=0; t<inst.t; t++){								
					//~ z[v][j][t].setBounds(round(zValue[v][j][t]), round(zValue[v][j][t]));					
				//~ }
			//~ }
			//~ #endif
			//~ //Fixing X variables
			//~ for(int a=0; a<x[v].getSize(); a++){				
				//~ x[v][a].setBounds(round(xValue[v][a]), round(xValue[v][a]));			   
			//~ }						
			//~ #ifndef NOperateAndGo2
			//~ for(int t=0; t<inst.t; t++){
				//~ w[v][t].setBounds(round(wValue[v][t]), round(wValue[v][t]));
			//~ }
			//~ #endif
		//~ }		
	//~ }else{//If the method return a infeasible solution
		//~ //TODO
		//~ cout << "No  feasile solution" << endl;
		//~ exit(1);
	//~ }
//~ }

//~ void Model::unFixInterval(Instance inst, const int& tS, const int& tF){
	//~ //Unfixing z
	//~ for(int v=0; v<inst.speed.getSize(); v++){
		//~ #ifndef NFixZvar
		//~ for(int j=0; j<inst.numTotalPorts; j++){							
			//~ for(int t=tS; t<tF; t++){								
				//~ z[v][j][t].setBounds(0, 1);					
			//~ }
		//~ }
		//~ #endif
		//~ //Unfixing X variables
		//~ for(int a=0; a<x[v].getSize(); a++){
			//~ int t1,t2, arcType;		
			//~ arcType = inst.travelArcType(inst.arcs[v][a],t1,t2);	
			//~ if( (t2 >= tS && t2 < tF) ||  // Arc times are between relaxed times and before 				
				//~ ((arcType == 3 || arcType == 4) && (t1 >= tS && t1 < tF)) ){ //Is a sink arc in the interval 
				//~ x[v][a].setBounds(0, 1);
		   //~ }
		//~ }		
	//~ }
//~ }
/*Pre-req: Fix All solution*/
//~ void Model::improvementPhase(IloEnv& env, Instance inst, const double& mIntervals, const double& timePerIter, const double& gap, const double& overlap, Timer<chrono::milliseconds>& timer_cplex,float& opt_time,
//~ const double& timeLimit, float& elapsed_time, double& incumbent){	
	//~ double prevObj = 1.0e+10;
	//~ double objValue = incumbent;
	//~ int i;			
	//~ Timer<chrono::milliseconds> timer_LS;
	//~ timer_LS.start();
	//~ cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
	//~ if (gap > 1e-04)
		//~ cplex.setParam(IloCplex::EpGap, gap/100);
	//~ int tS, tF;
	//~ while((prevObj - objValue > 0.0001) && elapsed_time/1000 < timeLimit){				
		//~ for(i=1;i<=ceil(mIntervals*(1+overlap/100));i++){
			//~ if(i==1)
				//~ tS = inst.t/mIntervals*(i-1);							
			//~ else
				//~ tS = inst.t/mIntervals*(i-1)*(1-overlap/100);						
			//~ tF = min(tS+inst.t/mIntervals, (double)inst.t);
			//~ cout << "Unfixing " << tS << "..." << tF << endl;
			//~ unFixInterval(inst, tS, tF);
			//~ 
			//~ cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
			//~ timer_cplex.start();
			//~ cplex.solve();
			//~ opt_time += timer_cplex.total();
			//~ elapsed_time += timer_LS.total();			
			//~ incumbent = cplex.getObjValue();
			//~ if(elapsed_time/1000 >= timeLimit){
				//~ cout << "Reached LS time limit " << timeLimit << ": " << elapsed_time/1000 << endl;
				//~ break;
			//~ }
			//~ timer_LS.start();						
			//~ cout << "Objective: " << incumbent << endl;
			//~ // //~ if (i < ceil(mIntervals*(1+overlap/100))){
				//~ cout << "Re-fixing " << tS << "..." << tF;
				//~ fixSolution(env, inst, tS, tF,0);	 //Does not fix the last interval either the sink arcs
				//~ cout << ". Done! " << endl;
			//~ // //~ }			
		//~ }
		//~ prevObj = objValue;
		//~ objValue = incumbent;		
		//~ if(elapsed_time/1000 >= timeLimit){
			//~ cout << "Reached LS time limit " << timeLimit << ": " << elapsed_time/1000 << endl;
			//~ break;
		//~ }				
	//~ }
//~ }

//~ void Model::getSolVals(IloEnv& env, const Instance& inst){
	//~ if((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){				
		//~ for(int v=0;v<inst.speed.getSize();v++){
			//~ xValue[v] = IloNumArray(env, xValue.getSize());			
			//~ cplex.getValues(x[v], xValue[v]);			
			//~ for(int j=0; j<inst.numTotalPorts; j++){				
				//~ zValue[v][j] = IloNumArray(env, zValue[v][j].getSize());				
				//~ cplex.getValues(z[v][j], zValue[v][j]);
			//~ }
			//~ #ifndef NOperateAndGo2
			//~ //wValue[v] = IloNumArray(env, wValue.getSize());			
			//~ // cplex.getValues(w[v], wValue[v]);			
			//~ #endif
		//~ }
	//~ }else{
		//~ cout << "Impossible to get feasible solution values!" << endl;
		//~ exit(1);
	//~ }
//~ }
void Model::getSolValsW(IloEnv& env, Instance inst, const int& tS, const int& tF){
	if((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){
		int v, i, j, t;
        int V = inst.speed.getSize();
        int J = inst.numTotalPorts;
        int N = J+2;
        for(v=0;v<V;v++){
            for(i=1;i<=J;i++){ //Not necessary consider source node either sink node
                for (j=i+1;j<=J;j++){ //Consider only ports as the sink arcs are not fixed
                    //~ xValue[v][i][j] = IloNumArray(env, inst.t+1); //Reset values
                    for(t=1;t<=tF;t++){
                        if(hasArc[v][i][j][t]==1){
                            int t2 = t + inst.travelTime[v][i-1][j-1];
                            if (t2>=tS && t2<=tF) //If is traveling arc is in the integer block or crossing a fixed block and integer block
                                xValue[v][i][j][t] = cplex.getValue(x[v][i][j][t]);                            
                        }
                    }
                }
                //Get the values of last index of previous block
                if(tS > 1){ //If is not the first iteration 
                    if(hasEnteringArc1st[v][i][tS-1] == 1){
						#ifndef WaitAfterOperate                        
                        oBValue[v][i][tS-1] = cplex.getValue(oB[v][i][tS-1]);                        
						#endif
						#ifdef WaitAfterOperate
						wBValue[v][i][tS-1] = cplex.getValue(wB[v][i][tS-1]);                        
						#endif
                        wValue[v][i][tS-1] = cplex.getValue(w[v][i][tS-1]);                        
					}
                }
                for(t=tS;t<=tF;t++){
                    if(hasEnteringArc1st[v][i][t] == 1){
						zValue[v][i][t] = cplex.getValue(z[v][i][t]);
                        //~ cout << "z_["<<v<<"]["<<i<<"]["<<t<<"] = " << zValue[v][i][t] << endl;
						
						#ifndef WaitAfterOperate
						oAValue[v][i][t] = cplex.getValue(oA[v][i][t]);
                        //~ cout << "oA_["<<v<<"]["<<i<<"]["<<t<<"] = " << oAValue[v][i][t] << endl;
                        #endif
                        
                        if(t<tF){
                            #ifndef WaitAfterOperate
                            oBValue[v][i][t] = cplex.getValue(oB[v][i][t]);
                            //~ cout << "oB_["<<v<<"]["<<i<<"]["<<t<<"] = " << oBValue[v][i][t] << endl;
                            #endif
                            #ifdef WaitAfterOperate
                            wBValue[v][i][t] = cplex.getValue(wB[v][i][t]);
                            //~ cout << "wB_["<<v<<"]["<<i<<"]["<<t<<"] = " << wBValue[v][i][t] << endl;
                            #endif                            
                            wValue[v][i][t] = cplex.getValue(w[v][i][t]);
                            //~ cout << "w_["<<v<<"]["<<i<<"]["<<t<<"] = " << wValue[v][i][t] << endl;
						}
					}
                }       
            }
        }
	}else{
		cout << "Impossible to get feasible solution values!" << endl;
		exit(1);
	}
}

//~ void Model::regionLocalSearch(IloEnv env, Instance inst, const double& timePerIter, const int& gap, Timer<chrono::milliseconds>& timer_cplex,float& opt_time,
//~ const double& timeLimit, float& elapsed_time, double& incumbent){
	//~ Timer<chrono::milliseconds> timer_LS;
	//~ timer_LS.start();	
	//~ unsigned int v,t,a, idPort;
	//~ int i,r,j;
	//~ double objValue = incumbent;
	//~ double prevObj = 1.0e+10;
	//~ 
	//~ cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
	//~ if (gap > 1e-04)
		//~ cplex.setParam(IloCplex::EpGap, gap/100);
	//~ 
	//~ //Get solution values
	//~ //getSolVals(env, inst);
	//~ 
	//~ //FixSolution
	//~ //fixAllSolution(env, inst);
	//~ 
	//~ while((prevObj - objValue > 0.0001) && elapsed_time/1000 < timeLimit){
		//~ prevObj = objValue;
		//~ for(i=1;i>=0;i--){ //Type region (first allow discharging region)
			//~ cout << "TYPE REGION = " << i << endl;
			//~ for(r=0;r<inst.identifyPort[i].getSize();r++){ //For each region				
				//~ for(j=0;j<inst.identifyPort[i][r].getSize();j++){ //Port					
					//~ idPort = inst.identifyPort[i][r][j];						
					//~ #ifndef NFixZvar
					//~ //Allow Z variables to be solved					
					//~ for(v=0;v<inst.speed.getSize();v++){
						//~ for(t=0;t<inst.t;t++){							
							//~ z[v][idPort][t].setBounds(0,1);
						//~ }
					//~ }				
					//~ #endif						
					//~ //Allow X variables to be solved
					//~ for(v=0;v<inst.speed.getSize();v++){
						//~ for(a=0;a<x[v].getSize();a++){
							//~ int j1,j2,t1,t2,type;
							//~ type = inst.getArcType(inst.arcs[v][a], t1, t2, j1, j2);
							//~ //if (type != 0 && ((i==0 && type != 2 && type != 4) || (i==1 && type != 1 && type != 3)) ) //Do not allow incoming regional type arcs
							//~ if (type != 0 && (j1 == idPort || j2 == idPort))
								//~ x[v][a].setBounds(0,1);							
						//~ }
					//~ }
				//~ }					
			//~ }
			//~ //Solve
			//~ cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
			//~ timer_cplex.start();
			//~ cplex.solve();		
			//~ opt_time += timer_cplex.total();
			//~ elapsed_time += timer_LS.total();
			//~ objValue = cplex.getObjValue();
			//~ cout << "Objective: " << objValue << endl;
			//~ if(elapsed_time/1000 >= timeLimit){				
				//~ break;
			//~ }			
			//~ timer_LS.start();
						//~ 
			//~ //Get solution values
			//~ getSolVals(env, inst);			
			//~ 
			//~ //Re-fix ports
			//~ for(r=0;r<inst.identifyPort[i].getSize();r++){ //For each region
				//~ for(j=0;j<inst.identifyPort[i][r].getSize();j++){ //Port					
					//~ idPort = inst.identifyPort[i][r][j];					
					//~ #ifndef NFixZvar
					//~ for(v=0;v<inst.speed.getSize();v++){
						//~ for(t=0;t<inst.t;t++){
							//~ z[v][idPort][t].setBounds(round(zValue[v][idPort][t]),round(zValue[v][idPort][t]));
						//~ }
					//~ }
					//~ #endif
					//~ //Fix X variables										
					//~ for(v=0;v<inst.speed.getSize();v++){						
						//~ for(a=0;a<inst.arcs[v].size();a++){		
							//~ int j1,j2,t1,t2,type;
							//~ type = inst.getArcType(inst.arcs[v][a], t1, t2, j1, j2);
							//~ if (type != 0 && (j1 == idPort || j2 == idPort)){								
								//~ x[v][a].setBounds(round(xValue[v][a]), round(xValue[v][a]));								
							//~ }
						//~ }
					//~ }					
				//~ }				
			//~ }			
		//~ }		
	//~ }
	//~ incumbent = objValue;
//~ }

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
		
		//Relax-and-fix
		double p = T/nIntervals*(1-overLap/100); // Units of t that are add at each iteration to the model.
		int s = T-(T/nIntervals*endBlock); 		 // Last t (relaxed) of model when starting relax-and-fix.
		int sizeInterval = T/nIntervals;		 // Last t of integer block (always starting with 1 integer interval)
        int maxIt = ceil((T-sizeInterval)/p);        
        for (v=1; v <= maxIt; v++){
			timer_cplex.start();
			cout << "Iteration: " << v << "/" << maxIt << " - Solving..." << endl;
			if(!model.cplex.solve())
				cout << model.cplex.getStatus() << endl;
			opt_time += timer_cplex.total();
			cout << "Solution Status " << model.cplex.getStatus() << " Value: " << model.cplex.getObjValue() << endl;
			
			t3S = ceil(T/nIntervals * (v-1) * (1-overLap/100))+1;
			t3F = min(ceil(T/nIntervals * v * (1-overLap/100)),(double)T);

			if (t3S < T){
                cout << "Fixing interval [" << t3S << "," << t3F << "]\n";
                model.fixSolution(env, inst, t3S, t3F, 0);
                cout << "Fixed! \n";
            }
			
			t2S = ceil(s+p*(v-1))+1;
			t2F = min(ceil(s+p*v), (double)T); 
			
			if(t2S <= T){
				cout << "Adding to the model [" << t2S << "," << t2F << "]\n";
				model.decreaseEndBlock (env, inst, nIntervals, t2S, t2F);
			}
			
			t1S = ceil(sizeInterval+p*(v-1)+1);
			t1F = min(ceil(sizeInterval+p*v),(double)T); 
			
			if (t1S <= T){
				cout << "Integralizing: [" << t1S << "," << t1F << "]\n";
				model.reIntegralize(env, inst, t1S, t1F);                
            }
            model.cplex.exportModel("mip_R&F.lp");
            
            double newGap = max(0.001, (gapFirst - (gapFirst/ceil(T-sizeInterval/p))*v) / 100);
            cout << "New GAP " << newGap*100 << " %" << endl;
            model.cplex.setParam(IloCplex::EpGap, newGap);            
		}
		
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
		//~ model.fixAllSolution(env, inst);
		
		//~ double tLimit=0;

		//~ model.fixAndOptTW(env, inst, mIntervals, timePerIterSecond, gapSecond, overlap2, timer_cplex, opt_time, timeLimit/2, elapsed_time, incumbent);
		
		//~ model.improvementPhase(env, inst, mIntervals, timePerIterSecond, gapSecond, overlap2, timer_cplex, opt_time, timeLimit/3, elapsed_time, incumbent);

		//~ tLimit = (timeLimit - elapsed_time/1000)/2;
		//~ model.fixAndOptmizeH(env, inst, timePerIterSecond, gapSecond, incumbent, timer_cplex, opt_time, tLimit, elapsed_time);
		
		//~ model.regionLocalSearch(env, inst,timePerIterSecond, gapSecond, timer_cplex, opt_time, timeLimit, elapsed_time, incumbent);
		
		global_time += timer_global.total();		
		time2ndPhase = elapsed_time;//global_time - time1stPhase;
		obj2ndPhase	= incumbent;//model.cplex.getObjValue();
		
		//For getting information about solution
		//~ model.cplex.setParam(IloCplex::TiLim, 1000);
		//~ model.cplex.setParam(IloCplex::NodeLim, 1);
		//~ model.cplex.solve();
		//~ #ifndef NBetas
		//~ double totalBeta=0;
		//~ IloArray<IloNumArray> betaVals(env, J); 
		//~ for(j=0;j<J;j++){
			//~ betaVals[j] = IloNumArray(env, T);
			//~ model.cplex.getValues(model.beta[j], betaVals[j]);
			//~ totalBeta += IloSum(betaVals[j]);
		//~ }		
		//~ #endif
		
		model.printSolution(env, inst);
		//~ cout << endl
		//~ << nIntervals << "\t"
		//~ << endBlock << "\t"
		//~ << overLap << "\t" 
		//~ << timePerIterFirst << "\t" 
		//~ << gapFirst << "\t"
		//~ << mIntervals << "\t"
		//~ << timePerIterSecond << "\t" 
		//~ << timeLimit << "\t"
		//~ << gapSecond << "\t"
		//~ << opt_time/1000 << "\t"
		//~ << global_time/1000 << "\t"
		//~ << time1stPhase/1000 << "\t"
		//~ << time2ndPhase/1000 << "\t"
		//~ << obj1stPhase << "\t"
		//~ << obj2ndPhase << "\t"
		//~ << abs((obj2ndPhase/obj1stPhase - 1)*100) << "\t"
		//~ //<< model.cplex.getBestObjValue() << "\t"		
		//~ << endl;
		//~ #ifndef NBetas
		//~ cout << "Total betas = " << totalBeta;
		//~ #endif
		
		cout << "First phase -> n : " << nIntervals << " GAP: " << gapFirst << "; f : " << f << "; Overlap: " << overLap << "%; |EndBlock| " << endBlock << "; Time per block " << timePerIterFirst << endl
		//~ << "Second phase -> m: " << mIntervals << "; GAP " << gapSecond << "; Time per local search " << timePerIterSecond << endl
		<< "CPLEX time: " << opt_time/1000 << endl << "Other times: " << (global_time-opt_time)/1000 << endl
		<< "Total time: " << global_time/1000 << endl
		<< "1st phase time: " << time1stPhase/1000 << endl;
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
