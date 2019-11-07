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
#define NBranching
//~ #define NAlpha
//~ #define NBetas  //Negative of alpha
//~ #define NThetas //Plus of alpha
#define WaitAfterOperate		//If defined, allows a vessel to wait after operated at a port.
#define NRelaxation				// Relax all intervals (end block should be 0) for obtaing the lower bound
#define NWWCCReformulation
#define NStartUPValidInequalities
#define NSimplifyModel				//Remove arcs between port i and j for vessel v if min_f_i + min_f_j > Q_v
#define NFixSinkArc

//~ #define FixedGAP
#define NImprovementPhase
#define NRandomTimeInterval				//Defined: Sequential selection of time intervals in the improvementPhase_timeIntervals, otherwise random selection

#define PENALIZATION 1000

ILOSTLBEGIN

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloNumArray> NumNumMatrix;

using namespace std;
using mirp::Model;
using mirp::Instance;

/* Build a model for the first iteration of fix and relax
 * All variables are created (including that are in end block)
 * Just 1st interval have integer variables. 
 * Transition variables (x, w, wB, oB) are relaxed if the next time period belongs to a relaxed interval
 * Variables that are in the end block are not included in any constraint and objective function
 * Constraints of arcs/nodes in the and block are created, but are empty 
 * Params bool (1 = turn on, 0=turn of)
 * - validIneq: Use knapsack valid inequalities
 * - addConstr: Forces a vessel with Q_v < F^Max_it to operante once and depart from i
 * - tightenInvConstr: Tights the inventory constraints of the ports in the last time period when part of the model is ommited
 * - proportionalAlpha: Forbid usign all available alpha before ommiting part of the model
 * - reduceArcs: Does not consider the arcs between port of same type and different regions, or ports i,j of the same region with F^Min_i+F^Min_j > Q_v
 * - tightenFlow: For flow variables fWB(fOB) - for DP change the upper bound, for LP create a lower bound
 */
void Model::buildFixAndRelaxModel(IloEnv& env, Instance inst, const double& nIntervals, const int& endBlock, 
		const bool& validIneq, const bool& addConstr, const bool& tightenInvConstr, const bool& proportionalAlpha,
		const bool& reduceArcs, const bool& tightenFlow){	
	int i, j,t,v,a;	
	int timePerIntervalId = nIntervals; //Number of time periods in each interval
	int J = inst.numTotalPorts;
	int T = inst.t+1;
	double intPart;
	int tOEB = inst.t - (timePerIntervalId*max(0.0,endBlock-modf(nIntervals, &intPart))) ; //Time periods out of End Block (index tOEB+1 is the first in the end block)
	int V = inst.speed.getSize(); //# of vessels
    int N = J + 2;
	cout << "Building model...\n Integer Block: [0," << timePerIntervalId << "]\n Relaxed block: [" << timePerIntervalId+1 << "," << tOEB << "]\n End block: [" << tOEB+1 << "," << T-1 << "]\n";
	
    ///Variables, converters and arrays to storage the values
    #ifndef NAlpha
    alpha = NumVarMatrix(env, N);    
    #endif
	#ifndef NBetas
	beta = NumVarMatrix(env, N);
	#endif
    #ifndef NThetas
    theta = NumVarMatrix(env, N);
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
	sPValue = IloArray<IloNumArray>(env,N);
	IloExpr expr_opd(env);
    if(addConstr){
		operateAndDepart = IloArray<IloArray<IloRangeArray> >(env,V);
	}
    
    #ifndef NWWCCReformulation
	wwcc_relaxation = IloArray<IloRangeArray>(env, N-1);
    #endif
    
    #ifndef NBranching
    //~ sumX = IntVarMatrix(env, V);
    //~ priorityX = IloArray<IloRangeArray>(env, V);
    sumOA = IntVarMatrix(env, V) ;
    priorityOA = IloArray<IloRangeArray>(env, V);
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
        
        if(addConstr){
			operateAndDepart[v] = IloArray<IloRangeArray> (env,N);
        }
        
        #ifndef NBranching
        //~ sumX[v] = IloIntVarArray(env, J+1,0,IloInfinity);
        //~ priorityX[v] = IloRangeArray(env, J+1);
        sumOA[v] = IloIntVarArray(env, J+1, 0,IloInfinity);
        priorityOA[v] = IloRangeArray(env, J+1);
        #endif

		for(j=0;j<N;j++){
			x[v][j] = IloArray<IloBoolVarArray>(env,N);
            convertX[v][j] = IloArray<IloArray<IloConversion> >(env,N);
            xValue[v][j] = IloArray<IloNumArray>(env,N);
			hasArc[v][j] = IloArray<IloIntArray>(env,N);
			arcCost[v][j] = IloArray<IloNumArray>(env,N);
			fX[v][j] = IloArray<IloNumVarArray>(env, N);
            
            for(i=0;i<N;i++){ 
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
					#ifndef NRelaxation
					convertX[v][j][i][t] = IloConversion(env, x[v][j][i][t], ILOFLOAT);
					model.add(convertX[v][j][i][t]);
					#endif
					ss.str(string());
                    //~ if(j > 0 && j <= J){ //If j is a port
                        //~ #ifdef NRelaxation
                        //~ if (i > 0 && i <= J){ //If i is a port
                            //~ if(t + inst.travelTime[v][j-1][i-1] > timePerIntervalId){ 	//If the arrive time is out of first interval, relax it                                
                                //~ convertX[v][j][i][t] = IloConversion(env, x[v][j][i][t], ILOFLOAT);
                                //~ model.add(convertX[v][j][i][t]);
                            //~ }
                        //~ }else if (i == N-1){ // If i is the sink node
                            //~ if(t > timePerIntervalId){ 
                                //~ convertX[v][j][i][t] = IloConversion(env, x[v][j][i][t], ILOFLOAT);
                                //~ model.add(convertX[v][j][i][t]);
                            //~ }
                        //~ }
                        //~ #endif                      
                    //~ }
					ss << "fX_"<<v<<","<<j<<","<<i<<","<<t;
					fX[v][j][i][t] = IloNumVar(env, 0, inst.q_v[v], ss.str().c_str());                    
				}
			}

			if(j>0 && j<=J){ //Only considering Ports (v,j)
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
                
                if(addConstr){
					operateAndDepart[v][j] = IloRangeArray (env,T);
                }
				
				for(t=1;t<T;t++){
					stringstream ss;
					ss << "f_(" << j << "," << t << ")," << v;
					int maxAmountOperation = min(inst.q_v[v], inst.f_max_jt[j-1][t-1]);
					f[v][j][t] = IloNumVar(env, 0, maxAmountOperation, ss.str().c_str());
					ss.str(string());
					ss << "fOA_(" << j << "," << t << ")," << v;
					fOA[v][j][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
					//~ fOA[v][j][t] = IloNumVar(env, 0, inst.q_v[v], ss.str().c_str());
					ss.str(string());
					ss << "z_(" << j << "," << t << ")," << v;
					z[v][j][t].setName(ss.str().c_str());
                    #ifdef NRelaxation
                    if(t > timePerIntervalId){
                        convertZ[v][j][t] = IloConversion(env, z[v][j][t], ILOFLOAT);
                        model.add(convertZ[v][j][t]);
                    }
                    #endif
                    #ifndef NRelaxation
					convertZ[v][j][t] = IloConversion(env, z[v][j][t], ILOFLOAT);
					model.add(convertZ[v][j][t]);
                    #endif
                
                    
					#ifndef WaitAfterOperate
					ss.str(string());
					ss << "oA_(" << j << "," << t << ")," << v;
					oA[v][j][t].setName(ss.str().c_str());
                    #ifdef NRelaxation
                    if(t > timePerIntervalId){
                        convertOA[v][j][t] = IloConversion(env, oA[v][j][t], ILOFLOAT);
                        model.add(convertOA[v][j][t]);
                    }
                    #endif
                    #ifndef NRelaxation
                    convertOA[v][j][t] = IloConversion(env, oA[v][j][t], ILOFLOAT);
					model.add(convertOA[v][j][t]);
                    #endif
					#endif
					if(t<T-1){ //Variables that haven't last time index
						ss.str(string());
						ss << "fW_(" << j << "," << t << ")," << v;
						fW[v][j][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
						//~ fW[v][j][t] = IloNumVar(env, 0, inst.q_v[v], ss.str().c_str());
                        
                        ss.str(string());
						ss << "w_(" << j << "," << t << ")," << v;
						w[v][j][t].setName(ss.str().c_str());
                        #ifdef NRelaxation
                        if(t >= timePerIntervalId){       //Note the use of >= when considering 'horizontal' transition arcs
                            convertW[v][j][t] = IloConversion(env, w[v][j][t], ILOFLOAT);
                            model.add(convertW[v][j][t]);
                        }
						#endif
						#ifndef NRelaxation
						convertW[v][j][t] = IloConversion(env, w[v][j][t], ILOFLOAT);
						model.add(convertW[v][j][t]);
						#endif
						
						#ifndef WaitAfterOperate
                        ss.str(string());
						ss << "fOB_(" << j << "," << t << ")," << v;
						fOB[v][j][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
						//~ fOB[v][j][t] = IloNumVar(env, 0, inst.q_v[v], ss.str().c_str());
                        
						ss.str(string());
						ss << "oB_(" << j << "," << t << ")," << v;
						oB[v][j][t].setName(ss.str().c_str());
                        #ifdef NRelaxation
                        if(t >= timePerIntervalId){
                            convertOB[v][j][t] = IloConversion(env, oB[v][j][t], ILOFLOAT);
                            model.add(convertOB[v][j][t]);
                        }
                        #endif
                        #ifndef NRelaxation
                        convertOB[v][j][t] = IloConversion(env, oB[v][j][t], ILOFLOAT);
						model.add(convertOB[v][j][t]);
						#endif
						#endif
						
						#ifdef WaitAfterOperate
                        ss.str(string());
                        ss << "wB_(" << j << "," << t << ")," << v;
                        wB[v][j][t].setName(ss.str().c_str());
                        #ifdef NRelaxation
                        if(t >= timePerIntervalId){
                            convertWB[v][j][t] = IloConversion(env, wB[v][j][t], ILOFLOAT);
                            model.add(convertWB[v][j][t]);
                        }
                        #endif                        
                        #ifndef NRelaxation
                        convertWB[v][j][t] = IloConversion(env, wB[v][j][t], ILOFLOAT);
						model.add(convertWB[v][j][t]);
						#endif
                        ss.str(string());
                        ss << "fWB_(" << j << "," << t << ")," << v;
                        fWB[v][j][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
                        //~ fWB[v][j][t] = IloNumVar(env, 0, inst.q_v[v], ss.str().c_str());
                        #endif
					}
				}
			}
		}
	}
    //Port time (i,t)
	for(i=1;i<N-1;i++){
		#ifndef NAlpha
		alpha[i] = IloNumVarArray(env,T);
		#endif
		sP[i] = IloNumVarArray(env, T);	
		sPValue[i] = IloNumArray(env,T);
        #ifndef NBetas
        beta[i] = IloNumVarArray(env,T);
        #endif
        #ifndef NThetas
        theta[i] = IloNumVarArray(env,T);
        #endif
        #ifndef NWWCCReformulation
        //Common for loading and discharging ports        
        int num_combinations = (t-1)*t/2; //Number of combinations for k <= j, k,j \in T
        wwcc_relaxation[i] = IloRangeArray(env, num_combinations);        
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
				//Unlimited port capacities in on leg
                //~ if(inst.typePort[i-1]==0)
                    //~ sP[i][t] = IloNumVar(env, -IloInfinity, inst.sMax_jt[i-1][0], ss.str().c_str());
                //~ else
                    //~ sP[i][t] = IloNumVar(env, inst.sMin_jt[i-1][0], IloInfinity, ss.str().c_str());
                //Unlimited port capacity at all
                //~ sP[i][t] = IloNumVar(env, -IloInfinity, IloInfinity, ss.str().c_str());
                
				#ifndef NAlpha
				ss.str(string());
				ss << "alpha_(" << i << "," << t << ")";
				alpha[i][t] = IloNumVar(env, 0, inst.alp_max_jt[i-1][t-1], ss.str().c_str());
				#endif
                #ifndef NBetas
                ss.str(string());
				ss << "beta_(" << i << "," << t << ")";
				//~ beta[i][t] = IloNumVar(env, 0, inst.d_jt[i-1][t-1], ss.str().c_str());
				beta[i][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
                #endif 
                #ifndef NThetas
                ss.str(string());
				ss << "theta_(" << i << "," << t << ")";
				//~ theta[i][t] = IloNumVar(env, 0, inst.d_jt[i-1][t-1], ss.str().c_str());
				theta[i][t] = IloNumVar(env, 0, IloInfinity, ss.str().c_str());
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
				bool addArcIJ = false;
				if(i != j && 
					((inst.delta[i-1] == inst.delta[j-1] && //if i and j are of the same type 
                    inst.f_min_jt[i-1][t-1] + inst.f_min_jt[j-1][t-1] <= inst.q_v[v]) || // and min operation sum is lesser than capacity op vessel
					inst.delta[i-1] != inst.delta[j-1])){ //or ports are of different types
						addArcIJ = true;
					}
				if(reduceArcs && inst.delta[i-1] == inst.delta[j-1] && inst.idRegion[i-1] != inst.idRegion[j-1]){
					addArcIJ = false;
				}
				int t2;				
				double arc_cost;
				if(addArcIJ){ ///Builds arcs between (i,t) to (j,t2) - implicit hasEnteringArc (j,t2). Builds a mirror (j,t)->(i,t2) when possible
					t2 = t + inst.travelTime[v][i-1][j-1]; 
					if (t2<T){ 	//If exists time to reach port j 
						if(hasArc[v][i][j][t] != 1){ //If not added yet
							if (inst.typePort[i-1]==1 && inst.typePort[j-1]==0){ //If port i is consuming and j is a producer port
								arc_cost = inst.costKm[v]*inst.distanceMatrix[i-1][j-1]*(1-inst.trav_empt[v]) + inst.portFee[j-1];
							}else{
								arc_cost = inst.costKm[v]*inst.distanceMatrix[i-1][j-1] + inst.portFee[j-1];
							}
							hasArc[v][i][j][t] = 1;	
							arcCost[v][i][j][t] = arc_cost;                        
							if(t2 <= tOEB)      //If arriving node is in the model                            
								expr1 += arcCost[v][i][j][t]*x[v][i][j][t];
							
							//Relaxing 
							#ifdef NRelaxation
							if(t2 > timePerIntervalId){
								convertX[v][i][j][t] = IloConversion(env, x[v][i][j][t], ILOFLOAT);
								model.add(convertX[v][i][j][t]);
							}
							#endif

							//Check if not set yet
							if(hasEnteringArc1st[v][j][t2] != 1){
								hasEnteringArc1st[v][j][t2] = 1;
								//~ #ifdef NRelaxation
								//~ convertZ[v][j][t2] = IloConversion(env, z[v][j][t2], ILOFLOAT);
								//~ model.add(convertZ[v][j][t2]);
								
								//~ convertW[v][j][t2] = IloConversion(env,w[v][j][t2], ILOFLOAT);
								//~ model.add(convertW[v][j][t2]);
								
								//~ #ifdef WaitAfterOperate
								//~ convertWB[v][j][t2] = IloConversion(env,wB[v][j][t2], ILOFLOAT);
								//~ model.add(convertWB[v][j][t2]);
								//~ #endif
								
								//~ #ifndef WaitAfterOperate
								//~ convertOB[v][j][t2] = IloConversion(env,oB[v][j][t2], ILOFLOAT);
								//~ model.add(convertOB[v][j][t2]);
								//~ convertOA[v][j][t2] = IloConversion(env,oA[v][j][t2], ILOFLOAT);
								//~ model.add(convertOA[v][j][t2]);
								//~ #endif  
								//~ #endif
							}						
						}						
						//when time t reach a time that can be built a mirror between j and i
						if(hasArc[v][j][i][t] != 1){
							if (t >= inst.firstTimeAv[v]+1 + inst.travelTime[v][i-1][j-1]){ 
								if (inst.typePort[j-1]==1 && inst.typePort[i-1]==0){ 		  //If port j is consuming and i is a producer port
									arc_cost = inst.costKm[v]*inst.distanceMatrix[j-1][i-1]*(1-inst.trav_empt[v]) + inst.portFee[i-1];		
								}else{
									arc_cost = inst.costKm[v]*inst.distanceMatrix[j-1][i-1] + inst.portFee[i-1];
								}
								hasArc[v][j][i][t] = 1;
								arcCost[v][j][i][t] = arc_cost;
								
								if(t2 <= tOEB)      //If arriving node is in the model
									expr1 += arcCost[v][j][i][t]*x[v][j][i][t];
								//Relaxing 
								#ifdef NRelaxation
								if(t2 > timePerIntervalId){ 
									convertX[v][j][i][t] = IloConversion(env, x[v][j][i][t], ILOFLOAT);
									model.add(convertX[v][j][i][t]);
									if(hasEnteringArc1st[v][i][t2] != 1){
										hasEnteringArc1st[v][i][t2] = 1; 
										//~ convertZ[v][i][t2] = IloConversion(env, z[v][i][t2], ILOFLOAT);
										//~ model.add(convertZ[v][i][t2]);
										
										//~ convertW[v][i][t2] = IloConversion(env,w[v][i][t2], ILOFLOAT);
										//~ model.add(convertW[v][i][t2]);
										
										//~ #ifdef WaitAfterOperate
										//~ convertWB[v][i][t2] = IloConversion(env,wB[v][i][t2], ILOFLOAT);
										//~ model.add(convertWB[v][i][t2]);
										//~ #endif
										
										//~ #ifndef WaitAfterOperate
										//~ convertOB[v][i][t2] = IloConversion(env,oB[v][i][t2], ILOFLOAT);
										//~ model.add(convertOB[v][i][t2]);
										//~ convertOA[v][i][t2] = IloConversion(env,oA[v][i][t2], ILOFLOAT);
										//~ model.add(convertOA[v][i][t2]);
										//~ #endif         
									}
								}
								#endif
							}
						}
					}
				}		
				///For cases where i and j are not compatible, but j is can be compatible with other j2
				///Builds waitJ, sinkJ, (j,t2)->(j2,t3), enteringJ2, 
				if(addArcIJ || hasEnteringArc1st[v][j][t] == 1){
					if(!addArcIJ){ //When i and j are not compatible
						t2 = t;
					}
					if(t2<T){					
						if (t2+1<T-1){ 			//Waiting arc of j
							if(hasEnteringArc1st[v][j][t2+1] != 1){
								hasEnteringArc1st[v][j][t2+1] = 1;
								//~ #ifdef NRelaxation
								//~ convertZ[v][j][t2+1] = IloConversion(env, z[v][j][t2+1], ILOFLOAT);
								//~ model.add(convertZ[v][j][t2+1]);
								
								//~ convertW[v][j][t2+1] = IloConversion(env,w[v][j][t2+1], ILOFLOAT);
								//~ model.add(convertW[v][j][t2+1]);
								
								//~ #ifdef WaitAfterOperate
								//~ convertWB[v][j][t2+1] = IloConversion(env,wB[v][j][t2+1], ILOFLOAT);
								//~ model.add(convertWB[v][j][t2+1]);
								//~ #endif
								
								//~ #ifndef WaitAfterOperate
								//~ convertOB[v][j][t2+1] = IloConversion(env,oB[v][j][t2+1], ILOFLOAT);
								//~ model.add(convertOB[v][j][t2+1]);
								//~ convertOA[v][j][t2+1] = IloConversion(env,oA[v][j][t2+1], ILOFLOAT);
								//~ model.add(convertOA[v][j][t2+1]);
								//~ #endif  
								//~ #endif
							}
						}
						//Sink arc from port j 
						if (hasArc[v][j][N-1][t2] != 1){
							hasArc[v][j][N-1][t2] = 1;
							arcCost[v][j][N-1][t2] = -(T-t2-1)*inst.perPeriodRewardForFinishingEarly;
							if(t2 <= tOEB)      //If arriving node is in the model
								expr1 += arcCost[v][j][N-1][t2]*x[v][j][N-1][t2];
							//Relaxing 
							#ifdef NRelaxation
							if(t2 > timePerIntervalId){ 
								convertX[v][j][N-1][t2] = IloConversion(env, x[v][j][N-1][t2], ILOFLOAT);
								model.add(convertX[v][j][N-1][t2]);
							}
							#endif
						}						
						
						//Create arc from j,t2 to others ports (j2) in time t3
						for(int j2=1;j2<=J;j2++){
							bool addArcJJ2 = false;
							if(j2 != i && j2 != j && 
								((inst.delta[j-1] == inst.delta[j2-1] && //if j and j2 are of the same type 
								inst.f_min_jt[j-1][0] + inst.f_min_jt[j2-1][0] <= inst.q_v[v]) || // and min operation sum is lesser than capacity op vessel
								inst.delta[j-1] != inst.delta[j2-1])){ //or ports are of different types
									addArcJJ2 = true;
								}
							if(reduceArcs && inst.delta[j-1] == inst.delta[j2-1] && inst.idRegion[j-1] != inst.idRegion[j2-1]){ //If they are of the same type but from different regions, no add the arc
								addArcJJ2 = false;
							}			
							
							if(addArcJJ2){
								int t3 = t2+inst.travelTime[v][j-1][j2-1];  
								if(hasArc[v][j][j2][t2] != 1){
									if(t3<T){
										if (inst.typePort[j-1]==1 && inst.typePort[j2-1]==0){
											arc_cost = inst.costKm[v]*inst.distanceMatrix[j-1][j2-1]*(1-inst.trav_empt[v]) + inst.portFee[j2-1];		
										}else{
											arc_cost = inst.costKm[v]*inst.distanceMatrix[j-1][j2-1] + inst.portFee[j2-1];
										}
										hasArc[v][j][j2][t2] = 1;
										arcCost[v][j][j2][t2] = arc_cost;
										if(t3 <= tOEB)      //Ir arriving node is in the model
											expr1 += arcCost[v][j][j2][t2]*x[v][j][j2][t2];
										//Relaxing 
										#ifdef NRelaxation
										if(t3 > timePerIntervalId){ 
											convertX[v][j][j2][t2] = IloConversion(env, x[v][j][j2][t2], ILOFLOAT);
											model.add(convertX[v][j][j2][t2]);
										}
										#endif
										if(hasEnteringArc1st[v][j2][t3] != 1){
											hasEnteringArc1st[v][j2][t3] = 1;									
											//~ #ifdef NRelaxation
											//~ convertZ[v][j2][t3] = IloConversion(env, z[v][j2][t3], ILOFLOAT);
											//~ model.add(convertZ[v][j2][t3]);
											
											//~ convertW[v][j2][t3] = IloConversion(env,w[v][j2][t3], ILOFLOAT);
											//~ model.add(convertW[v][j2][t3]);
											
											//~ #ifdef WaitAfterOperate
											//~ convertWB[v][j2][t3] = IloConversion(env,wB[v][j2][t3], ILOFLOAT);
											//~ model.add(convertWB[v][j2][t3]);
											//~ #endif
											
											//~ #ifndef WaitAfterOperate
											//~ convertOB[v][j2][t3] = IloConversion(env,oB[v][j2][t3], ILOFLOAT);
											//~ model.add(convertOB[v][j2][t3]);
											//~ convertOA[v][j2][t3] = IloConversion(env,oA[v][j2][t3], ILOFLOAT);
											//~ model.add(convertOA[v][j2][t3]);
											//~ #endif  
											//~ #endif
										}
									}
								}
							}
						}
					}
				}
			}
			if(t+1<T){ //Waitinng i
				if(hasEnteringArc1st[v][i][t+1] != 1){
					hasEnteringArc1st[v][i][t+1] = 1;
					
					//~ #ifdef NRelaxation
					//~ convertZ[v][i][t+1] = IloConversion(env, z[v][i][t+1], ILOFLOAT);
					//~ model.add(convertZ[v][i][t+1]);
					
					//~ convertW[v][i][t+1] = IloConversion(env,w[v][i][t+1], ILOFLOAT);
					//~ model.add(convertW[v][i][t+1]);
					
					//~ #ifdef WaitAfterOperate
					//~ convertWB[v][i][t+1] = IloConversion(env,wB[v][i][t+1], ILOFLOAT);
					//~ model.add(convertWB[v][i][t+1]);
					//~ #endif
					
					//~ #ifndef WaitAfterOperate
					//~ convertOB[v][i][t+1] = IloConversion(env,oB[v][i][t+1], ILOFLOAT);
					//~ model.add(convertOB[v][i][t+1]);
					//~ convertOA[v][i][t+1] = IloConversion(env,oA[v][i][t+1], ILOFLOAT);
					//~ model.add(convertOA[v][i][t+1]);
					//~ #endif  
					//~ #endif
				}
			}
			//Sink arc from port i
			if(hasArc[v][i][N-1][t] != 1){
				hasArc[v][i][N-1][t] = 1;
				arcCost[v][i][N-1][t] = -(T-t-1)*inst.perPeriodRewardForFinishingEarly;
				if(t <= tOEB)      //If arriving node is in the model
					expr1 += arcCost[v][i][N-1][t]*x[v][i][N-1][t];
				//Relaxing 
				#ifdef NRelaxation
				if(t > timePerIntervalId){ 
					convertX[v][i][N-1][t] = IloConversion(env, x[v][i][N-1][t], ILOFLOAT);
					model.add(convertX[v][i][N-1][t]);
				}
				#endif
			}
		}

		for(i=1;i<N-1;i++){		//Only considering ports
			for(t=1;t<=tOEB;t++){
				if(hasEnteringArc1st[v][i][t]==1){
					expr += inst.r_jt[i-1][t-1]*f[v][i][t];									//1st term
					expr1 += (t-1)*inst.attemptCost*z[v][i][t];								//3rd term
				}
			}
		}
	}
	
	for(j=1;j<N-1;j++){
		for(t=1;t<=tOEB;t++){			
			#ifndef NAlpha
			expr1 += inst.p_jt[j-1][t-1]*alpha[j][t];									//4rd term
			#endif
            #ifndef NBetas
            expr1 += PENALIZATION*beta[j][t]; //100 may be enough
            #endif
            #ifndef NThetas
            expr1 += PENALIZATION*theta[j][t];
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
	#ifndef NAlpha
	cumSlack = IloRangeArray(env,N-1);
	IloExpr expr_cumSlack(env);
	#endif
	
	//Limits on operation values
	operationLowerLimit = IloArray<IloArray<IloRangeArray> >(env,V);
	operationUpperLimit = IloArray<IloArray<IloRangeArray> >(env,V);
	
	//Flow upper limits
	flowCapacityX = IloArray<IloArray<IloArray<IloRangeArray> > > (env,V);
	flowCapacityOA = IloArray<IloArray<IloRangeArray> >(env,V);
	flowCapacityW = IloArray<IloArray<IloRangeArray> >(env,V);
	#ifndef WaitAfterOperate
	flowCapacityOB = IloArray<IloArray<IloRangeArray> >(env,V);
	#endif
	#ifdef WaitAfterOperate
	flowCapacityWB = IloArray<IloArray<IloRangeArray> >(env,V);
	#endif
    //Flow lower limits
    if(tightenFlow){
		flowMinCapacityX = IloArray<IloArray<IloArray<IloRangeArray> > > (env,V);
		flowMinCapacityOA = IloArray<IloArray<IloRangeArray> >(env,V);
		flowMinCapacityW = IloArray<IloArray<IloRangeArray> >(env,V);
		#ifndef WaitAfterOperate	
		flowMinCapacityOB = IloArray<IloArray<IloRangeArray> >(env,V);
		#endif
		#ifdef WaitAfterOperate	
		flowMinCapacityWB = IloArray<IloArray<IloRangeArray> >(env,V);
		#endif
	}
    //WWCC reformulation
    #ifndef NWWCCReformulation
    IloExpr expr_sumF(env);
    IloExpr expr_sumO(env);
    IloExpr expr_wwcc(env);
    int it_kt;
    #endif
    
    #ifndef NStartUPValidInequalities
    startup_sumStartIfOperate = IloArray<IloRangeArray> (env, J+1);
    startup_sumForceStartup = IloArray<IloRangeArray> (env, J+1);
    startup_dlsccs = IloArray<IloRangeArray> (env, J);
    startup_validInequality = IloArray<IloRangeArray> (env, J+1);
    
    IloExpr expr_sum_vInV_OA;
    IloExpr expr_sum_vInV_O;
    IloExpr expr_sum_vInV_O_t_1;
    #endif
	
	
	if(validIneq){
		knapsack_P_1 = IloArray<IloRangeArray> (env, J+1);			//Altough created for all J ports, it is only considered for loading(production) or consumption(discharging) ports. We create J+1 as index j starts with 1
		knapsack_P_2 = IloArray<IloRangeArray> (env, J+1);	
		knapsack_D_1 = IloArray<IloRangeArray> (env, J+1);	
		knapsack_D_2 = IloArray<IloRangeArray> (env, J+1);	
		#ifndef WaitAfterOperate
		knapsack_D_3 = IloArray<IloRangeArray> (env, J+1);	
		#endif
	}
	unsigned int num_combinations = (T-3)*2;
	for(i=1;i<=J;i++){
        #ifndef NWWCCReformulation
        it_kt = 0;
        #endif
		stringstream ss1;
		berthLimit[i] = IloRangeArray(env,T);
		#ifndef NAlpha
		expr_cumSlack.clear();
		#endif
		portInventory[i] = IloRangeArray(env,T);
		if(validIneq){
			if (inst.typePort[i-1] == 0){ 	//Loading port
				knapsack_P_1[i] = IloRangeArray(env, num_combinations);		//(T-3)*2 = Number of combinations 1,...t + t,...,T for all t \in T. Using T-3 because de increment of T = inst.t+1
				knapsack_P_2[i] = IloRangeArray(env, num_combinations);
			}else{							//Discharging port
				knapsack_D_1[i] = IloRangeArray(env, num_combinations);
				knapsack_D_2[i] = IloRangeArray(env, num_combinations);
				#ifndef WaitAfterOperate
				knapsack_D_3[i] = IloRangeArray(env, num_combinations);
				#endif
			}        
			int it,it2=0,k,l;
			IloExpr expr_kP1_LHS(env), expr_kP2_LHS(env), expr_kD1_LHS(env), expr_kD2_LHS(env);		
			for(it=0;it< num_combinations-(T-1-tOEB);it++){	//For each valid inequality - limited to the constraint that uses the interval [tOEB-1, tOEB]
				double kP1_RHS=0, kP2_RHS=0, kD1_RHS=0, kD2_RHS=0, sum_alphaMax=0, alphaUB=0;
				expr_kP1_LHS.clear();
				expr_kP2_LHS.clear(); //Also used for the equivalent kD3
				expr_kD1_LHS.clear();
				expr_kD2_LHS.clear();
				//Definining the size of set T =[l,k] 
				l=1;
				k=tOEB;
				if(it<tOEB-2){
					k = it+2;
					it2++;		//It2 gets the same value of it until it reach the value tOEB-2
				}else if(it == tOEB-2){
					it = T-3;	//it jumps to half of array of the IloRangeArray (starting the constraints of range [t,...|T|], t \in T) - whent T=45, starts in 43
					l = it2+4-k;
					it2++;
				}
				else{ //After passed the half of array
					l = it2+4-k;
					it2++;
				}
				//~ if(i==1)
					//~ //~ //~ cout it << " [" << l << "," << k << "]\n";
				for(v=0;v<V;v++){
					#ifdef WaitAfterOperate
					if(k<tOEB && hasEnteringArc1st[v][i][k]){
						expr_kP1_LHS += wB[v][i][k];
					}
					if(l>1 && hasEnteringArc1st[v][i][l-1]==1)
						expr_kD1_LHS += w[v][i][l-1] + wB[v][i][l-1];
					#endif
					#ifndef WaitAfterOperate
					expr_kP1_LHS += oB[v][i][k];
					if(l>1 && hasEnteringArc1st[v][i][l-1]==1)
						expr_kD1_LHS += w[v][i][l-1] + oB[v][i][l-1];
					#endif
					#ifndef WaitAfterOperate
					if(l>1 && hasEnteringArc1st[v][i][l-1]==1)
						expr_kD2_LHS += oB[v][i][l-1];
					#endif
					for(t=l;t<=k;t++){
						if(hasEnteringArc1st[v][i][t]==1){
							expr_kP2_LHS += z[v][i][t];	
							#ifndef WaitAfterOperate
							expr_kD2_LHS += oA[v][i][t];
							#endif
							#ifdef WaitAfterOperate
							expr_kD2_LHS += z[v][i][t];
							#endif
						}
						for(j=0;j<N;j++){
							if(j==0 && inst.initialPort[v]+1 == i && t==inst.firstTimeAv[v]+1) 	//Source arc
								expr_kD1_LHS += x[v][j][i][0];
							else if (j == N-1 && hasArc[v][i][j][t] == 1)						//Sink arc
								expr_kP1_LHS += x[v][i][j][t];
							else if (j > 0 && j <= J){											//Port arcs	
								if(i != j){
									//~ if(hasArc[v][i][j][t] == 1)  //Ignoring the rule that the arriving node must be in the model
									if(hasArc[v][i][j][t] == 1 && t+inst.travelTime[v][i-1][j-1] <= tOEB)  //If arc exists and arrives at port j in the integer or relaxed block
										expr_kP1_LHS += x[v][i][j][t];
									if(t - inst.travelTime[v][j-1][i-1] > 0){					//Only if it is possible an entering arc due the travel time
										if(hasArc[v][j][i][t-inst.travelTime[v][j-1][i-1]] == 1)
											expr_kD1_LHS += x[v][j][i][t-inst.travelTime[v][j-1][i-1]];
									}
								}
							}
						}
						//Port time iterator (i,t)
						if (v==0){
							kP1_RHS += inst.d_jt[i-1][t-1];
							kD1_RHS += inst.d_jt[i-1][t-1];
							#ifndef NAlpha
							sum_alphaMax += inst.alp_max_jt[i-1][t-1];
							#endif
						}
					}
				}
				if(l==1){
					kP1_RHS += -inst.sMax_jt[i-1][0] + inst.s_j0[i-1];
					kD1_RHS += -inst.s_j0[i-1] + inst.sMin_jt[i-1][k-1];
					
				}else{
					kP1_RHS += -inst.sMax_jt[i-1][0] + inst.sMin_jt[i-1][0];
					kD1_RHS += -inst.sMax_jt[i-1][l-2] + inst.sMin_jt[i-1][k-1]; //equivalent 'l-1' and 'k' if index starts in 1                
				}
				///If considering alpha parameters, otherwise comment the above 2 lines
				#ifndef NAlpha
				kP1_RHS += - min(sum_alphaMax, inst.alp_max_j[i-1]);
				kD1_RHS += - min(sum_alphaMax, inst.alp_max_j[i-1]);
				#endif
				
				kP2_RHS = max(0.0,ceil(kP1_RHS/inst.f_max_jt[i-1][0]));
				kP1_RHS = max(0.0,ceil(kP1_RHS/inst.maxVesselCapacity));
				kD2_RHS = max(0.0,ceil(kD1_RHS/inst.f_max_jt[i-1][0]));
				kD1_RHS = max(0.0,ceil(kD1_RHS/inst.maxVesselCapacity));
				
				stringstream ss, ss1, ss2;
				if(inst.typePort[i-1] == 0){
					ss << "knpasackP1_" << i << "_(" << l << "," << k << ")";
					knapsack_P_1[i][it] = IloRange(env, kP1_RHS, expr_kP1_LHS, IloInfinity, ss.str().c_str());
					ss1 << "knapsackP2_" << i << "_(" << l << "," << k << ")";
					knapsack_P_2[i][it] = IloRange(env, kP2_RHS, expr_kP2_LHS, IloInfinity, ss1.str().c_str());
					model.add(knapsack_P_1[i][it]);
					model.add(knapsack_P_2[i][it]);
				}else{
					ss << "knapsackD1_" << i << "_(" << l << "," << k << ")";
					knapsack_D_1[i][it] = IloRange(env, kD1_RHS, expr_kD1_LHS, IloInfinity, ss.str().c_str());
					ss1 << "knapsackD2_" << i << "_(" << l << "," << k << ")";
					#ifndef WaitAfterOperate
					knapsack_D_2[i][it] = IloRange(env, kD1_RHS, expr_kD2_LHS, IloInfinity, ss1.str().c_str());
					#endif
					#ifdef WaitAfterOperate
					knapsack_D_2[i][it] = IloRange(env, kD2_RHS, expr_kD2_LHS, IloInfinity, ss1.str().c_str());
					#endif
					model.add(knapsack_D_1[i][it]);
					model.add(knapsack_D_2[i][it]);
					#ifndef WaitAfterOperate
					ss2 << "knapsackD3_" << i << "_(" << l << "," << k << ")";
					knapsack_D_3[i][it] = IloRange(env, kD2_RHS, expr_kP2_LHS, IloInfinity, ss2.str().c_str());
					model.add(knapsack_D_3[i][it]);
					#endif
				}
			}
		}
        
        #ifndef WaitAfterOperate
        #ifndef NStartUPValidInequalities
        startup_sumStartIfOperate[i] = IloRangeArray(env, T);
        startup_sumForceStartup[i] = IloRangeArray(env, T);
        startup_dlsccs[i] = IloRangeArray(env, T);
        startup_validInequality[i] = IloRangeArray(env);
        
        //TODO after tested, implement in optimized format
        for(t=1;t<=tOEB;t++){
            expr_sum_vInV_O.clear();
            expr_sum_vInV_O_t_1.clear();
            expr_sum_vInV_OA.clear();
            for(v=0;v<V;v++){
                expr_sum_vInV_O += z[v][i][t];
                expr_sum_vInV_OA += oA[v][i][t];
                if(t>1)
                    expr_sum_vInV_O_t_1 += z[v][i][t-1];
            }
            stringstream ss1,ss2,ss3;
            ss1 << "startUp_sumStartIfOperate_(" << i << "," << t << ")";
            ss2 << "startUp_sumForceStartUp_(" << i << "," << t << ")";
            ss3 << "startUp_DLSCCS_(" << i << "," << t << ")";
            startup_sumStartIfOperate[i][t] = IloRange(env, -IloInfinity, expr_sum_vInV_OA - expr_sum_vInV_O, 0, ss1.str().c_str());
            model.add(startup_sumStartIfOperate[i][t]);
            startup_sumForceStartup[i][t] = IloRange(env, -IloInfinity, expr_sum_vInV_O - expr_sum_vInV_O_t_1 - expr_sum_vInV_OA, 0, ss2.str().c_str());
            model.add(startup_sumForceStartup[i][t]);
            
            expr_sum_vInV_O.clear();
            for(int u=1;u<=t;u++){
                for(v=0;v<V;v++)
                    expr_sum_vInV_O += z[v][i][t];
            }
            startup_dlsccs[i][t] = IloRange(env, inst.lb_oper_jt[i-1][t-1], expr_sum_vInV_O, IloInfinity, ss3.str().c_str());
            model.add(startup_dlsccs[i][t]);
            
            //~ cplex.addUserCut();
        }
        #endif
        #endif
		float thighetInventoryValue = 0;
		for(t=1;t<=tOEB;t++){
			expr_berth.clear();
			expr_invBalancePort.clear();
			stringstream ss, ss2;
			ss << "berthLimit_(" << i << "," << t << ")";
			bool emptyExpr = true;
            #ifndef NWWCCReformulation
            expr_sumF.clear();
            expr_sumO.clear();
            #endif
			for(v=0;v<V;v++){
				if (hasEnteringArc1st[v][i][t]==1){ //Only if exists an entering arc in the node				
					expr_berth += z[v][i][t];
					emptyExpr = false;
                    #ifndef NWWCCReformulation                    
                    expr_sumF += f[v][i][t];
                    expr_sumO += z[v][i][t];                    
                    #endif
					expr_invBalancePort += -f[v][i][t];
				}
			}
			if(!emptyExpr){ //Only if there are some Z variable in the expr
				berthLimit[i][t] = IloRange(env,-IloInfinity, expr_berth, inst.b_j[i-1], ss.str().c_str());
				model.add(berthLimit[i][t]);
			}
			
			#ifndef NAlpha
			expr_cumSlack += alpha[i][t];
			expr_invBalancePort += -alpha[i][t];
			#endif
			#ifndef NBetas
            expr_invBalancePort += beta[i][t];
            #endif
            #ifndef NThetas
            expr_invBalancePort += -theta[i][t];
            #endif
            
			ss2 << "invBalancePort_(" << i << "," << t << ")";
			portInventory[i][t] = IloRange(env, inst.delta[i-1]*inst.d_jt[i-1][t-1],
				sP[i][t]-sP[i][t-1]-inst.delta[i-1]*expr_invBalancePort,
				inst.delta[i-1]*inst.d_jt[i-1][t-1], ss2.str().c_str());
			model.add(portInventory[i][t]);
					
            #ifndef NWWCCReformulation            
            stringstream ss3;
            //Common for both loading and discharging ports
            for(int k=1;k<=t;k++){
                ss3.str(string());
                ss3 << it_kt << "_wwcc_relaxation_" << i << "(" << k << "," << t << ")";                    
                expr_wwcc.clear();
                int wwcc_rhs=0;
                bool hasExpr=false;
                for(int u=k;u<=t;u++){
                    for (v=0;v<V;v++){
                        if(hasEnteringArc1st[v][i][u]==1){
                            expr_wwcc += z[v][i][u];
                            hasExpr = true;
                        }
                    }
                    if (inst.typePort[i-1] == 1) //Discharging port
                        wwcc_rhs += inst.dM_jt[i-1][u];
                    else    //Loading port
                        wwcc_rhs += inst.d_jt[i-1][u-1];
                }
                expr_wwcc *= inst.f_max_jt[i-1][t-1];
                if(k==1){ //When k=1 and k-1=0, net inventory variable is 0
                    if(inst.typePort[i-1] == 0){ //Loading port
                        wwcc_rhs += inst.s_j0[i-1] - inst.sMax_jt[i-1][1];
                    }
                }else{ 
                    if(inst.typePort[i-1] == 1){ //Discharging port
                        expr_wwcc += sP[i][k-1];
                        wwcc_rhs += inst.sMinM_jt[i-1][k-1];
                    }else{  //Loading port
                        expr_wwcc -= sP[i][k-1];
                        wwcc_rhs -= inst.sMax_jt[i-1][1]; // == -\overline{S}_i
                    }                    
                }                
                if(hasExpr || k>1){                    
                    //~ //~ //~ cout " " <<  expr_wwcc << " >= " << wwcc_rhs << endl;
                    wwcc_relaxation[i][it_kt] = IloRange(env, wwcc_rhs, expr_wwcc, IloInfinity, ss3.str().c_str());                    
                    model.add(wwcc_relaxation[i][it_kt]);
                }
                //~ //~ //~ cout ss3.str() << endl;
                it_kt++;
            }
            #endif
		}
		#ifndef NAlpha
		ss1 << "cum_slack("<<i<<")";		
		if(proportionalAlpha && nIntervals > 1 && endBlock > 0){
			cumSlack[i] = IloRange(env, expr_cumSlack, inst.alp_max_j[i-1]/(T-1)*tOEB, ss1.str().c_str());
		}else{
			cumSlack[i] = IloRange(env, expr_cumSlack, inst.alp_max_j[i-1], ss1.str().c_str());
		}
		model.add(cumSlack[i]);  
		#endif
		if(tightenInvConstr && nIntervals > 1 && endBlock > 0){
			if(inst.delta[i-1] == 1){
				sP[i][tOEB].setUB(inst.sMax_jt[i-1][0]-inst.sMax_jt[i-1][0]*0.1);
			}else{
				sP[i][tOEB].setLB(inst.sMin_jt[i-1][0] + inst.sMax_jt[i-1][0]*0.1);
			}
		}    
	}
    //~ //~ //~ cout "it_kt = " << it_kt << endl;
	#ifndef NAlpha
	expr_cumSlack.end();
	#endif
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
		
        if(tightenFlow){
			flowMinCapacityX[v] = IloArray<IloArray<IloRangeArray> >(env,N);
			flowMinCapacityOA[v] = IloArray<IloRangeArray> (env,N-1);
			#ifndef WaitAfterOperate
			flowMinCapacityOB[v] = IloArray<IloRangeArray> (env,N-1);
			#endif
			flowMinCapacityW[v] = IloArray<IloRangeArray> (env,N-1);
			#ifdef WaitAfterOperate
			flowMinCapacityWB[v] = IloArray<IloRangeArray> (env,N-1);
			#endif
		}
        
		for(i=0;i<N;i++){
			if(i>0 && i <= J){ //Only considering ports				
				firstLevelBalance[v][i] = IloRangeArray(env,T,0,0); //Id 0 is not used
				secondLevelBalance[v][i] = IloRangeArray(env,T,0,0); 
				firstLevelFlow[v][i] = IloRangeArray(env,T,0,0); 
				secondLevelFlow[v][i] = IloRangeArray(env,T,0,0); 
				#ifndef WaitAfterOperate
				linkBalance[v][i] = IloRangeArray(env,T,-IloInfinity,0); 
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
                                if(hasArc[v][i][j][t] == 1 &&
									t + inst.travelTime[v][i-1][j-1] <= tOEB){ //Only if arc departs and arrives in the integer or relaxed block
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
								if(hasArc[v][i][j][t]==1){
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
				flowCapacityOA[v][i] = IloRangeArray(env,T);
				#ifndef WaitAfterOperate
				flowCapacityOB[v][i] = IloRangeArray(env,T);
				#endif
				flowCapacityW[v][i] = IloRangeArray(env,T);
				#ifdef WaitAfterOperate
				flowCapacityWB[v][i] = IloRangeArray(env,T);
				#endif
                
                if(tightenFlow){
					flowMinCapacityOA[v][i] = IloRangeArray(env,T);
					#ifndef WaitAfterOperate
					flowMinCapacityOB[v][i] = IloRangeArray(env,T);
					#endif
					flowMinCapacityW[v][i] = IloRangeArray(env,T);
					#ifdef WaitAfterOperate
					flowMinCapacityWB[v][i] = IloRangeArray(env,T);
					#endif
                }
                #ifndef NBranching
                //~ IloExpr expr_sumEnteringX(env);
                IloExpr expr_sumOA(env);                    
                #endif
				
				for(t=1;t<=tOEB;t++){
					stringstream ss, ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8, ss9, ss10, ss11;
					ss << "First_level_balance_" << v << "(" << i << "," << t << ")";
					ss1 << "Second_level_balance_" << v << "(" << i << "," << t << ")";
					ss2 << "link_balance_" << v << "(" << i << "," << t << ")";
					ss3 << "First_level_flow_" << v << "(" << i << "," << t << ")";
					ss4 << "Second_level_flow_" << v << "(" << i << "," << t << ")";
                    ss11 << "operate_and_depart_" << v << "(" << i << "," << t << ")";
					
					if(hasArc[v][i][N-1][t]==1)
						expr_sinkFlow += x[v][i][N-1][t];
					
					expr_1stLevel.clear();
					expr_1stFlow.clear();
					expr_2ndLevel.clear();
					expr_2ndFlow.clear();

                    if(addConstr){
						expr_opd.clear();
					}
                    
					for(j=0;j<N;j++){
						if(j<N-1){ //No consider sink arc (first level balance)
							//If j is the source node and reach i in time t
							if(j==0 && inst.initialPort[v]+1 == i && t==inst.firstTimeAv[v]+1){
								expr_1stLevel += x[v][j][i][0];
								expr_1stFlow += fX[v][j][i][0];
								fX[v][j][i][0].setBounds(inst.s_v0[v], inst.s_v0[v]); //Fixing initial inventory 
                                #ifndef NBranching
                                //~ expr_sumEnteringX += x[v][j][i][0];
                                #endif
							}
							else if (j>0){ //When j is a port
								if (t - inst.travelTime[v][j-1][i-1] >= 0){ //If is possible to exist an arc from j to i
									if(hasArc[v][j][i][t-inst.travelTime[v][j-1][i-1]] == 1){ //If the arc exists
										expr_1stLevel += x[v][j][i][t-inst.travelTime[v][j-1][i-1]];
										expr_1stFlow += fX[v][j][i][t-inst.travelTime[v][j-1][i-1]];
                                        #ifndef NBranching
                                        //~ expr_sumEnteringX += x[v][j][i][t-inst.travelTime[v][j-1][i-1]];
                                        #endif
									}
								}                                
							}
						}
						if(j>0){ //No consider source arc (second level balance)
                            if (hasArc[v][i][j][t] == 1){
                                if(j == N-1){ //If j is the sink node
                                    expr_2ndLevel += x[v][i][j][t];
                                    expr_2ndFlow += fX[v][i][j][t];
                                    if(addConstr){
										if (inst.q_v[v] <= inst.f_max_jt[i-1][0]){ //Only if a vessel can load(unload) in the port in just 1 time period
											expr_opd += x[v][i][j][t];
										}
                                    }
                                }else if (t + inst.travelTime[v][i-1][j-1] <= tOEB){ //If j is a port, it is then necessary that the arrival is in the model
                                    expr_2ndLevel += x[v][i][j][t];
                                    expr_2ndFlow += fX[v][i][j][t];
                                    if(addConstr){
										if (inst.q_v[v] <= inst.f_max_jt[i-1][0] && inst.typePort[i-1] != inst.typePort[j-1]){ //Only if a vessel can load(unload) in the port in just 1 time period and i and j are ports of different types
											expr_opd += x[v][i][j][t];
										}
                                    }
                                }
							}
						}
					}
                    if(addConstr){
                        if(hasEnteringArc1st[v][i][t]==1 && inst.q_v[v] <= inst.f_max_jt[i-1][t-1]){
							#ifdef WaitAfterOperate                        
								expr_opd += -z[v][i][t];
							#endif
							#ifndef WaitAfterOperate                        
								expr_opd += -oA[v][i][t];
							#endif
							operateAndDepart[v][i][t] = IloRange(env, 0, expr_opd, 0,ss11.str().c_str());                    
							model.add(operateAndDepart[v][i][t]);
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
					if(hasEnteringArc1st[v][i][t]==1){
						model.add(firstLevelBalance[v][i][t]);
											
						secondLevelBalance[v][i][t].setName(ss1.str().c_str());
						model.add(secondLevelBalance[v][i][t]);
					}
					
					
					
					#ifndef WaitAfterOperate                    
					linkBalance[v][i][t].setName(ss2.str().c_str());
					model.add(linkBalance[v][i][t]);
					#endif
					
					if(hasEnteringArc1st[v][i][t]==1){
					
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
						stringstream ss7_1;
						ss7_1 << "flowMinLimitOA_"<<v<<","<<i<<","<<t;
						#ifndef WaitAfterOperate
						flowCapacityOA[v][i][t] = IloRange(env, -IloInfinity, fOA[v][i][t]-oA[v][i][t]*inst.q_v[v], 0, ss7.str().c_str());
						if(tightenFlow){
							if(inst.typePort[i-1] == 1){ //Minimum flow in discharging ports
								flowMinCapacityOA[v][i][t] = IloRange(env, -IloInfinity, oA[v][i][t]*inst.f_min_jt[i-1][t-1]-fOA[v][i][t], 0, ss7_1.str().c_str());
								model.add(flowMinCapacityOA[v][i][t]);
							}else{ //Change upper bound for loading ports
								flowCapacityOA.setExpr(fOA[v][i][t] - oA[v][i][t]*(inst.q_v[v]-inst.f_min_jt[i-1][t-1]));
							}
						}
						#endif
						#ifdef WaitAfterOperate
						flowCapacityOA[v][i][t] = IloRange(env, -IloInfinity, fOA[v][i][t]-z[v][i][t]*inst.q_v[v], 0, ss7.str().c_str());
						if(tightenFlow){
							if(inst.typePort[i-1] == 1){ //Minimum flow in discharging port
								flowMinCapacityOA[v][i][t] = IloRange(env, -IloInfinity, z[v][i][t]*inst.f_min_jt[i-1][t-1]-fOA[v][i][t], 0, ss7_1.str().c_str());
								model.add(flowMinCapacityOA[v][i][t]);
							}else{ //Change upper bound for loading ports
								flowCapacityOA[v][i][t].setExpr(fOA[v][i][t] - z[v][i][t]*(inst.q_v[v]-inst.f_min_jt[i-1][t-1]));
							}
						}                    
						#endif
						model.add(flowCapacityOA[v][i][t]);
						
						if(t<tOEB && hasEnteringArc1st[v][i][t]==1){ //Constraints with no last time index
							#ifndef WaitAfterOperate
							ss8 << "flowLimitOB_"<<v<<","<<i<<","<<t;
							stringstream ss8_1;
							ss8 << "flowMinLimitOB_"<<v<<","<<i<<","<<t;
							flowCapacityOB[v][i][t] = IloRange(env, -IloInfinity, fOB[v][i][t]-oB[v][i][t](inst.q_v[v]), 0, ss8.str().c_str());
							model.add(flowCapacityOB[v][i][t]);
							if(tightenFlow){
								//Equal for discharging and loading ports
								flowMinCapacityOB[v][i][t] = IloRange(env, -IloInfinity, oB[v][i][t]*inst.f_min_jt[i-1][t-1] - fOB[v][i][t], 0, ss8_1.str().c_str());
								model.add(flowMinCapacityOB[v][i][t]);
								flowCapacityOB[v][i][t].setExpr(oB[v][i][t](inst.q_v[v]-inst.f_min_jt[i-1][t-1]) - fOB[v][i][t]);
							}
							#endif

							#ifdef WaitAfterOperate
							ss10 << "flowLimitWB_"<<v<<","<<i<<","<<t;
							stringstream ss10_1;
							ss10_1 << "flowMinLimitWB_"<<v<<","<<i<<","<<t;
							
							flowCapacityWB[v][i][t] = IloRange(env, -IloInfinity, fWB[v][i][t]-inst.q_v[v]*wB[v][i][t], 0, ss10.str().c_str());
							model.add(flowCapacityWB[v][i][t]);
							if(tightenFlow){
								//Equal for discharging and loading ports
								flowMinCapacityWB[v][i][t] = IloRange(env, -IloInfinity, wB[v][i][t]*inst.f_min_jt[i-1][t-1]-fWB[v][i][t], 0, ss10_1.str().c_str());
								model.add(flowMinCapacityWB[v][i][t]);
								flowCapacityWB[v][i][t].setExpr(wB[v][i][t]*(inst.q_v[v]-inst.f_min_jt[i-1][t-1])-fWB[v][i][t]);							
							}
							#endif
							ss9 << "flowLimitW_"<<v<<","<<i<<","<<t;
							stringstream ss9_1;
							ss9_1 << "flowLimitW_"<<v<<","<<i<<","<<t;
							flowCapacityW[v][i][t] = IloRange(env, -IloInfinity, fW[v][i][t]-inst.q_v[v]*w[v][i][t], 0, ss9.str().c_str());
							model.add(flowCapacityW[v][i][t]);
							if(tightenFlow){
								if(inst.typePort[i-1] == 1){ //Lower bound on flow for discharging ports
									flowMinCapacityW[v][i][t] = IloRange(env, -IloInfinity, w[v][i][t]*inst.f_min_jt[i-1][t-1]-fW[v][i][t], 0, ss9_1.str().c_str());
									model.add(flowMinCapacityW[v][i][t]);
								}else{ //Upper bound on loading ports
									flowCapacityW[v][i][t].setExpr(fW[v][i][t]-w[v][i][t]*(inst.q_v[v]-inst.f_min_jt[i-1][t-1]));
								}
							}
						}
                    }
                    #ifndef NBranching                    
                        expr_sumOA += oA[v][i][t];
                    #endif
				}
                flowCapacityX[v][i] = IloArray<IloRangeArray>(env,N);
                if(tightenFlow){
					flowMinCapacityX[v][i] = IloArray<IloRangeArray>(env,N);
				}	
                
                #ifndef NBranching                
                //~ priorityX[v][i] = IloRange(env, 0, expr_sumEnteringX - sumX[v][i], 0);
                //~ model.add(priorityX[v][i]);
                
                priorityOA[v][i] = IloRange(env, 0, expr_sumOA - sumOA[v][i], 0);
                model.add(priorityOA[v][i]);   
                #endif 
				for(j=1;j<N;j++){ //No arriving arc in source arc j=0
					if(i != j && i < N-1){ //There is no departing arc from sink node
						flowCapacityX[v][i][j] = IloRangeArray(env,T);
						if(tightenFlow){
							flowMinCapacityX[v][i][j] = IloRangeArray(env,T);
                        }
						for(t=0;t<=tOEB;t++){
							if(hasArc[v][i][j][t]==1){
								stringstream ss,ss1;
								ss << "flowLimitX_" << v << "," << i << "," << j << "," << t;
								ss1 << "flowMinLimitX_" << v << "," << i << "," << j << "," << t;
								if(j == N-1){ //If j is sink node                            
									flowCapacityX[v][i][j][t] = IloRange(env, -IloInfinity, fX[v][i][j][t] - inst.q_v[v]*x[v][i][j][t], 0, ss.str().c_str());
									model.add(flowCapacityX[v][i][j][t]);
								}else if (i>0 && t + inst.travelTime[v][i-1][j-1] <= tOEB){ //If x variable reach j in the time period
									//Original flow upper limit
									flowCapacityX[v][i][j][t] = IloRange(env, -IloInfinity, fX[v][i][j][t] - inst.q_v[v]*x[v][i][j][t], 0, ss.str().c_str());
									model.add(flowCapacityX[v][i][j][t]);
									if(tightenFlow){ //Thighter upper bound if traveling between 2 loading ports
										if(inst.typePort[i-1] == inst.typePort[j-1]){ //Only for cases where ports are of the same type
											if (inst.typePort[i-1] == 0){   //Loading ports 
												flowCapacityX[v][i][j][t].setExpr(fX[v][i][j][t] - x[v][i][j][t]*(inst.q_v[v]-inst.f_min_jt[j-1][0]));
												flowMinCapacityX[v][i][j][t] = IloRange(env, -IloInfinity, x[v][i][j][t]*inst.f_min_jt[i-1][0] - fX[v][i][j][t] , 0, ss1.str().c_str());
											}else{
												flowCapacityX[v][i][j][t].setExpr(fX[v][i][j][t] - x[v][i][j][t]*(inst.q_v[v]-inst.f_min_jt[i-1][0]));
												flowMinCapacityX[v][i][j][t] = IloRange(env, -IloInfinity, x[v][i][j][t]*inst.f_min_jt[j-1][0] - fX[v][i][j][t] , 0, ss1.str().c_str());
											}
											model.add(flowMinCapacityX[v][i][j][t]);
										}
									}
								}
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
    for(v=0;v<V;v++){
        for(i=1;i<=J;i++){
            //~ cplex.setPriority(sumX[v][i],2);
            cplex.setPriority(sumOA[v][i],1);
        }
    }
	#endif
	
	#ifndef NOperateAndGo
	
	#endif	
	
	#ifndef N2PortNorevisit

	#endif
    
    
}

/* Param p = 0 If algorithm can extract var solution values; 1 otherwise*/
void Model::fixSolution(IloEnv& env, Instance inst, const int& tS, const int& tF,const int& p, const bool& fixSinkArc){
	int J = inst.numTotalPorts;
	int N = J+1;
	int t0;
	
	if(p == 0){
		getSolValsW(env, inst, tS, tF, fixSinkArc);	//Get only values of the interval
	}	
	//Fix solution - If p == 1, consider the values from the previous p==0
	for(int v=0; v<inst.speed.getSize(); v++){
		for(int i=1; i<=J; i++){
            //Fixing values of last index of previous fixed block (if is not first iteration)
            if(tS>1){
                if(hasEnteringArc1st[v][i][tS-1]==1){
                    w[v][i][tS-1].setBounds(round(wValue[v][i][tS-1]), round(wValue[v][i][tS-1]));
                    #ifdef WaitAfterOperate
                    wB[v][i][tS-1].setBounds(round(wBValue[v][i][tS-1]), round(wBValue[v][i][tS-1]));
                    #endif
                    #ifndef WaitAfterOperate
                    oB[v][i][tS-1].setBounds(round(oBValue[v][i][tS-1]), round(oBValue[v][i][tS-1]));
                    #endif
                }
            }
			for(int j=1;j<=J;j++){ //Sink arc are considered separately
				t0 = max(1,tS-(int)inst.travelTime[v][i-1][j-1]);
				for(int t=t0; t<=tF; t++){
					if(i==1){	//(v,j,t) iterator
						if(hasEnteringArc1st[v][j][t]==1){
							if(t>=tS){ //Fixing values only of the current block                        
								z[v][j][t].setBounds(round(zValue[v][j][t]), round(zValue[v][j][t]));
								if(t<tF){
									w[v][j][t].setBounds(round(wValue[v][j][t]), round(wValue[v][j][t]));
									#ifdef WaitAfterOperate
									wB[v][j][t].setBounds(round(wBValue[v][j][t]), round(wBValue[v][j][t]));
									#endif
									#ifndef WaitAfterOperate
									oB[v][j][t].setBounds(round(oBValue[v][j][t]), round(oBValue[v][j][t]));
									#endif
								}
								#ifndef WaitAfterOperate
								oA[v][j][t].setBounds(round(oAValue[v][j][t]), round(oAValue[v][j][t]));
								#endif
                                if(fixSinkArc)
                                    x[v][j][N][t].setBounds(round(xValue[v][j][N][t]), round(xValue[v][j][N][t]));
							}
						}
					}
					if(i != j){ //For arcs
						if(hasArc[v][i][j][t]==1){
							int t2 = t + inst.travelTime[v][i-1][j-1];
							if (t2 >= tS && t2 <= tF) 
								x[v][i][j][t].setBounds(round(xValue[v][i][j][t]), round(xValue[v][i][j][t]));
						}
					}
				}
			}
		}
	}
}

void Model::reIntegralize(IloEnv& env, Instance inst, const int& tS, const int& tF){
	int i,j,t,v;
	int J = inst.numTotalPorts;	
	int V = inst.speed.getSize(); //# of vessels
    int N = J + 1;
	
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
            
            for(j=1;j<=N;j++){
                int t0;
                if (j<=J)
                    t0 = max(1, tS-(int)inst.travelTime[v][i-1][j-1]);
                else
                    t0 = tS;
                for(t=t0; t<=tF; t++){
                    ///Vessel port time (v,i,t)
                    if(j==1){
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
                    }
                    ///End vessel port time (v,i,t)
                    if(j!=i){
                        if (j < N){   		//If arrives at a port AND
                            int t2 = t + inst.travelTime[v][i-1][j-1];
                            if(t2 >= tS && t2 <= tF){ //Traveling arc is in the integer block OR crosses previos integer block and the integer block
                                model.remove(convertX[v][i][j][t]);
                            }
                        }else{ //Sink arc
                            if (t >= tS)
                                model.remove(convertX[v][i][j][t]);
                        }
                    }
                }
            }            
		}
	}
}

void Model::getSolValsW(IloEnv& env, Instance inst, const int& tS, const int& tF, const bool& fixSinkArc){
	if((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){
		int v, i, j, t, t0;
        int V = inst.speed.getSize();
        int J = inst.numTotalPorts;
        int N = J+1;
        for(v=0;v<V;v++){
            for(i=1;i<=J;i++){ //Not needed to consider source node either sink node
                for (j=1;j<=J;j++){ //Consider only ports, as the sink arcs are considered separately
					t0 = max(1,tS-(int)inst.travelTime[v][i-1][j-1]);
					for(t=t0;t<=tF;t++){
						if(hasArc[v][i][j][t]==1){
							int t2 = t + inst.travelTime[v][i-1][j-1];
							if (t2>=tS && t2<=tF) //If is traveling arc is in the integer block or crossing a fixed block and integer block
								xValue[v][i][j][t] = cplex.getValue(x[v][i][j][t]); 
						}
						if (i == 1){ //(v,j,t) iterator
							if(hasEnteringArc1st[v][j][t] == 1) {
								zValue[v][j][t] = cplex.getValue(z[v][j][t]);
								//~ //~ //~ cout "z_["<<v<<"]["<<i<<"]["<<t<<"] = " << zValue[v][j][t] << endl;
								
								#ifndef WaitAfterOperate
								oAValue[v][j][t] = cplex.getValue(oA[v][j][t]);
								//~ //~ //~ cout "oA_["<<v<<"]["<<i<<"]["<<t<<"] = " << oAValue[v][j][t] << endl;
								#endif
								
								if(t<tF){
									#ifndef WaitAfterOperate
									oBValue[v][j][t] = cplex.getValue(oB[v][j][t]);
									//~ //~ //~ cout "oB_["<<v<<"]["<<i<<"]["<<t<<"] = " << oBValue[v][j][t] << endl;
									#endif
									#ifdef WaitAfterOperate
									wBValue[v][j][t] = cplex.getValue(wB[v][j][t]);
									//~ //~ //~ cout "wB_["<<v<<"]["<<i<<"]["<<t<<"] = " << wBValue[v][j][t] << endl;
									#endif                            
									wValue[v][j][t] = cplex.getValue(w[v][j][t]);
									//~ //~ //~ cout "w_["<<v<<"]["<<i<<"]["<<t<<"] = " << wValue[v][j][t] << endl;
								}								
								if(fixSinkArc){
                                    xValue[v][j][N][t] = cplex.getValue(x[v][j][N][t]);  //Gets the value of sink arc
                                    //~ //~ //~ cout "x_["<<v<<"]["<<j<<"]["<<N<<"]["<<t<<"] = " << xValue[v][j][N][t] << endl;
                                    //~ //~ //~ cout "fx_["<<v<<"]["<<j<<"]["<<N<<"]["<<t<<"] = " << cplex.getValue(fX[v][j][N][t]) << endl;
                                }
							}
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
                        //~ //~ //~ cout "wb["<<v<<"]["<<i<<"["<<tS-1<<"] = " << wBValue[v][i][tS-1] << endl;
						#endif
                        wValue[v][i][tS-1] = cplex.getValue(w[v][i][tS-1]);                        
                        //~ //~ //~ cout "w["<<v<<"]["<<i<<"["<<tS-1<<"] = " << wValue[v][i][tS-1] << endl;
					}
                }  
                //Gets the node inventory - for feeding local valid inequalities
                for(t=1;t<=tF;t++){
					sPValue[i][t] = cplex.getValue(sP[i][t]);					  
					//~ cout << "sP[" << i << "," << t << "] " << sPValue[i][t] << endl;
				}
            }
        }
	}else{
		cout << "Impossible to get feasible solution values!" << endl;
		exit(1);
	}
}
void Model::getSolution(IloEnv& env, Instance inst){
	if((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){
		int v, i, j, t, t0;
        int V = inst.speed.getSize();
        int J = inst.numTotalPorts;
        int N = J+1;
        for(v=0;v<V;v++){
            for(i=1;i<=J;i++){ //Not needed to consider source node either sink node
                for (j=1;j<=J;j++){ //Consider only ports, as the sink arcs are considered separately					
					for(t=1;t<=inst.t;t++){
						if(hasArc[v][i][j][t]==1){
							xValue[v][i][j][t] = cplex.getValue(x[v][i][j][t]); 
						}
						if (i == 1){ //(v,j,t) iterator
							if(hasEnteringArc1st[v][j][t] == 1) {
								zValue[v][j][t] = cplex.getValue(z[v][j][t]);
							
								#ifndef WaitAfterOperate
								oAValue[v][j][t] = cplex.getValue(oA[v][j][t]);
								#endif
								
								if(t<inst.t){
									#ifndef WaitAfterOperate
									oBValue[v][j][t] = cplex.getValue(oB[v][j][t]);
									#endif
									#ifdef WaitAfterOperate
									wBValue[v][j][t] = cplex.getValue(wB[v][j][t]);
									#endif     
									wValue[v][j][t] = cplex.getValue(w[v][j][t]);
								}				
		                         xValue[v][j][N][t] = cplex.getValue(x[v][j][N][t]);  //Gets the value of sink arc
							}
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
void Model::getSolutionVesselPair(IloEnv& env, Instance inst, const unsigned int& v1, const unsigned int& v2){
	if((cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible)){
		unsigned int v, idV, i, j, t, t0;
		array<unsigned int,2> VIDs = {v1,v2};
        
        int J = inst.numTotalPorts;
        int N = J+1;
        for(idV=0;idV<VIDs.size();idV++){
			v = VIDs[idV];
            for(i=1;i<=J;i++){ //Not needed to consider source node either sink node
                for (j=1;j<=J;j++){ //Consider only ports, as the sink arcs are considered separately					
					for(t=1;t<=inst.t;t++){
						if(hasArc[v][i][j][t]==1){
							xValue[v][i][j][t] = cplex.getValue(x[v][i][j][t]); 
						}
						if (i == 1){ //(v,j,t) iterator
							if(hasEnteringArc1st[v][j][t] == 1) {
								zValue[v][j][t] = cplex.getValue(z[v][j][t]);
							
								#ifndef WaitAfterOperate
								oAValue[v][j][t] = cplex.getValue(oA[v][j][t]);
								#endif
								
								if(t<inst.t){
									#ifndef WaitAfterOperate
									oBValue[v][j][t] = cplex.getValue(oB[v][j][t]);
									#endif
									#ifdef WaitAfterOperate
									wBValue[v][j][t] = cplex.getValue(wB[v][j][t]);
									#endif                            
									wValue[v][j][t] = cplex.getValue(w[v][j][t]);
								}								
                                xValue[v][j][N][t] = cplex.getValue(x[v][j][N][t]);  //Gets the value of sink arc
							}
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

void Model::modifyModel(IloEnv& env, Instance inst, const int& nIntervals, const int& tS_fix, const int& tF_fix, 
    const int& tS_add, const int& tF_add, const int& tS_int, const int& tF_int, const bool& validIneq,
    const bool& addConstr, const bool& tightenInvConstr, const bool& proportionalAlpha, const bool& tightenFlow){
	int i,j,t,v,a, t0;
	int T = inst.t;
	//~ int timePerIntervalId = T/nIntervals; //Number of time periods in each interval
	int timePerIntervalId = nIntervals; //Number of time periods in each interval
	int J = inst.numTotalPorts;	
	double intPart;	
	int V = inst.speed.getSize(); //# of vessels
    int N = J + 1;    	
	//First gets the solution values for fixing (when needed)
	if (tS_fix < T){		
		cout <<  "Fixing interval [" << tS_fix << "," << tF_fix << "]\n";
		getSolValsW(env, inst, tS_fix, tF_fix, false);	//Get only values of the interval for fixing
	}
	
	if (tS_int <= T){
		cout << "Integralizing: [" << tS_int << "," << tF_int << "]\n";
	}	
	///Removing end-block
	if(tS_add <= T){
		cout << "Adding to the model [" << tS_add << "," << tF_add << "]\n";
	}
		
    ///Objective function	
	IloExpr expr_obj_current = cplex.getObjective().getExpr(); //Gets the current objective function
	IloExpr expr_obj_cost(env);
	IloExpr expr_obj_revenue(env);
			
	#ifndef NAlpha
	IloExpr expr_cumSlack(env);
	#endif
	IloExpr expr_berth(env), expr_invBalancePort(env);    
	IloExpr expr_sinkFlow(env), expr_1stLevel(env), expr_2ndLevel(env), expr_1stFlow(env), expr_2ndFlow(env);
	IloExpr expr_opd(env);
    
    
    #ifndef NWWCCReformulation
    IloExpr expr_sumF(env);
    IloExpr expr_sumO(env);
    IloExpr expr_wwcc(env);
    int it_kt;
    #endif
    
	for(v=0;v<V;v++){
		if(tS_add <= T){
			expr_sinkFlow.clear();
			expr_sinkFlow = sinkNodeBalance[v].getExpr();
		}
		for(i=1;i<=J;i++){			
			//Fixing values of last index of previous fixed block (if is not first iteration)
            if(tS_fix>1){
                if(hasEnteringArc1st[v][i][tS_fix-1]==1){
                    w[v][i][tS_fix-1].setBounds(round(wValue[v][i][tS_fix-1]), round(wValue[v][i][tS_fix-1]));
                    #ifdef WaitAfterOperate
                    wB[v][i][tS_fix-1].setBounds(round(wBValue[v][i][tS_fix-1]), round(wBValue[v][i][tS_fix-1]));
                    #endif
                    #ifndef WaitAfterOperate
                    oB[v][i][tS_fix-1].setBounds(round(oBValue[v][i][tS_fix-1]), round(oBValue[v][i][tS_fix-1]));
                    #endif
                }
            }
            
			//Integralizing the last index of previous block
			if (tS_int <= T){
				//~ if(hasEnteringArc1st[v][i][tS_int-1]==1){
					model.remove(convertW[v][i][tS_int-1]);
					#ifndef WaitAfterOperate
					model.remove(convertOB[v][i][tS_int-1]);
					#endif
					#ifdef WaitAfterOperate
					model.remove(convertWB[v][i][tS_int-1]);
					#endif
				//~ }
			}
            
            ///Updating flows and limits of last time period of previous interval
			if(tS_add <= T){
				expr_1stLevel.clear();
				expr_2ndLevel.clear();
				expr_1stFlow.clear();
				expr_2ndFlow.clear();
				
				if(hasEnteringArc1st[v][i][tS_add-1]==1){ //Only necessary if there is entering arc in the node
					expr_1stLevel = firstLevelBalance[v][i][tS_add-1].getExpr();
					expr_2ndLevel = secondLevelBalance[v][i][tS_add-1].getExpr();
					expr_1stFlow = firstLevelFlow[v][i][tS_add-1].getExpr();
					expr_2ndFlow = secondLevelFlow[v][i][tS_add-1].getExpr();
					
					expr_1stLevel += - w[v][i][tS_add-1];
					expr_1stFlow += -fW[v][i][tS_add-1];
					
					#ifndef WaitAfterOperate                
					expr_2ndLevel += oB[v][i][tS_add-1];
					expr_2ndFlow += fOB[v][i][tS_add-1];
					#endif
					
					#ifdef WaitAfterOperate
					expr_2ndLevel += wB[v][i][tS_add-1];
					expr_2ndFlow += fWB[v][i][tS_add-1];
					#endif
					
					firstLevelBalance[v][i][tS_add-1].setExpr(expr_1stLevel);
					firstLevelFlow[v][i][tS_add-1].setExpr(expr_1stFlow);
					secondLevelBalance[v][i][tS_add-1].setExpr(expr_2ndLevel);
					secondLevelFlow[v][i][tS_add-1].setExpr(expr_2ndFlow);
				
					//Add constraints on the flow variables of last time period of previous block.
					stringstream ss8,ss9,ss10;
					#ifndef WaitAfterOperate
					ss8 << "flowLimitOB_"<<v<<","<<i<<","<<tS_add-1;
					flowCapacityOB[v][i][tS_add-1] = IloRange(env, -IloInfinity, fOB[v][i][tS_add-1]-inst.q_v[v]*oB[v][i][tS_add-1], 0, ss8.str().c_str());
					model.add(flowCapacityOB[v][i][tS_add-1]);                            
					#endif
					ss9 << "flowLimitW_"<<v<<","<<i<<","<<tS_add-1;
					flowCapacityW[v][i][tS_add-1] = IloRange(env, -IloInfinity, fW[v][i][tS_add-1]-inst.q_v[v]*w[v][i][tS_add-1], 0, ss9.str().c_str());
					model.add(flowCapacityW[v][i][tS_add-1]);
					#ifdef WaitAfterOperate
					ss10 << "flowLimitWB_"<<v<<","<<i<<","<<tS_add-1;
					flowCapacityWB[v][i][tS_add-1] = IloRange(env, -IloInfinity, fWB[v][i][tS_add-1]-inst.q_v[v]*wB[v][i][tS_add-1], 0, ss10.str().c_str());
					model.add(flowCapacityWB[v][i][tS_add-1]);
					#endif				
				}

			}
			///end updating    
						
			///port iterator (i)
			if(v==0){
				if(tS_add <= T){
					#ifndef NAlpha
					expr_cumSlack.clear();
					expr_cumSlack = cumSlack[i].getExpr();
					#endif
				}						
					if(validIneq){ 
						int it,it2=0,k,l;
						IloExpr expr_kP1_LHS(env), expr_kP2_LHS(env), expr_kD1_LHS(env), expr_kD2_LHS(env);
						for(it=0;it < ((T-2)*2 - (T-tF_add));it++){	//'Revise' existing inequalities (x variables) and add the new until [tF-1, tF]
							double kP1_RHS=0, kP2_RHS=0, kD1_RHS=0, kD2_RHS=0, sum_alphaMax=0, alphaUB=0;
							expr_kP1_LHS.clear();
							expr_kP2_LHS.clear();
							expr_kD1_LHS.clear();
							expr_kD2_LHS.clear();
							//Definining the size of set T =[l,k]
							l=1;
							k=tF_add;
							if(it<tF_add-2){ //First part of array - varying k
								k = it+2;
								it2++;
							}else if(it==tF_add-2){ //If reaches the 'limit' of the model - Starting to vary l
								if(tF_add-2 < T-3) //If it is not the last part of adding the model (otherwise it does not need to be changed)
									it = T-3;
								l = 2 + (it-(T-2));
								it2++;
							}else{ // when it is > T-2
								l = 2 + (it-(T-2));
							}  
							//valid inequality is created from scratch - If it is in interval [1...[tS-2..tF-1]], or If it is in interval [[2..tS-1]...tF]
							if( (it >= tS_add-3 && it < tF_add-2) || 
								( (it >= T-2 + tS_add-5 && it < T-2+tF_add-2) || (tF_add == T && it > T-2 + tS_add-5) ) ||
								  (it >= T-2 && it < T-2 + tS_add-4)){ //Case of extending which is equivalent to add
								//~ if(i==1){
									//~ cout << "Adding new: " << it << " [" << l << "," << k << "]";
									//~ if (it >= T-2 && it < T-2-1 + tS_add-2) //When updating the VI1
										//~ //~ //~ cout " - Replacing";
									//~ //~ //~ cout endl;
								//~ }
								for(int v1=0;v1<V;v1++){
									#ifdef WaitAfterOperate
									if(hasEnteringArc1st[v1][i][k]==1 && k < tF_add){
										expr_kP1_LHS += wB[v1][i][k];
									}
									if(l>1 && hasEnteringArc1st[v1][i][l-1]==1)      
										expr_kD1_LHS += w[v1][i][l-1] + wB[v1][i][l-1];
									#endif                                
									#ifndef WaitAfterOperate
									if(hasEnteringArc1st[v1][i][k]==1 && k < tF_add){
										expr_kP1_LHS += oB[v1][i][k];
									}
									if(l>1 && hasEnteringArc1st[v1][i][l-1]==1){
										expr_kD1_LHS += w[v1][i][l-1] + oB[v1][i][l-1];
										expr_kD2_LHS += oB[v1][i][l-1];
									}
									#endif       
									for(t=l;t<=k;t++){
										if(hasEnteringArc1st[v1][i][t]==1){                       
											expr_kP2_LHS += z[v1][i][t];	//It is used for both loading and discharging ports
											#ifndef WaitAfterOperate
											expr_kD2_LHS += oA[v1][i][t];
											#endif
											#ifdef WaitAfterOperate
											expr_kD2_LHS += z[v1][i][t];
											#endif
										}
										for(j=0;j<=N;j++){						
											if(j==0 && inst.initialPort[v1]+1 == i && t==inst.firstTimeAv[v1]+1) 	//Source arc
												expr_kD1_LHS += x[v1][j][i][0];								
											else if (j == N && hasArc[v1][i][j][t] == 1){						 //Sink arc
												expr_kP1_LHS += x[v1][i][j][t];
											}
											else if (j > 0 && j <= J){											//"Normal" arcs
												if(i != j){
													//~ if(hasArc[v][i][j][t] == 1)  //Ignoring the rule that the arriving node must be in the model
													if(hasArc[v1][i][j][t] == 1 && t+inst.travelTime[v1][i-1][j-1] <= tF_add){  //If arc exists and arrives at port j in the integer or relaxed block
														expr_kP1_LHS += x[v1][i][j][t];														
													}
													if(t - inst.travelTime[v1][j-1][i-1] > 0){					//Only if it is possible an entering arc due the travel time
														if(hasArc[v1][j][i][t-inst.travelTime[v1][j-1][i-1]] == 1)
															expr_kD1_LHS += x[v1][j][i][t-inst.travelTime[v1][j-1][i-1]];
													}
												}
											}
										}
										//Port-time loop (i,t)
										if(v1==0){
											kP1_RHS += inst.d_jt[i-1][t-1];
											kD1_RHS += inst.d_jt[i-1][t-1];
											#ifndef NAlpha
											sum_alphaMax += inst.alp_max_jt[i-1][t-1];
											#endif
										}
									}
								}
								if(l==1){
									kP1_RHS += -inst.sMax_jt[i-1][0] + inst.s_j0[i-1];
									kD1_RHS += -inst.s_j0[i-1] + inst.sMin_jt[i-1][0];
								}else{
									kP1_RHS += -inst.sMax_jt[i-1][0];
									kD1_RHS += inst.sMin_jt[i-1][0];
									if(l <= tS_fix){ //Local valid inequality
										kP1_RHS += sPValue[i][l-1];
										kD1_RHS += -sPValue[i][l-1];
									}else{
										kP1_RHS += inst.sMin_jt[i-1][0];
										kD1_RHS += -inst.sMax_jt[i-1][0];
									}
								}                                
								///If considering alpha parameters, otherwise comment the above 2 lines
								#ifndef NAlpha
								kP1_RHS += - min(sum_alphaMax, inst.alp_max_j[i-1]);
								kD1_RHS += - min(sum_alphaMax, inst.alp_max_j[i-1]);
								#endif
								
								kP2_RHS = max(0.0,ceil(kP1_RHS/inst.f_max_jt[i-1][0]));
								kP1_RHS = max(0.0,ceil(kP1_RHS/inst.maxVesselCapacity));
								kD2_RHS = max(0.0,ceil(kD1_RHS/inst.f_max_jt[i-1][0]));
								kD1_RHS = max(0.0,ceil(kD1_RHS/inst.maxVesselCapacity));
								
								stringstream ss, ss1, ss2;
								if(inst.typePort[i-1] == 0){
									ss << "knpasackP1_" << i << "_(" << l << "," << k << ")";
									ss1 << "knapsackP2_" << i << "_(" << l << "," << k << ")";
									if (it >= T-2 && it < T-2-1 + tS_add-2){ //When updating the VI1
										//~ cout << "RHS = (P) " << kP1_RHS << endl;
										knapsack_P_1[i][it].setBounds(kP1_RHS,IloInfinity);
										knapsack_P_1[i][it].setExpr(expr_kP1_LHS);
										knapsack_P_1[i][it].setName(ss.str().c_str());
										knapsack_P_2[i][it].setBounds(kP2_RHS,IloInfinity);
										knapsack_P_2[i][it].setExpr(expr_kP2_LHS);
										knapsack_P_2[i][it].setName(ss1.str().c_str());
									}else{ //When adding a new										
										knapsack_P_1[i][it] = IloRange(env, kP1_RHS, expr_kP1_LHS, IloInfinity, ss.str().c_str());
										knapsack_P_2[i][it] = IloRange(env, kP2_RHS, expr_kP2_LHS, IloInfinity, ss1.str().c_str());
										model.add(knapsack_P_1[i][it]);
										model.add(knapsack_P_2[i][it]);
									}
								}else{
									ss << "knpasackD1_" << i << "_(" << l << "," << k << ")";
									ss1 << "knpasackD2_" << i << "_(" << l << "," << k << ")";
									knapsack_D_1[i][it] = IloRange(env, kD1_RHS, expr_kD1_LHS, IloInfinity, ss.str().c_str());
									if (it >= T-2 && it < T-2-1 + tS_add-2){ //When updating the VI1
										knapsack_D_2[i][it].setExpr(expr_kD2_LHS);
										knapsack_D_2[i][it].setName(ss1.str().c_str());
										#ifndef WaitAfterOperate                                    
										knapsack_D_2[i][it].setBounds(kD1_RHS, IloInfinity);
										#endif
										#ifdef WaitAfterOperate                                    
										knapsack_D_2[i][it].setBounds(kD2_RHS, IloInfinity);                                    
										#endif
										#ifndef WaitAfterOperate
										ss2 << "knpasackD3_" << i << "_(" << l << "," << k << ")";
										knapsack_D_3[i][it].setBounds(kD2_RHS, IloInfinity);
										knapsack_D_3[i][it].setExpr(expr_kP2_LHS);
										knapsack_D_3[i][it].setName(ss2.str().c_str());                                    
										#endif
									}else{ //When adding a new
										#ifndef WaitAfterOperate
										knapsack_D_2[i][it] = IloRange(env, kD1_RHS, expr_kD2_LHS, IloInfinity, ss1.str().c_str());
										#endif
										#ifdef WaitAfterOperate
										knapsack_D_2[i][it] = IloRange(env, kD2_RHS, expr_kD2_LHS, IloInfinity, ss1.str().c_str());
										#endif
										model.add(knapsack_D_1[i][it]);
										model.add(knapsack_D_2[i][it]);
										#ifndef WaitAfterOperate
										ss2 << "knpasackD3_" << i << "_(" << l << "," << k << ")";
										knapsack_D_3[i][it] = IloRange(env, kD2_RHS, expr_kP2_LHS, IloInfinity, ss2.str().c_str());
										model.add(knapsack_D_3[i][it]);
										#endif
									}
								}
							}
							else if (it >= 0 && it < tS_add-3){		//Only modify the previously added valid inequalities of type [1..j] (only needed add x variables in the case of loading ports)
								int t0 = max(1, tS_add-inst.maxTravelTimeInstance);
								if(k >= t0){
									//~ if(i==1)
										//~ cout << "Updating  " << it << " ["<< l << "," << k << "]\n";
									//Get the current expr values                 
									if (inst.typePort[i-1] == 0)
										expr_kP1_LHS = knapsack_P_1[i][it].getExpr();
									for(int v1=0;v1<V;v1++){										
										for(t=l;t<=k;t++){
											for(j=1;j<=J;j++){													
												//"Normal" arcs	
												if(i != j){
													int t2 = t+inst.travelTime[v1][i-1][j-1];
													if(hasArc[v1][i][j][t] == 1 && (t2>=tS_add && t2<=tF_add) )  //If arc exists, was not added in the previous iteration and arrives at port j in the integer or relaxed block
														expr_kP1_LHS += x[v1][i][j][t];												
												}										
											}
										}
									}
									if(inst.typePort[i-1] == 0)
										knapsack_P_1[i][it].setExpr(expr_kP1_LHS);												
								}
							}						
						}
					}	
				//~ }
                #ifndef NWWCCReformulation
                it_kt= (tS_add-1)*tS_add/2;
                #endif
			}
						
			for(j=1;j<=N;j++){//Considering ports and sink node								 
				if(j<=J)
					t0 = max(1,tS_fix-(int)inst.travelTime[v][i-1][j-1]);
				else
					t0 = tS_fix;
				float thighetInventoryValue = 0;
				for(t=t0;t<=tF_add;t++){
					///Port-time iterator (i,t)
					if(v==0 && j==1){
						//For decreaseEndBlock
						if(tS_add <= T){
							if(t>=tS_add){
								expr_berth.clear();
								expr_invBalancePort.clear();
								stringstream ss, ss2;
								ss << "berthLimit_(" << i << "," << t << ")";
								bool emptyExpr = true;
								for(int v1=0;v1<V;v1++){
									if (hasEnteringArc1st[v1][i][t]==1){ //Only if exists an entering arc in the node
										expr_berth += z[v1][i][t];
										emptyExpr = false;
									
										expr_invBalancePort += -f[v1][i][t];
									}
								}
								if(!emptyExpr){ //Only if there are some Z variable in the expr
									berthLimit[i][t] = IloRange(env,-IloInfinity, expr_berth, inst.b_j[i-1], ss.str().c_str());
									model.add(berthLimit[i][t]);
								}
								
								#ifndef NAlpha
								expr_cumSlack += alpha[i][t];
								expr_invBalancePort += -alpha[i][t];
								#endif
								#ifndef NBetas
								expr_invBalancePort += beta[i][t];
								#endif
                                #ifndef NThetas
                                expr_invBalancePort += -theta[i][t];
								#endif
								ss2 << "invBalancePort_(" << i << "," << t << ")";
								portInventory[i][t] = IloRange(env, inst.delta[i-1]*inst.d_jt[i-1][t-1],
									sP[i][t]-sP[i][t-1]-inst.delta[i-1]*expr_invBalancePort,
									inst.delta[i-1]*inst.d_jt[i-1][t-1], ss2.str().c_str());
								model.add(portInventory[i][t]);
								
								//Obj
								#ifndef NAlpha
								expr_obj_cost += inst.p_jt[i-1][t-1]*alpha[i][t];								//4rd term
								#endif
								#ifndef NBetas
								expr_obj_cost += PENALIZATION*beta[i][t];										//Auxiliary variables
								#endif
                                #ifndef NThetas
                                expr_obj_cost += PENALIZATION*theta[i][t];										//Auxiliary variables
								#endif
                                
                                #ifndef NWWCCReformulation
                                expr_sumF.clear();
                                expr_sumO.clear();
                                #endif
							}
						}
					}
					///end Port-time iterator
					
					///Vessel port time iterator (v,i,t)
					if(j==1){
						//For fixing
						if(t<=tF_fix){
							if(hasEnteringArc1st[v][i][t]==1){
								if(t>=tS_fix){ //Fixing values only of the current block                        
									z[v][i][t].setBounds(round(zValue[v][i][t]), round(zValue[v][i][t]));
									if(t<tF_fix){
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
									
									#ifndef NFixSinkArc
									x[v][i][N][t].setBounds(round(xValue[v][i][N][t]), round(xValue[v][i][N][t]));
									#endif
								}							
							}
						}
						
						//For integralizing block
						if(tS_int<=T){
							//~ if(hasEnteringArc1st[v][i][t]==1){
								if(t>=tS_int && t<=tF_int){
									model.remove(convertZ[v][i][t]);
									#ifndef WaitAfterOperate
									model.remove(convertOA[v][i][t]);
									#endif
									if(t<tF_int){
										model.remove(convertW[v][i][t]);
										#ifndef WaitAfterOperate
										model.remove(convertOB[v][i][t]);
										#endif
										#ifdef WaitAfterOperate
										model.remove(convertWB[v][i][t]);
										#endif
									}
								}
							//~ }
						}
						//For decreaseEndBlock
						if(tS_add <= T){
							if(t>=tS_add){ //New components
								if(hasEnteringArc1st[v][i][t]==1){
									expr_obj_revenue += inst.r_jt[i-1][t-1]*f[v][i][t];						//1st term
									expr_obj_cost += (t-1)*inst.attemptCost*z[v][i][t];						//3rd term
                                    #ifndef NWWCCReformulation
                                    expr_sumF += f[v][i][t];
                                    expr_sumO += z[v][i][t];
                                    #endif
								}
							}
						}						
					}					
					///end vessel port time iterator (v,i,t)
					stringstream ss, ss1, ss2, ss3;
					ss2 << "flowLimitX_" << v << "," << i << "," << j << "," << t;
					ss3 << "flowMinLimitX_" << v << "," << i << "," << j << "," << t;
					if (j<=J){				//When j is a port
						int t2 = t + inst.travelTime[v][i-1][j-1];
						//For decreaseEndBlock
						if(tS_add <= T){
							if(t2>=tS_add && t2<=tF_add){ // Arc in the new added block or in the intersection with new block and previous interval
								if(hasArc[v][i][j][t] == 1){
									//Obj - traveling arcs
									expr_obj_cost += arcCost[v][i][j][t]*x[v][i][j][t];	
									//~ cout << x[v][i][j][t].getName() << "Travel Cost " << arcCost[v][i][j][t] << endl;
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
									if (i != j){									
										//Normal upper limit flow
										flowCapacityX[v][i][j][t] = IloRange(env, -IloInfinity, fX[v][i][j][t] - inst.q_v[v]*x[v][i][j][t], 0, ss2.str().c_str());
										if(tightenFlow){ //Thighter upper bound if traveling between 2 loading ports
											if(inst.typePort[i-1] == inst.typePort[j-1]){ //Only for cases where ports are of the same type
												if (inst.typePort[i-1] == 0){   //Loading ports 
													flowCapacityX[v][i][j][t].setExpr(fX[v][i][j][t] - x[v][i][j][t]*(inst.q_v[v]-inst.f_min_jt[j-1][0]));
													flowMinCapacityX[v][i][j][t] = IloRange(env, -IloInfinity, x[v][i][j][t]*inst.f_min_jt[i-1][0] - fX[v][i][j][t] , 0, ss1.str().c_str());
												}else{  // Discharging ports
													flowCapacityX[v][i][j][t].setExpr(fX[v][i][j][t] - x[v][i][j][t]*(inst.q_v[v]-inst.f_min_jt[i-1][0]));
													flowMinCapacityX[v][i][j][t] = IloRange(env, -IloInfinity, x[v][i][j][t]*inst.f_min_jt[j-1][0] - fX[v][i][j][t] , 0, ss1.str().c_str());
												}
												model.add(flowMinCapacityX[v][i][j][t]);
											}
										}
										model.add(flowCapacityX[v][i][j][t]);
									}
								}
							}
						}
						
						//For integralizingBlock  (including not existing arcs, hasArc[v][i][j][t]=0)
						if(hasArc[v][i][j][t] == 1){
							if(tS_int <= T){
								if(t2>=tS_int && t2<=tF_int){
									model.remove(convertX[v][i][j][t]);
								}
							}
						}
						//For fixing block
						if(hasArc[v][i][j][t] == 1){
							if (t2>= tS_fix && t2 <= tF_fix){
								x[v][i][j][t].setBounds(round(xValue[v][i][j][t]), round(xValue[v][i][j][t]));
							}
						}
						
					}else{ // when j is the sink node
						//For decreaseEndBlock
						if(tS_add <= T){
							if(t >= tS_add){ //Only for nodes in the current added block
								if(hasArc[v][i][j][t]==1){
									//Obj
									expr_obj_cost += arcCost[v][i][j][t]*x[v][i][j][t];
									//~ cout << x[v][i][j][t].getName() << "Sink Cost " << arcCost[v][i][j][t] << endl;
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
									
									flowCapacityX[v][i][j][t] = IloRange(env, -IloInfinity, fX[v][i][j][t] - inst.q_v[v]*x[v][i][j][t], 0, ss2.str().c_str());
									model.add(flowCapacityX[v][i][j][t]);
								}
							}
						}
						//For integralizingBlock (including not existing arcs, hasArc[v][i][j][t]=0)
						if(hasArc[v][i][j][t] == 1){
							if(tS_int <= T){
								if (t >= tS_int && t <= tF_int){
									model.remove(convertX[v][i][j][t]);
								}
							}
						}
					}
				}
				if(tightenInvConstr){
					if(v==0 && j == 1 && tS_add <= T){
						//Update the previous thighted inventory and thight the new last interval
						if(inst.delta[i-1]==1){
							sP[i][tS_add-1].setUB(inst.sMax_jt[i-1][0]);							
							if(tF_add < T){
								sP[i][tF_add].setUB(inst.sMax_jt[i-1][0]-inst.sMax_jt[i-1][0]*0.1);
							}
						}else{
							sP[i][tS_add-1].setLB(inst.sMin_jt[i-1][0]);
							if(tF_add < T){
								sP[i][tF_add].setLB(inst.sMin_jt[i-1][0]+inst.sMax_jt[i-1][0]*0.1);
							}
						}
					}
				}
			}
			///Need another time iterator (v,i,t) - only for decrease endBlock
			if(tS_add <= T){
                #ifndef NBranching
                //~ IloExpr expr_sumEnteringX = priorityX[v][i].getExpr();
                IloExpr expr_sumOA = priorityOA[v][i].getExpr();
                #endif
				t0 = max(1, tS_add-(int)inst.max_travelTime[v]);
				for(t=t0;t<=tF_add;t++){
					if(t<tS_add){   ///Needed for update outgoing flow balance arcs (just 2nd level)
						expr_2ndLevel = secondLevelBalance[v][i][t].getExpr();
						expr_2ndFlow = secondLevelFlow[v][i][t].getExpr();
						if(addConstr){
							if(hasEnteringArc1st[v][i][t]==1 && inst.q_v[v] <= inst.f_max_jt[i-1][0]){
								expr_opd.clear();
								expr_opd = operateAndDepart[v][i][t].getExpr();
							}
						}
						for(int j1=1;j1<=J;j1++){
							if (hasArc[v][i][j1][t] == 1){
								if (t + inst.travelTime[v][i-1][j1-1] >= tS_add && t + inst.travelTime[v][i-1][j1-1] <= tF_add){ //If arc departs from previos interval (t<tS) and arrive in the current interval
									expr_2ndLevel += x[v][i][j1][t];
									expr_2ndFlow += fX[v][i][j1][t];
									if(addConstr){
										if (inst.q_v[v] <= inst.f_max_jt[i-1][0] && inst.typePort[i-1] != inst.typePort[j1-1])
											expr_opd += x[v][i][j1][t];
									}
								}
							}
						}
						//Updating previous flow 2nd level
						secondLevelBalance[v][i][t].setExpr(expr_2ndLevel);
						secondLevelFlow[v][i][t].setExpr(expr_2ndFlow);
						if(addConstr){
							if (hasEnteringArc1st[v][i][t]==1 && inst.q_v[v] <= inst.f_max_jt[i-1][0]){							
								operateAndDepart[v][i][t].setExpr(expr_opd);							
							}
						}
					}else{ /// if t>= tS_add : Only considering the time-periods of the added interval
						stringstream ss, ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8, ss9, ss10, ss11;
						ss << "First_level_balance_" << v << "(" << i << "," << t << ")";
						ss1 << "Second_level_balance_" << v << "(" << i << "," << t << ")";
						ss2 << "link_balance_" << v << "(" << i << "," << t << ")";
						ss3 << "First_level_flow_" << v << "(" << i << "," << t << ")";
						ss4 << "Second_level_flow_" << v << "(" << i << "," << t << ")";
						
						if(hasArc[v][i][N][t]==1){
							expr_sinkFlow += x[v][i][N][t];
						}
						
						expr_1stLevel.clear();
						expr_1stFlow.clear();
						expr_2ndLevel.clear();
						expr_2ndFlow.clear();
						if(addConstr) expr_opd.clear();
						
						for(int j1=0;j1<=N;j1++){
							if(j1> 0 && j1<N){ //No consider sink arc (first level balance)
								if (t - inst.travelTime[v][j1-1][i-1] >= 0){ //If it is possible to exist an arc from j1 to i
									if(hasArc[v][j1][i][t-inst.travelTime[v][j1-1][i-1]] == 1){ //If the arc exists
										expr_1stLevel += x[v][j1][i][t-inst.travelTime[v][j1-1][i-1]];
										expr_1stFlow += fX[v][j1][i][t-inst.travelTime[v][j1-1][i-1]];
                                        #ifndef NBranching
                                        //~ expr_sumEnteringX += x[v][j1][i][t-inst.travelTime[v][j1-1][i-1]];
                                        #endif
									}
								}
							}
							if(j1>0){ //No consider source arc (second level balance)
								if (hasArc[v][i][j1][t] == 1){
									if(j1 == N){ //If j1 is the sink node, add to expr
										expr_2ndLevel += x[v][i][j1][t];
										expr_2ndFlow += fX[v][i][j1][t];
										if(addConstr){
											if (inst.q_v[v] <= inst.f_max_jt[i-1][0])
												expr_opd += x[v][i][j1][t];
										}
									}else if (t + inst.travelTime[v][i-1][j1-1] <= tF_add){ //If j1 is a port, it is necessary that the arrival is in the model
										expr_2ndLevel += x[v][i][j1][t];
										expr_2ndFlow += fX[v][i][j1][t];
										if(addConstr){
											if(inst.q_v[v] <= inst.f_max_jt[i-1][0] && inst.typePort[i-1] != inst.typePort[j1-1])
												expr_opd += x[v][i][j1][t];
										}
									}
								}
							}
						}
						if(addConstr){
							ss11 << "operate_and_depart_" << v << "(" << i << "," << t << ")";
							if(hasEnteringArc1st[v][i][t]==1 && inst.q_v[v] <= inst.f_max_jt[i-1][0]){
								#ifdef WaitAfterOperate                        
									expr_opd += -z[v][i][t];
								#endif
								#ifndef WaitAfterOperate                        
									expr_opd += -oA[v][i][t];
								#endif
								operateAndDepart[v][i][t] = IloRange(env, 0, expr_opd, 0,ss11.str().c_str());                    
								model.add(operateAndDepart[v][i][t]);
							}
						}
						IloExpr expr_link;
						#ifndef WaitAfterOperate 
						if (t < tF_add && hasEnteringArc1st[v][i][t-1]==0){ //Not last time period nor entering arc in the previous time period
							expr_1stLevel += - w[v][i][t] - oA[v][i][t];
							expr_2ndLevel += -oA[v][i][t] + oB[v][i][t];
							expr_1stFlow+= -fW[v][i][t] - fOA[v][i][t];
							expr_link = oA[v][i][t] - z[v][i][t];
							expr_2ndFlow += -fOA[v][i][t] - inst.delta[i-1]*f[v][i][t] + fOB[v][i][t];
						}else if (t==tF_add){ //Last time period
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
						if (t < tF_add && hasEnteringArc1st[v][i][t-1]==0){ //Not last time period nor entering arc in the previous time period
							expr_1stLevel += - w[v][i][t] - z[v][i][t];
							expr_2ndLevel += -z[v][i][t] + wB[v][i][t];	
							expr_1stFlow+= -fW[v][i][t] - fOA[v][i][t];	
							expr_2ndFlow += -fOA[v][i][t] - inst.delta[i-1]*f[v][i][t] + fWB[v][i][t];
						}else if (t==tF_add){ //Last time period
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
						if(hasEnteringArc1st[v][i][t]==1){
							model.add(firstLevelBalance[v][i][t]);
						
							secondLevelBalance[v][i][t].setName(ss1.str().c_str());
							model.add(secondLevelBalance[v][i][t]);						
						}						
						
						#ifndef WaitAfterOperate
						linkBalance[v][i][t].setName(ss2.str().c_str());
						model.add(linkBalance[v][i][t]);
						#endif
						
						if(hasEnteringArc1st[v][i][t]==1){
							
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
							stringstream ss7_1;
							ss7_1 << "flowMinLimitOA_"<<v<<","<<i<<","<<t;
							#ifndef WaitAfterOperate
							if(tightenFlow){
								if(inst.typePort[i-1] == 1){ //Minimum flow in discharging ports
									flowMinCapacityOA[v][i][t] = IloRange(env, -IloInfinity, oA[v][i][t]*inst.f_min_jt[i-1][t-1]-fOA[v][i][t], 0, ss7_1.str().c_str());
									model.add(flowMinCapacityOA[v][i][t]);
								}else{ //Change upper bound for loading ports
									flowCapacityOA.setExpr(fOA[v][i][t] - oA[v][i][t]*(inst.q_v[v]-inst.f_min_jt[i-1][t-1]));
								}
							}
							#endif
							#ifdef WaitAfterOperate
							flowCapacityOA[v][i][t] = IloRange(env, -IloInfinity, fOA[v][i][t]-inst.q_v[v]*z[v][i][t], 0, ss7.str().c_str());
							if(tightenFlow){
								if(inst.typePort[i-1] == 1){ //Minimum flow in discharging port
									flowMinCapacityOA[v][i][t] = IloRange(env, -IloInfinity, z[v][i][t]*inst.f_min_jt[i-1][t-1]-fOA[v][i][t], 0, ss7_1.str().c_str());
									model.add(flowMinCapacityOA[v][i][t]);
								}else{ //Change upper bound for loading ports
									flowCapacityOA[v][i][t].setExpr(fOA[v][i][t] - z[v][i][t]*(inst.q_v[v]-inst.f_min_jt[i-1][t-1]));
								}
							}   
							#endif
							model.add(flowCapacityOA[v][i][t]);
							
							if(t<tF_add){ //Constraints with no last time index
								#ifndef WaitAfterOperate
								ss8 << "flowLimitOB_"<<v<<","<<i<<","<<t;
								stringstream ss8_1;
								ss8_1 << "flowMinLimitOB_"<<v<<","<<i<<","<<t;
								flowCapacityOB[v][i][t] = IloRange(env, -IloInfinity, fOB[v][i][t]-inst.q_v[v]*oB[v][i][t], 0, ss8.str().c_str());
								model.add(flowCapacityOB[v][i][t]);
								if(tightenFlow){
									if(inst.typePort[i-1] == 1){
										flowMinCapacityOB[v][i][t] = IloRange(env, -IloInfinity, inst.f_min_jt[i-1][t-1]*oB[v][i][t] - fOB[v][i][t], 0, ss8_1.str().c_str());
										model.add(flowMinCapacityOB[v][i][t]);
									}
								}
								#endif
								ss9 << "flowLimitW_"<<v<<","<<i<<","<<t;
								stringstream ss9_1;
								ss9_1 << "flowMinLimitW_"<<v<<","<<i<<","<<t;
								flowCapacityW[v][i][t] = IloRange(env, -IloInfinity, fW[v][i][t]-inst.q_v[v]*w[v][i][t], 0, ss9.str().c_str());
								model.add(flowCapacityW[v][i][t]);
								if(tightenFlow){
									if(inst.typePort[i-1] == 1){ //Lower bound on flow for discharging ports
										flowMinCapacityW[v][i][t] = IloRange(env, -IloInfinity, w[v][i][t]*inst.f_min_jt[i-1][t-1]-fW[v][i][t], 0, ss9_1.str().c_str());
										model.add(flowMinCapacityW[v][i][t]);
									}else{ //Upper bound on loading ports
										flowCapacityW[v][i][t].setExpr(fW[v][i][t]-w[v][i][t]*(inst.q_v[v]-inst.f_min_jt[i-1][t-1]));
									}
								}
								#ifdef WaitAfterOperate
								ss10 << "flowLimitWB_"<<v<<","<<i<<","<<t;
								stringstream ss10_1;
								ss10_1 << "flowMinLimitWB_"<<v<<","<<i<<","<<t;
								flowCapacityWB[v][i][t] = IloRange(env, -IloInfinity, fWB[v][i][t]-inst.q_v[v]*wB[v][i][t], 0, ss10.str().c_str());
								model.add(flowCapacityWB[v][i][t]);
								if(tightenFlow){
									//Equal for discharging and loading ports
									flowMinCapacityWB[v][i][t] = IloRange(env, -IloInfinity, wB[v][i][t]*inst.f_min_jt[i-1][t-1]-fWB[v][i][t], 0, ss10_1.str().c_str());
									model.add(flowMinCapacityWB[v][i][t]);
									flowCapacityWB[v][i][t].setExpr(wB[v][i][t]*(inst.q_v[v]-inst.f_min_jt[i-1][t-1])-fWB[v][i][t]);							
								}
								#endif
							}
						}
                        ///Another port-time iteraror (i,t)
                        if(v==0){
                            #ifndef NWWCCReformulation //TODO - Inventory balance and upper and lower bound are note necessary
                            stringstream ss12;
                            //Common for both loading and discharging ports
                            for(int k=1;k<=t;k++){
                                ss12.str(string());
                                ss12 << it_kt << "_wwcc_relaxation_" << i << "(" << k << "," << t << ")";                    
                                expr_wwcc.clear();
                                int wwcc_rhs=0;
                                bool hasExpr=false;                                
                                for(int u=k;u<=t;u++){
                                    for (int v1=0;v1<V;v1++){
                                        if(hasEnteringArc1st[v1][i][u]==1){
                                            expr_wwcc += z[v1][i][u];
                                            hasExpr = true;
                                        }
                                    }
                                    if (inst.typePort[i-1] == 1) //Discharging port
                                        wwcc_rhs += inst.dM_jt[i-1][u];
                                    else    //Loading port
                                        wwcc_rhs += inst.d_jt[i-1][u-1];
                                }
                                expr_wwcc *= inst.f_max_jt[i-1][t-1];                                
                                if(k==1){ //When k=1 and k-1=0, net inventory variable is 0
                                    if(inst.typePort[i-1] == 0){ //Loading port
                                        wwcc_rhs += inst.s_j0[i-1] - inst.sMax_jt[i-1][1];
                                    }                                    
                                }else{ 
                                    if(inst.typePort[i-1] == 1){ //Discharging port
                                        expr_wwcc += sP[i][k-1];
                                        wwcc_rhs += inst.sMinM_jt[i-1][k-1];
                                    }else{  //Loading port
                                        expr_wwcc -= sP[i][k-1];
                                        wwcc_rhs -= inst.sMax_jt[i-1][1]; // == -\overline{S}_i
                                    }                                    
                                }                
                                if(hasExpr || k>1){                                    
                                    //~ //~ //~ cout " " <<  expr_wwcc << " >= " << wwcc_rhs << endl;
                                    wwcc_relaxation[i][it_kt] = IloRange(env, wwcc_rhs, expr_wwcc, IloInfinity, ss12.str().c_str());                    
                                    model.add(wwcc_relaxation[i][it_kt]); 
                                }
                                it_kt++;
                            }
                            #endif
                        }
                        #ifndef NBranching
                        expr_sumOA += oA[v][i][t];
                        #endif
					}
				}
				#ifndef NAlpha
				if(v==0){
					cumSlack[i].setExpr(expr_cumSlack);
					if(proportionalAlpha){
						cumSlack[i].setUB(inst.alp_max_j[i-1]/T*tF_add);
					}
				}
				#endif
                #ifndef NBranching
                //~ priorityX[v][i].setExpr(expr_sumEnteringX);
                priorityOA[v][i].setExpr(expr_sumOA);
                #endif
			}
		}
		if(tS_add <= T)
			sinkNodeBalance[v].setExpr(expr_sinkFlow);
	}
	//New objective
	if(tS_add <= T)
		obj.setExpr(expr_obj_current + expr_obj_cost - expr_obj_revenue);
	///ending removing endblock

}
void Model::improvementPhase_timeIntervals(IloEnv& env, Instance inst, const double& mIntervals, 
    const double& timePerIter, const double& gap, const double& overlap, Timer<std::chrono::milliseconds>& timer_cplex,
    float& opt_time,const double& timeLimit, float& elapsed_time, double& incumbent){
    double prevObj = 1.0e+10;
	double objValue = incumbent;
	double prevObj1 = incumbent;
	int i;
    float timeCoef = 0.3; //Percentage in which the time is increased or decreased for interval
	Timer<chrono::milliseconds> timer_LS;
	timer_LS.start();
	cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
	if (gap > 1e-04)
		cplex.setParam(IloCplex::EpGap, gap/100);
	int tS, tF;
    int totalIntervals = ceil(mIntervals*(1+overlap/100));
    unsigned int itCount=1; //Count the macro iterarions and is used for feed the random seed 
    vector <unsigned int> interVector; //Used for selecting intervals at random without repetitions
    vector <pair <unsigned int,float> > intervalCoef;  // Controls the time coefficient according the evaluation of each interval <interval, coef>
    
    //Feeding the vectors of intervals and time coefs.
    for(i=1;i<=totalIntervals;i++){
        intervalCoef.push_back(make_pair(i,1.0));
    }
    
	while((prevObj - objValue >= 0.1) && elapsed_time/1000 < timeLimit){
		//Feed the vector of intervals
		#ifndef NRandomTimeInterval			
        for(i=1;i<=totalIntervals;i++){
            interVector.push_back(i);
        }
        #endif
    
        srand(++itCount);
        for(i=1;i<=totalIntervals;i++){            			
			int rIntervalId;
			#ifndef NRandomTimeInterval			
            rIntervalId = iRand(0,interVector.size()-1); // index on intervals vector (0,m)
            int intervalNumber = interVector[rIntervalId]; // number of the interval (1,m+1)
            interVector.erase(interVector.begin()+rIntervalId); //remove the current interval
            #endif
            #ifdef NRandomTimeInterval
            rIntervalId = i;
            int intervalNumber = rIntervalId;
            #endif
    
            if(intervalNumber==1){
				tS = 1;
			}else{ //Considering the overlap for itrations > 1
				tS = inst.t/mIntervals*(intervalNumber-1)*(1-overlap/100);
			}
			tF = min(tS+inst.t/mIntervals, (double)inst.t);
			//~ cout "Unfixing interval " << intervalNumber << " = " << tS << "..." << tF << endl;
			unFixInterval(inst, tS, tF);
	            //~ //~ cout "Solving...\n";
			cplex.setParam(IloCplex::TiLim, min(max(0.0,timePerIter*intervalCoef[intervalNumber-1].second), max(timeLimit-elapsed_time/1000,0.0)));
			timer_cplex.start();
			cplex.solve();
			opt_time += timer_cplex.total();
			elapsed_time += timer_LS.total();
			incumbent = cplex.getObjValue();
            
            //Increase or the decrease the time for the next iteration if improved the solution
            if(prevObj1-incumbent >= 0.1){
                intervalCoef[intervalNumber-1].second = intervalCoef[intervalNumber-1].second*(1+timeCoef);
            }else{
                intervalCoef[intervalNumber-1].second = intervalCoef[intervalNumber-1].second*(1-timeCoef);
			}
			prevObj1 = incumbent;
			
            if(elapsed_time/1000 >= timeLimit){
				//~ //~ cout "Reached LS time limit " << timeLimit << ": " << elapsed_time/1000 << endl;
				break;
			}
			timer_LS.start();
			//~ cout "Objective: " << incumbent << endl;
            //~ //~ cout "Re-fixing " << tS << "..." << tF;
            fixSolution(env, inst, tS, tF,0,true);
            //~ //~ cout ". Done! " << endl;

		}
		prevObj = objValue;
		objValue = incumbent;
		if(elapsed_time/1000 >= timeLimit){			
			break;
		}
	}
}
/*
 * VND strategy: if solving an interval improved the solution, re-start the search from the first interval
 */
void Model::improvementPhaseVND_timeIntervals(IloEnv& env, Instance inst, const double& mIntervals, 
    const double& timePerIter, const double& gap, const double& overlap, Timer<std::chrono::milliseconds>& timer_cplex,
    float& opt_time,const double& timeLimit, float& elapsed_time, double& incumbent,  
    unsigned int& stopsByGap,unsigned int& stopsByTime){
    double prevObj = incumbent;	
	int i=1;
    
	Timer<chrono::milliseconds> timer_LS;
	timer_LS.start();	
	if (gap > 1e-04)
		cplex.setParam(IloCplex::EpGap, gap/100);
	int tS, tF;
    int totalIntervals = ceil(mIntervals*(1+overlap/100));

	while(i <= totalIntervals){
	   if(i==1){
			tS = 1;				
			tF = ceil(inst.t/mIntervals);
		}else{ //Considering the overlap for itrations > 1
			tS = ceil(inst.t/mIntervals*(i-1)*(1-overlap/100)+1);
			tF = min(tS+inst.t/mIntervals-1, (double)inst.t);
		}		
		//~ cout << "Unfixing interval " << i << " = " << tS << "..." << tF << endl;
		unFixInterval(inst, tS, tF);
		
		optimize:
		cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
		timer_cplex.start();
		cplex.solve();
		if(cplex.getMIPRelativeGap() > 1e-04){
			++stopsByTime;
		}else{
			++stopsByGap;
		}		
		opt_time += timer_cplex.total();
		elapsed_time += timer_LS.total();
		timer_LS.start();
		incumbent = cplex.getObjValue();		
		if(prevObj - incumbent >= 0.1){ //Improved the solution
			//~ cout "Improved " << prevObj << " -> " << incumbent << endl;
			prevObj = incumbent;			
			if(i!=1){ //If is not the first iteration need to re-fix the interval				
				//~ cout "Re-fixing " << tS << "..." << tF << endl;
				fixSolution(env, inst, tS, tF,0,true);
				i=1;
			}else{ //Continue optimizing first interval
				//~ cout "Solve again \n";
				goto optimize;
			}            
			
		}else{ //Not improved solution
			i++;
			//~ cout "Not improved Re-fixing " << tS << "..." << tF << endl;
			fixSolution(env, inst, tS, tF,0,true);
		}
		if(elapsed_time/1000 >= timeLimit){			
			break;
		}
	}
}

void Model::unFixInterval(Instance inst, const int& tS, const int& tF){
    int i,v,j,t;
    int T=inst.t;
    int V=inst.speed.getSize();
    int J=inst.numTotalPorts;
    int N=J+1;
    for(v=0;v<V;v++){
        for(i=1;i<=J;i++){
            for(j=1;j<=N;j++){
                for(t=tS;t<=tF;t++){
                    if(hasArc[v][i][j][t] == 1){
						if(i != j){
							if(j < N){ //If j is a port
								int t2 = t + inst.travelTime[v][i-1][j-1];
								if (t2<=tF)
									x[v][i][j][t].setBounds(0,1);
							}else{ //If j is the sink node
								x[v][i][j][t].setBounds(0,1);
							}
						}
					}
                    if(j==1){ //(v,i,t) iterator
                        if(hasEnteringArc1st[v][i][t]==1){
							z[v][i][t].setBounds(0,1);
							#ifndef WaitAfterOperate
							oA[v][i][t].setBounds(0,1);
							#endif
							if(t<tF){
								w[v][i][t].setBounds(0,1);
								#ifdef WaitAfterOperate
								wB[v][i][t].setBounds(0,1);
								#endif
								#ifndef WaitAfterOperate
								oB[v][i][t].setBounds(0,1);
								#endif
							}
						}
                    }
                }
            }
        }
    }
}

void Model::improvementPhase_typePortsLS(IloEnv env, Instance inst, const double& timePerIter, const int& gap, Timer<chrono::milliseconds>& timer_cplex,float& opt_time,
const double& timeLimit, float& elapsed_time, double& incumbent, unsigned int& stopsByGap,unsigned int& stopsByTime){
	Timer<chrono::milliseconds> timer_LS;
	timer_LS.start();
	int v,t,a, idPort;
	int i,r,j,j1;
    int T = inst.t;
    int J = inst.numTotalPorts;
    int N = J+1;
    int V = inst.speed.getSize();
	double objValue = incumbent;
	double prevObj = 1.0e+10;
	
	cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
	if (gap > 1e-04)
		cplex.setParam(IloCplex::EpGap, gap/100);
	
	while((prevObj - objValue >= 0.1) && elapsed_time/1000 < timeLimit){
		prevObj = objValue;
		for(i=1;i>=0;i--){ //Type region (first allow discharging region)
			//~ cout << "Unfixing TYPE REGION = " << i << endl;
			///Allow region variables
            for(r=0;r<inst.identifyPort[i].getSize();r++){ //For each region
				for(j=0;j<inst.identifyPort[i][r].getSize();j++){ //Port
					idPort = inst.identifyPort[i][r][j]+1;                    
					//Allow vessels/port-time variables
					for(v=0;v<V;v++){
						for(t=1;t<=T;t++){
							if(hasEnteringArc1st[v][idPort][t]==1){
								z[v][idPort][t].setBounds(0,1);
								if(t<T){
									w[v][idPort][t].setBounds(0,1);
									 #ifdef WaitAfterOperate
									wB[v][idPort][t].setBounds(0,1);
									 #endif
									#ifndef WaitAfterOperate
									oB[v][idPort][t].setBounds(0,1);                                    
									#endif
								}
								#ifndef WaitAfterOperate
								oA[v][idPort][t].setBounds(0,1);                                
								#endif
							}
						}
					}
					//Allow X variables to be solved
					for(v=0;v<inst.speed.getSize();v++){
						for(j1=1;j1<=N;j1++){
							if(idPort != j1){
                                for(t=1;t<=T;t++){
									if(hasArc[v][idPort][j1][t]==1){
										x[v][idPort][j1][t].setBounds(0,1);
									}
                                }
                            }
						}
					}
				}
			}
			//Solve
			cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
			timer_cplex.start();
			cplex.solve();
			if(cplex.getMIPRelativeGap() > 1e-04){
				++stopsByTime;
			}else{
				++stopsByGap;
			}
			opt_time += timer_cplex.total();
			elapsed_time += timer_LS.total();
			incumbent = cplex.getObjValue();
			//~ //~ cout "Objective: " << objValue << endl << endl;
			if(elapsed_time/1000 >= timeLimit){				
                break;
			}
			timer_LS.start();
			/////TODO - Non-optimized version. 
			///Get solution values
			for(r=0;r<inst.identifyPort[i].getSize();r++){ //For each region
				for(j=0;j<inst.identifyPort[i][r].getSize();j++){ //Port
					idPort = inst.identifyPort[i][r][j]+1;
					//Get vessels/port-time variables values
					for(v=0;v<inst.speed.getSize();v++){
						for(t=1;t<=T;t++){
							if(hasEnteringArc1st[v][idPort][t]==1){
								zValue[v][idPort][t] = cplex.getValue(z[v][idPort][t]);
								if(t<T){
									wValue[v][idPort][t] = cplex.getValue(w[v][idPort][t]);
									#ifdef WaitAfterOperate
									wBValue[v][idPort][t] = cplex.getValue(wB[v][idPort][t]);
									#endif
									#ifndef WaitAfterOperate
									oBValue[v][idPort][t] = cplex.getValue(oB[v][idPort][t]);
									#endif
								}
								#ifndef WaitAfterOperate
								oAValue[v][idPort][t] = cplex.getValue(oA[v][idPort][t]);
								#endif
							}
						}
					}
					
					//Get x variables value
					for(v=0;v<inst.speed.getSize();v++){
						for(j1=1;j1<=J;j1++){
							if(idPort != j1){
                                for(t=1;t<=T;t++){
									if(hasArc[v][idPort][j1][t] == 1){
										xValue[v][idPort][j1][t] = cplex.getValue(x[v][idPort][j1][t]);
									}
                                }
                            }
						}
					}
				}
			}
            //~ cout << "Fixing TYPE REGION = " << i << endl;
            ///Re-fix the solved variables
			for(r=0;r<inst.identifyPort[i].getSize();r++){ //For each region
				for(j=0;j<inst.identifyPort[i][r].getSize();j++){ //Port
					idPort = inst.identifyPort[i][r][j]+1;
					//Set vessels/port-time variables values
					for(v=0;v<inst.speed.getSize();v++){
						for(t=1;t<=T;t++){
							if(hasEnteringArc1st[v][idPort][t]==1){
								z[v][idPort][t].setBounds(round(zValue[v][idPort][t]),round(zValue[v][idPort][t]));
								if(t<T){
									w[v][idPort][t].setBounds(round(wValue[v][idPort][t]),round(wValue[v][idPort][t]));
									#ifdef WaitAfterOperate
									wB[v][idPort][t].setBounds(round(wBValue[v][idPort][t]),round(wBValue[v][idPort][t]));
									#endif
									#ifndef WaitAfterOperate
									oB[v][idPort][t].setBounds(round(oBValue[v][idPort][t]),round(oBValue[v][idPort][t]));
									#endif
								}
								#ifndef WaitAfterOperate
								oA[v][idPort][t].setBounds(round(oAValue[v][idPort][t]),round(oAValue[v][idPort][t]));
								#endif
							}
						}
					}
					
					//Set x variables value
					for(v=0;v<inst.speed.getSize();v++){
						for(j1=1;j1<=J;j1++){
							if(idPort != j1){
                                for(t=1;t<=T;t++){
									if(hasArc[v][idPort][j1][t] == 1){
										x[v][idPort][j1][t].setBounds(round(xValue[v][idPort][j1][t]),round(xValue[v][idPort][j1][t]));
									}
                                }
                            }
						}
					}
				}
			}
		}
		objValue = incumbent;
	}

}

void Model::warmStart(IloEnv env, Instance inst, const double& timePerIter){
    int v,i,j,t;
    int T = inst.t;
    int J = inst.numTotalPorts;
    int N = J+1;
    int V = inst.speed.getSize();
    double objValue;
	//~ //~ cout "Warm starting..." << endl;    
    //Unfix all variables
    for(v=0;v<V;v++){
        for(i=1;i<=J;i++){
            for(j=1;j<=N;j++){
                for(t=1;t<=T;t++){
                    if(i != j){
                        x[v][i][j][t].setBounds(0,1);
                    }
                    if(j==1){ //(v,i,t) iterator
                        z[v][i][t].setBounds(0,1);
                        if(t<T){
                            w[v][i][t].setBounds(0,1);
                            #ifdef WaitAfterOperate
                            wB[v][i][t].setBounds(0,1);
                            #endif
                            #ifndef WaitAfterOperate
                            oB[v][i][t].setBounds(0,1);
                            #endif
                        }
                        #ifndef WaitAfterOperate
                        oA[v][i][t].setBounds(0,1);
                        #endif
                    }
                }
            }
        }
    }
        
    //Solve
    cplex.setParam(IloCplex::EpGap, 0.0001);
    cplex.setParam(IloCplex::TiLim, timePerIter);
    cplex.solve();    
    objValue = cplex.getObjValue();
    //~ //~ cout "Objective: " << objValue << endl;
}
void Model::fixAllSolution(IloEnv& env, const Instance& inst){
    int v,i,j,t;
    int T = inst.t;
    int J = inst.numTotalPorts;
    int N = J+1;
    int V = inst.speed.getSize();
    //Get the solution values
    getSolution(env,inst);
    for(v=0;v<V;v++){
        for(i=1;i<=J;i++){
            for(j=1;j<=N;j++){
                for(t=1;t<=T;t++){
                    if(i != j && hasArc[v][i][j][t]==1){
                        x[v][i][j][t].setBounds(round(xValue[v][i][j][t]),round(xValue[v][i][j][t]));
                    }
                    if(j==1){ //(v,i,t) iterator
                        if(hasEnteringArc1st[v][j][t] == 1) {
							z[v][i][t].setBounds(round(zValue[v][i][t]),round(zValue[v][i][t]));
							if(t<T){
								w[v][i][t].setBounds(round(wValue[v][i][t]),round(wValue[v][i][t]));
								#ifdef WaitAfterOperate
								wB[v][i][t].setBounds(round(wBValue[v][i][t]),round(wBValue[v][i][t]));
								#endif
								#ifndef WaitAfterOperate
								oB[v][i][t].setBounds(round(oBValue[v][i][t]),round(oBValue[v][i][t]));
								#endif
							}
							#ifndef WaitAfterOperate
							oA[v][i][t].setBounds(round(oAValue[v][i][t]),round(oAValue[v][i][t]));
							#endif
						}
                    }
                }
            }
        }
    }
}

void Model::unFixVessel(Instance inst, const int& v){
    int i,j,t;
    int T = inst.t;
    int J = inst.numTotalPorts;
    int N = J+1;
    for(i=1;i<=J;i++){
        for(j=1;j<=N;j++){
            for(t=1;t<=T;t++){
                if(i != j)
                    x[v][i][j][t].setBounds(0,1);
                
                if(j==1){ //(v,i,t) iterator
                    z[v][i][t].setBounds(0,1);
                    #ifndef WaitAfterOperate
                    oA[v][i][t].setBounds(0,1);
                    #endif
                    if(t<T){
                        w[v][i][t].setBounds(0,1);
                        #ifdef WaitAfterOperate
                        wB[v][i][t].setBounds(0,1);
                        #endif
                        #ifndef WaitAfterOperate
                        oB[v][i][t].setBounds(0,1);
                        #endif
                    }
                }
            }
        }
    }
}
void Model::fixVesselLessInterval(IloEnv env, Instance inst, const int& v, const int& tS, const int& tF){
    int i,j,t;
    int T = inst.t;
    int J = inst.numTotalPorts;
    int N = J+1;
    int V = inst.speed.getSize();
    getSolValsW(env, inst, 0, T, true); //TODO define a method that get the values of times out of interval [tS,tF]
    
    for(i=1;i<=J;i++){
        for(j=1;j<=N;j++){
            for(t=1;t<=T;t++){
                if(t < tS || t > tF){   //Only fix if t not belongs to the interval [tS,tF]
                    if(i != j){
                        if(hasArc[v][i][j][t]==1){
                            if(j<N){ //j is a port                            
                                int t2 = t + inst.travelTime[v][i-1][j-1];
                                if (t2 < tS || t2 > tF)
                                    x[v][i][j][t].setBounds(round(xValue[v][i][j][t]),round(xValue[v][i][j][t]));
                            }else //j is the sink node
                                x[v][i][j][t].setBounds(round(xValue[v][i][j][t]),round(xValue[v][i][j][t]));
                        }
                    }
                    if(j==1){ //(v,i,t) iterator
                        if(hasEnteringArc1st[v][j][t]==1){
                            z[v][i][t].setBounds(round(zValue[v][i][t]),round(zValue[v][i][t]));
                            if(t<T){ 
                                w[v][i][t].setBounds(round(wValue[v][i][t]),round(wValue[v][i][t]));
                                #ifdef WaitAfterOperate
                                wB[v][i][t].setBounds(round(wBValue[v][i][t]),round(wBValue[v][i][t]));
                                #endif
                                #ifndef WaitAfterOperate
                                oB[v][i][t].setBounds(round(oBValue[v][i][t]),round(oBValue[v][i][t]));
                                #endif
                            }
                            if(t-1 == tF && tF<T){ //Also fix these variables with coefficent tF
                                w[v][i][t-1].setBounds(round(wValue[v][i][t-1]),round(wValue[v][i][t-1]));
                                #ifdef WaitAfterOperate
                                wB[v][i][t-1].setBounds(round(wBValue[v][i][t-1]),round(wBValue[v][i][t-1]));
                                #endif
                                #ifndef WaitAfterOperate
                                oB[v][i][t-1].setBounds(round(oBValue[v][i][t-1]),round(oBValue[v][i][t-1]));
                                #endif
                            }
                            #ifndef WaitAfterOperate
                            oA[v][i][t].setBounds(round(oAValue[v][i][t]),round(oAValue[v][i][t]));
                            #endif
                        }
                    }
                }
            }
        }
    }
}

void Model::improvementPhaseVND_intervalVessel(IloEnv& env, Instance inst, const double& mIntervals, const double& timePerIter, 
    const double& gap, const double& overlap, Timer<chrono::milliseconds>& timer_cplex,float& opt_time, 
    const double& timeLimit, float& elapsed_time, double& incumbent, unsigned int& stopsByGap,unsigned int& stopsByTime){        		
	double prevObj = incumbent;
		
	int V = inst.speed.getSize();    
    int i,v, tS, tF, totalIntervals;    
    totalIntervals = ceil(mIntervals*(1+overlap/100));
    
	Timer<chrono::milliseconds> timer_LS;
	timer_LS.start();	
	if (gap > 1e-04)
		cplex.setParam(IloCplex::EpGap, gap/100);
        
    i=1;
    while(i<=totalIntervals){
		optmizeFromScratch:
		if(i==1){
			tS = 1;				
			tF = ceil(inst.t/mIntervals);
		}else{ //Considering the overlap for itrations > 1
			tS = ceil(inst.t/mIntervals*(i-1)*(1-overlap/100)+1);
			tF = min(tS+inst.t/mIntervals-1, (double)inst.t);
		}
		//~ cout << "Unfixing [" << tS << "..." << tF << "]" << endl;
        unFixInterval(inst, tS, tF);
		v=0;
		while(v<V){	
			//~ cout << "Unfixing vessel " << v << endl;
			unFixVessel(inst,v);			
			optmizeSameSubproblem:
			cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
			timer_cplex.start();
			cplex.solve();
			if(cplex.getMIPRelativeGap() > 1e-04){
				++stopsByTime;
			}else{
				++stopsByGap;
			}
			opt_time += timer_cplex.total();
			elapsed_time += timer_LS.total();
			timer_LS.start();			
			incumbent = cplex.getObjValue();
			if(prevObj-incumbent >= 0.1){ //Improved
				//~ cout << "Improved " << prevObj << " -> " << incumbent << endl;
				prevObj = incumbent;
				if(i==1 && v==0){ //Improvement on the first subproblem (does not need change the variables fixing/unfixing)
					//~ cout << "re-optimize i=1 v=0\n";
					goto optmizeSameSubproblem;
				}else if(i==1){
					fixVesselLessInterval(env,inst, v, tS, tF);
					//~ cout << "Fixing vessel " << v << endl;
					v=0;					
				}else{ //i>1 and v>0
					//~ cout << "Fixing vessel " << v << " and ";
					//~ cout << "Fixing " << tS << "..." << tF << endl;
					fixAllSolution(env,inst);
					v=0;
					i=1;
					goto optmizeFromScratch;				
				}
			}else{ //Not improved
				//~ cout << "No improvement\n";		
				fixVesselLessInterval(env,inst, v, tS, tF);
				//~ cout << "Fixing vessel " << v << endl;		
				v++;
			}		
			if(elapsed_time/1000 >= timeLimit){
				break;
			}
		}
		//~ cout << "Fixing " << tS << "..." << tF << endl;			
		fixSolution(env, inst, tS, tF,1,true);
		i++;
		if(elapsed_time/1000 >= timeLimit){
			break;	
		}			
	}
}
void Model::improvementPhase_intervalVessel(IloEnv& env, Instance inst, const double& mIntervals, const double& timePerIter, 
    const double& gap, const double& overlap, Timer<chrono::milliseconds>& timer_cplex,float& opt_time, 
    const double& timeLimit, float& elapsed_time, double& incumbent){
        
	//Stores the objective for each iteration (m*v subproblems)
	double objValue = incumbent;
	double prevObj = 1.0e+10;
	//Stores the objective of each subproblem
	double objValue1 = incumbent;
	double prevObj1 = 1.0e+10;
	
	int V = inst.speed.getSize();    
    int i,v, tS, tF;
    
    cplex.setParam(IloCplex::TiLim, timePerIter);
	Timer<chrono::milliseconds> timer_LS;
	timer_LS.start();
	
	if (gap > 1e-04)
		cplex.setParam(IloCplex::EpGap, gap/100);
        
    //For the random selection of vessels and variable processing time 
    float timeCoef = 0.3;    
    unsigned int sumIt = 1, rId, combId, vId;    
    vector <unsigned int> vesselsVector; // <v>  v = [0....V)
    vector <float> combCoefVector(ceil(mIntervals)*V); // <coef>
    fill(combCoefVector.begin(),combCoefVector.end(),1);
	
	while((prevObj - objValue >= 0.1) && elapsed_time/1000 <= timeLimit){
		//~ //~ cout "prevObj - objValue = " << prevObj - objValue << endl;
        srand(sumIt++);
        for(i=1;i<=ceil(mIntervals);i++){
			if(i==1)
				tS = 1;
			else
				tS = inst.t/mIntervals*(i-1)*(1-overlap/100);
			tF = min(tS+inst.t/mIntervals, (double)inst.t);
			//~ //~ cout "Unfixing interval " << tS << "..." << tF << endl;
			unFixInterval(inst, tS, tF);
            //~ for (v=0;v<V;v++){
                //~ vesselsVector.push_back(v);
            //~ }
			for (v=0;v<V;v++){
				///For random vessel selection
                //~ rId = iRand(0,vesselsVector.size()-1);
                //~ vId = vesselsVector[rId];
                //~ combId = (i-1)*V+vId;
                //~ //~ cout "Vessel id " << vId << " combId " << combId << endl;
                //~ vesselsVector.erase(vesselsVector.begin()+rId);
				//~ //~ cout "Unfixing vessel " << vId << " coef " << combCoefVector[combId] << endl;
				
				vId=v;
				combId = (i-1)*V+vId;
				unFixVessel(inst,vId);
				cplex.setParam(IloCplex::TiLim, min(max(0.0,timePerIter*combCoefVector[combId]), max(timeLimit-elapsed_time/1000,0.0)));
				timer_cplex.start();
				cplex.solve();
				opt_time += timer_cplex.total();
				elapsed_time += timer_LS.total();
				timer_LS.start();				
				incumbent = cplex.getObjValue();
				//~ //~ cout "Objective: " << objValue1 << endl;
                //Update the vector coef
                if(prevObj1-incumbent >= 0.1)
                    combCoefVector[combId] = combCoefVector[combId]*(1+timeCoef);
                else
                    combCoefVector[combId] = combCoefVector[combId]*(1-timeCoef);
				prevObj1 = incumbent;
				//~ //~ cout "Fixing vessel " << vId << endl;
				fixVesselLessInterval(env,inst, vId, tS, tF); //Get all solution values
				if(elapsed_time/1000 >= timeLimit){
					break;
                }
			}
            //~ //~ cout "Re-fixing " << tS << "..." << tF;
            fixSolution(env, inst, tS, tF, 1,true);
            //~ //~ cout ". Done! " << endl;

            if(elapsed_time/1000 >= timeLimit){
				break;
				//~ //~ cout "Stopped by time" << endl;
			}
		}
		prevObj = objValue;
		objValue = incumbent;
		if (elapsed_time/1000 >= timeLimit){
			//~ //~ cout "STOP by TIME \n";
			break;
		}
	}
}

void Model::fixVesselPair(IloEnv env, Instance inst, const int& v,const int& v1,const bool& getValues){
    int i,j,t;
    int T = inst.t;
    int J = inst.numTotalPorts;
    int N = J+1;
    //Get the values
    if(getValues){
		for(i=1;i<=J;i++){
			for(j=1;j<=N;j++){
				for(t=1;t<=T;t++){
					if(i != j){
						if(hasArc[v][i][j][t]==1)
							xValue[v][i][j][t] = cplex.getValue(x[v][i][j][t]);
						if(hasArc[v1][i][j][t]==1)
							xValue[v1][i][j][t] = cplex.getValue(x[v1][i][j][t]);
					}
					if(j==1){ //(v,i,t) iterator
						if(hasEnteringArc1st[v][i][t]==1){
							zValue[v][i][t] = cplex.getValue(z[v][i][t]);
							#ifndef WaitAfterOperate
							oAValue[v][i][t] = cplex.getValue(oA[v][i][t]);
							#endif
							if(t<T){
								wValue[v][i][t] = cplex.getValue(w[v][i][t]);
								#ifdef WaitAfterOperate
								wBValue[v][i][t] = cplex.getValue(wB[v][i][t]);
								#endif
								#ifndef WaitAfterOperate
								oBValue[v][i][t] = cplex.getValue(oB[v][i][t]);
								#endif
							}
						}
						if(hasEnteringArc1st[v1][i][t]==1){
							zValue[v1][i][t] = cplex.getValue(z[v1][i][t]);
							#ifndef WaitAfterOperate
							oAValue[v1][i][t] = cplex.getValue(oA[v1][i][t]);
							#endif
							if(t<T){                            
								wValue[v1][i][t] = cplex.getValue(w[v1][i][t]);
								#ifdef WaitAfterOperate                            
								wBValue[v1][i][t] = cplex.getValue(wB[v1][i][t]);
								#endif
								#ifndef WaitAfterOperate                            
								oBValue[v1][i][t] = cplex.getValue(oB[v1][i][t]);                        
								#endif
							}
						}
					}
				}
			}
		}
	}
    //Fix vessels
    for(i=1;i<=J;i++){
        for(j=1;j<=N;j++){
            for(t=1;t<=T;t++){
                if(i != j){
                    if(hasArc[v][i][j][t]==1)
                        x[v][i][j][t].setBounds(round(xValue[v][i][j][t]),round(xValue[v][i][j][t]));
                    if(hasArc[v1][i][j][t]==1)
                        x[v1][i][j][t].setBounds(round(xValue[v1][i][j][t]),round(xValue[v1][i][j][t]));
                }
                if(j==1){ //(v,i,t) iterator
                    //v
                    if(hasEnteringArc1st[v][i][t]==1){
                        z[v][i][t].setBounds(round(zValue[v][i][t]),round(zValue[v][i][t]));
                        #ifndef WaitAfterOperate
                        oA[v][i][t].setBounds(round(oAValue[v][i][t]),round(oAValue[v][i][t]));
                        #endif
                        if(t<T){
                            w[v][i][t].setBounds(round(wValue[v][i][t]),round(wValue[v][i][t]));
                            #ifdef WaitAfterOperate
                            wB[v][i][t].setBounds(round(wBValue[v][i][t]),round(wBValue[v][i][t]));
                            #endif
                            #ifndef WaitAfterOperate
                            oB[v][i][t].setBounds(round(oBValue[v][i][t]),round(oBValue[v][i][t]));
                            #endif
                        }
                    }
                    //v1
                    if(hasEnteringArc1st[v1][i][t]==1){
                        z[v1][i][t].setBounds(round(zValue[v1][i][t]),round(zValue[v1][i][t]));
                        #ifndef WaitAfterOperate
                        oA[v1][i][t].setBounds(round(oAValue[v1][i][t]),round(oAValue[v1][i][t]));
                        #endif
                        if(t<T){
                            w[v1][i][t].setBounds(round(wValue[v1][i][t]),round(wValue[v1][i][t]));
                            #ifdef WaitAfterOperate
                            wB[v1][i][t].setBounds(round(wBValue[v1][i][t]),round(wBValue[v1][i][t]));
                            #endif
                            #ifndef WaitAfterOperate
                            oB[v1][i][t].setBounds(round(oBValue[v1][i][t]),round(oBValue[v1][i][t]));
                            #endif
                        }
                    }
                }
            }
        }
    }
}
void Model::fixVessel(IloEnv env, Instance inst, const int& v,const bool& getValues){
    int i,j,t;
    int T = inst.t;
    int J = inst.numTotalPorts;
    int N = J+1;
    //Get the values
    if(getValues){
		for(i=1;i<=J;i++){
			for(j=1;j<=N;j++){
				for(t=1;t<=T;t++){
					if(i != j){
						if(hasArc[v][i][j][t]==1){
							xValue[v][i][j][t] = cplex.getValue(x[v][i][j][t]);
						}						
					}
					if(j==1){ //(v,i,t) iterator
						if(hasEnteringArc1st[v][i][t]==1){
							zValue[v][i][t] = cplex.getValue(z[v][i][t]);
							#ifndef WaitAfterOperate
							oAValue[v][i][t] = cplex.getValue(oA[v][i][t]);
							#endif
							if(t<T){
								wValue[v][i][t] = cplex.getValue(w[v][i][t]);
								#ifdef WaitAfterOperate
								wBValue[v][i][t] = cplex.getValue(wB[v][i][t]);
								#endif
								#ifndef WaitAfterOperate
								oBValue[v][i][t] = cplex.getValue(oB[v][i][t]);
								#endif
							}
						}						
					}
				}
			}
		}
	}
    //Fix vessels
    for(i=1;i<=J;i++){
        for(j=1;j<=N;j++){
            for(t=1;t<=T;t++){
                if(i != j){
                    if(hasArc[v][i][j][t]==1){
                        x[v][i][j][t].setBounds(round(xValue[v][i][j][t]),round(xValue[v][i][j][t]));
					}
                }
                if(j==1){ //(v,i,t) iterator
                    //v
                    if(hasEnteringArc1st[v][i][t]==1){
                        z[v][i][t].setBounds(round(zValue[v][i][t]),round(zValue[v][i][t]));
                        #ifndef WaitAfterOperate
                        oA[v][i][t].setBounds(round(oAValue[v][i][t]),round(oAValue[v][i][t]));
                        #endif
                        if(t<T){
                            w[v][i][t].setBounds(round(wValue[v][i][t]),round(wValue[v][i][t]));
                            #ifdef WaitAfterOperate
                            wB[v][i][t].setBounds(round(wBValue[v][i][t]),round(wBValue[v][i][t]));
                            #endif
                            #ifndef WaitAfterOperate
                            oB[v][i][t].setBounds(round(oBValue[v][i][t]),round(oBValue[v][i][t]));
                            #endif
                        }
                    }                    
                }
            }
        }
    }
}

void Model::improvementPhase_vessels(IloEnv& env, Instance inst, const double& timePerIter, const double& gap, double& incumbent, Timer<chrono::milliseconds>& timer_cplex,float& opt_time,
const double& timeLimit, float& elapsed_time){
	Timer<chrono::milliseconds> timer_LS;
	timer_LS.start();
	int v, v1, V = inst.speed.getSize();
	if (gap > 1e-04)
		cplex.setParam(IloCplex::EpGap, gap/100);
	cplex.setParam(IloCplex::TiLim, min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
	double currentObj = incumbent;
	double prevObj = 1e+20;
	double prevObj1 = 1e+20;
	float elapsed_local_time{0};
    float timeCoef = 0.3;
    unsigned int combId;
    unsigned int itCount=1; //Count the macro iteration and feed the random seed

	//Get random pairs of vessels
	int count, rId, v2, sumC =0, sumComb = 0;  
	vector <tuple<unsigned int,unsigned int, unsigned int> > vesselsComb;    // <id, v1, v2>
    vector <pair<unsigned int, float> > vesselsCombTimeCoef;                 //<id, timeCoef>
	srand(V);
    //Building vector of timeCoeficients
    sumC = 0;
    for(v=0;v<V-1;v++){
        for(v1=v+1;v1<V;v1++){
            vesselsCombTimeCoef.push_back(make_pair(sumC,1.0));
            sumC++;
        }
    }
	while ((prevObj - currentObj >= 0.1) && 
        elapsed_local_time/1000 < timeLimit){		
        //Build the combinations
		sumComb = 0;
		for(v=0;v<V-1;v++){
			for(v1=v+1;v1<V;v1++){
				vesselsComb.push_back(tuple<unsigned int,unsigned int,unsigned int> (sumComb,v,v1));
				sumComb++;
			}
		}
		for (int i=0;i<sumComb;i++){
			//~ //~ cout "Size: " << vesselsComb.size() << ": ";
			rId = iRand(0,vesselsComb.size()-1);
			combId = get<0>(vesselsComb[rId]); //Id of the combination -> used for iterate in the vesselsCombTimeCoef vector
            v1 = get<1>(vesselsComb[rId]);
			v2 = get<2>(vesselsComb[rId]);
			vesselsComb.erase(vesselsComb.begin()+rId);
           
			//~ //~ cout "Improving vessels " << v1 << " and " << v2 << " Coef " << vesselsCombTimeCoef[combId].second << endl;
			unFixVessel(inst,v1);
			unFixVessel(inst,v2);
			
			cplex.setParam(IloCplex::TiLim, min(max(0.0,timePerIter*vesselsCombTimeCoef[combId].second), max(timeLimit-elapsed_local_time/1000,0.0)));
			timer_cplex.start();
			cplex.solve();

			vesselsCombTimeCoef[combId].second = vesselsCombTimeCoef[combId].second*(1+timeCoef); //Increase the time timeCoeficient
		
			opt_time += timer_cplex.total();
			incumbent = cplex.getObjValue();
			if(prevObj1 - incumbent >= 0.1){			
				vesselsCombTimeCoef[combId].second = vesselsCombTimeCoef[combId].second*(1+timeCoef); //Increase the time timeCoeficient
			}else{			
				vesselsCombTimeCoef[combId].second = vesselsCombTimeCoef[combId].second*(1-timeCoef); //Decrease the time timeCoeficient
			}		
			prevObj1 = incumbent;	
			//~ //~ cout "Objective " << incumbent << endl;
			fixVesselPair(env, inst, v1,v2,true); //Fix the 2 vessels
			elapsed_local_time += timer_LS.total();
			if(elapsed_local_time/1000 >= timeLimit){
				break;
			}
			timer_LS.start();
		}
        //~ //~ cout endl;
		prevObj = currentObj;
		currentObj = incumbent;
	}
	elapsed_time += elapsed_local_time;
}
/*
 * Select itevatively a pair of vessel to be optimized.
 * If a subproblem improves the solution, the search is re-started from the first vessel pair
 */
void Model::improvementPhaseVND_vessels(IloEnv& env, Instance inst, const double& timePerIter, const double& gap, 
double& incumbent, Timer<chrono::milliseconds>& timer_cplex,float& opt_time, const double& timeLimit, float& elapsed_time,
	unsigned int& stopsByGap,unsigned int& stopsByTime){
	Timer<chrono::milliseconds> timer_LS;
	timer_LS.start();
	int v1, v2, V = inst.speed.getSize();
	if (gap > 1e-04)
		cplex.setParam(IloCplex::EpGap, gap/100);
		
	double prevObj = incumbent;    
    for(v1=0;v1<V;v1++){
		optimizeV1:
		//~ cout << "Unfixing vessel v1 " << v1 << endl;
		unFixVessel(inst,v1);
		for(v2=v1+1;v2<V;v2++){
			//~ cout << "unfixing vessel v2 " << v2 << endl;
			unFixVessel(inst,v2);
			
			optimizeV2:	
			if(elapsed_time/1000 >= timeLimit){				
				//~ cout << "Stop by time\n";
				break;					
			}		
			cplex.setParam(IloCplex::TiLim,min(timePerIter, max(timeLimit-elapsed_time/1000,0.0)));
			timer_cplex.start();
			cplex.solve();
			if(cplex.getMIPRelativeGap() > 1e-04){
				++stopsByTime;
			}else{
				++stopsByGap;
			}
			opt_time += timer_cplex.total();
			elapsed_time += timer_LS.total();
			timer_LS.start();			
			incumbent = cplex.getObjValue();
			getSolutionVesselPair(env,inst,v1,v2);
			if(prevObj-incumbent >= 0.1){ //Improved
				//~ cout << "Improved " << prevObj << " -> " << incumbent << endl;				
				prevObj = incumbent;					
				if(v1 == 0 && v2 == 1){ //Improvement on the first iteration
					//~ cout << "Re-optimize\n";
					goto optimizeV2;
				}else if(v1 == 0){ //When v1 = 0 Only fixes v2
					//~ cout << "Fixing v2 " << v2 << endl;
					fixVessel(env,inst,v2,false);
					v2=v1;						
				}else{ //Refix both vessels
					//~ cout << "Fixing " << v1 << " and " << v2 << endl;
					fixVesselPair(env, inst, v1,v2,false); 
					v1 = 0;
					v2 = 0;	
					goto optimizeV1;				
				}					
			}else{					
				//~ cout << "No improvement\n";	
				//~ cout << "Fixing v2 " << v2 << endl;
				fixVessel(env,inst,v2,false);
			}
		}		
		if(elapsed_time/1000 >= timeLimit){
			break;					
		}
		//~ cout << "Fixing v1 " << v1 << endl;
		fixVessel(env,inst,v1,false);		
	}
	//Use for verify if terminated by time (not change v1 and v2 index), or terminated normally (need change the index)
	if(v1==V) v1--;
	if(v2==V) v2--;
	//~ cout << "Fixing vessel v1 and v2 " << v1 << " " << v2 << endl;
	fixVesselPair(env,inst,v1,v2,false);
}
/*
 * Note: 07-11-19 - parameter nIntervals corresponds to the number of time periods per interval
 */  
void mirp::fixAndRelax(string file, string fixOptStr, const double& nIntervals, const double& gapFirst, const int& f, const double& overLap, const int& endBlock,
const int& timePerIterFirst, const double& mIntervals, const int& timePerIterSecond, const double& gapSecond,
 const double& overlap2, const double& timeLimit,  const bool& validIneq, const bool& addConstr, 
 const bool& tightenInvConstr, const bool& proportionalAlpha, const bool& reduceArcs, const bool& tightenFlow){
	///Time parameters
	Timer<chrono::milliseconds> timer_cplex;
	Timer<chrono::milliseconds> timer_global;
	Timer<chrono::milliseconds> timer_1stPhase;
	timer_global.start();
	float global_time {0};
	float opt_time {0};
	float elapsed_time {0};
	float elapsed_time_it {0};
	bool abortedRF = false;
	
	double obj1stPhase = 0, obj2ndPhase = 0, time1stPhase=0, time2ndPhase=0, intPart;
	
	///Read input files
	IloEnv env;
	
	try{
		Instance inst(env);
		inst.readInstance(env, file);
				
		int v, j, t, t1S, t1F, t2S, t2F, t3S, t3F;
		int J = inst.numTotalPorts;
		int T = inst.t;
		int V = inst.speed.getSize(); //# of vessels        

		/// NEW MODEL		
		Model model(env);
		model.buildFixAndRelaxModel(env,inst, nIntervals, endBlock, validIneq, addConstr, tightenInvConstr, 
			proportionalAlpha, reduceArcs, tightenFlow);
		model.setParameters(env, timePerIterFirst, gapFirst);
        //~ model.cplex.exportModel("mip_R&F.lp");
		//Relax-and-fix
		double p = nIntervals*(1-overLap/100); // Units of t that are add at each iteration to the model.
		int s = T-(nIntervals*endBlock); 		 // Last t (relaxed) of model when starting relax-and-fix.
		double sizeInterval = nIntervals;		 // Last t of integer block (always starting with 1 integer interval)		
        int maxIt = ceil((T-sizeInterval)/p);					
        
        //For write and reading a mst file
		size_t pos = file.find("LR");
		string str = file.substr(pos);
		str.erase(str.end()-1);
		stringstream ss;
		//ss << "test_no_beta/" << str << "_" << inst.t << ".mst";
		ss << str << "_" << inst.t << ".mst";
		//ss << "msts_no_form_adds/" << str << "_" << inst.t << ".mst";
		//ss << "msts_only_VI/" << str << "_" << inst.t << ".mst";
	    //ss << "msts_only_TightInventory/" << str << "_" << inst.t << ".mst";
	    //ss << "msts_only_AddConst/" << str << "_" << inst.t << ".mst";
		//ss << "msts_only_Simplification/" << str << "_" << inst.t << ".mst";
		//ss << "msts_only_PropAlpha/" << str << "_" << inst.t << ".mst";

		timer_1stPhase.start();
        if(mIntervals == 0){
			#ifdef NRelaxation
			for (v=1; v <= maxIt; v++){
				timer_cplex.start();
				//~ cout << "Iteration: " << v << "/" << maxIt << " - Solving..." << endl;
				if(!model.cplex.solve()){
					//~ cout << model.cplex.getStatus() << endl;
					opt_time += timer_cplex.total();
					time1stPhase += timer_1stPhase.total();
					abortedRF = true;
					obj1stPhase = 99999999;
					goto abortRF;
				}
				opt_time += timer_cplex.total();
				//~ cout "Solution Status " << model.cplex.getStatus() << " Value: " << model.cplex.getObjValue() << endl;
								
				t3S = ceil(sizeInterval * (v-1) * (1-overLap/100))+1;
				t3F = min(ceil(sizeInterval * v * (1-overLap/100)),(double)T);

				t2S = ceil(s+p*(v-1))+1;
				t2F = min(ceil(s+p*v), (double)T); 

				t1S = ceil(sizeInterval+p*(v-1)+1);
				t1F = min(ceil(sizeInterval+p*v),(double)T); 

				//~ //~ cout "Printing until time " << t2S-1 << endl;
				//~ model.printSolution(env, inst, t2S-1);

				model.modifyModel(env, inst, nIntervals, t3S, t3F, t2S, t2F, t1S, t1F, 
					validIneq, addConstr, tightenInvConstr, proportionalAlpha, tightenFlow);
				
				#ifndef FixedGAP
				double newGap = max(0.001, (gapFirst - gapFirst/maxIt*v)/100);
				//~ cout << "New GAP " << newGap*100 << " % \n \n";
				model.cplex.setParam(IloCplex::EpGap, newGap);            
				#endif
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
			abortRF:
			double incumbent = obj1stPhase;
			//~ //~ cout "Solution Status " << model.cplex.getStatus() << " Value: " << obj1stPhase << endl;
			
			//~ model.cplex.exportModel("mip_R&F.lp");
			if(!abortedRF){
				//~ cout << "Sucess!\n";
				//~ model.printSolution(env, inst, T);
				model.cplex.writeMIPStarts(ss.str().c_str());
			}
		}else{		
			//~ cout << "Reading MST file...\n";
			model.cplex.readMIPStarts(ss.str().c_str());
			model.cplex.setParam(IloCplex::TiLim, 1);        
			model.cplex.solve(); 
			obj1stPhase = model.cplex.getObjValue();       
        }
        double incumbent = obj1stPhase;
        unsigned int stopsByGap=0, stopsByTime=0;
        #ifndef NImprovementPhase	
		//~ //~ cout "\n\n\n\n IMPROVING SOLUTION... \n\n\n" << endl;
		model.cplex.setParam(IloCplex::EpGap, 1e-04);	//Set default GAP
        model.fixAllSolution(env, inst);
		double tLimit=timeLimit;
		for(string::iterator it=fixOptStr.begin();it!=fixOptStr.end();++it){
			tLimit = timeLimit;
			elapsed_time = 0;
			switch(*it){
				case 'a':
					//~ cout << *it << " - Interval Vessels " << tLimit << "\n";
					model.improvementPhaseVND_intervalVessel(env, inst, mIntervals, timePerIterSecond, gapSecond, overlap2, timer_cplex, opt_time, 
						tLimit, elapsed_time, incumbent, stopsByGap, stopsByTime);
					time2ndPhase += elapsed_time;
					break;
				case 'b':
					//~ cout << *it << " - Time Interval " << tLimit << "\n";
					model.improvementPhaseVND_timeIntervals(env, inst, mIntervals, timePerIterSecond, gapSecond, overlap2, timer_cplex, opt_time, 
						tLimit, elapsed_time, incumbent, stopsByGap, stopsByTime);
					time2ndPhase += elapsed_time;
					break;
				case 'c':
					//~ cout << *it << " - Vessel Pairs " << tLimit << "\n";
					model.improvementPhaseVND_vessels(env, inst, timePerIterSecond, gapSecond, incumbent, timer_cplex, opt_time, tLimit, elapsed_time, stopsByGap, stopsByTime);
					time2ndPhase += elapsed_time;
					break;
				case 'd':
					//~ cout << *it << " - Port Type " << tLimit << "\n";
					model.improvementPhase_typePortsLS(env, inst,timePerIterSecond, gapSecond, timer_cplex, opt_time, tLimit, elapsed_time, incumbent, stopsByGap, stopsByTime);
					time2ndPhase += elapsed_time;
					break;
				default:
					cerr << "No fix-and-optimize procedure for " << *it << endl;
					break;					
			}
		}		
						
        //~ model.improvementPhase_intervalVessel(env, inst, mIntervals, timePerIterSecond, gapSecond, overlap2, timer_cplex, opt_time, 
        //~ timeLimit/3, elapsed_time, incumbent);
                     
        //~ tLimit = (timeLimit - elapsed_time/1000);//2;        
        //~ //~ cout "Elapsed time: " << elapsed_time/1000 << " >> reaming: " << tLimit << endl;        
		
        //~ model.improvementPhase_timeIntervals(env, inst, mIntervals, timePerIterSecond, gapSecond, overlap2, timer_cplex, opt_time, 
        //~ tLimit, elapsed_time, incumbent); 
               
        //~ model.improvementPhase_vessels(env, inst, timePerIterSecond, gapSecond, incumbent, timer_cplex, opt_time, timeLimit, elapsed_time);
        		
        //~ model.warmStart(env,inst,timePerIterSecond*72);
        
        ///SOMENTE NESCESSARIO PARA OBTENÇÃO DE SOLUÇÃO COMPLETA
        model.cplex.setParam(IloCplex::TiLim, 1);        
        model.cplex.solve();
        ///
        #endif
        
		global_time += timer_global.total();
		obj2ndPhase	= incumbent;//model.cplex.getObjValue();
		
		//For getting information about solution
		//~ model.cplex.setParam(IloCplex::TiLim, 1000);
		//~ model.cplex.setParam(IloCplex::NodeLim, 1);
		//~ model.cplex.solve();
		#ifdef NRelaxation
		bool isInfeasible = false;
        if(!abortedRF){
			#ifndef NBetas
			double totalBeta=0;
			double totalTheta=0;
			IloArray<IloNumArray> betaVals(env, J+1); 
			IloArray<IloNumArray> thetaVals(env, J+1); 
			for(j=1;j<=J;j++){
				betaVals[j] = IloNumArray(env, T+1);
				thetaVals[j] = IloNumArray(env, T+1);
				for(t=1;t<=T;t++){
					betaVals[j][t] = model.cplex.getValue(model.beta[j][t]);
					totalBeta += betaVals[j][t];
					thetaVals[j][t] = model.cplex.getValue(model.theta[j][t]);
					totalTheta += thetaVals[j][t];
					if (betaVals[j][t] >= 0.01 || thetaVals[j][t] >= 0.01){
						isInfeasible = true;
					}
				}
			}
			#endif
		}
		
		//~ model.printSolution(env, inst, T);
		
		#ifndef NBetas
		//~ //~ cout "Total betas = " << totalBeta << endl;
		#endif
        #ifndef NThetas
		//~ //~ cout "Total thetas = " << totalTheta << endl;
		#endif
		
        /// Data (header in the main method)
        //~ cout << str << "\t" <<
				//~ maxIt << "\t" << 
				//~ time1stPhase/1000 << "\t" <<
				//~ obj1stPhase << "\t" << 
				//~ time2ndPhase/1000 << "\t" << 
				//~ obj2ndPhase << "\t" << 
				//~ opt_time/1000 << "\t" <<
				//~ (global_time-opt_time)/1000 << "\t" <<
				//~ abs((obj2ndPhase/obj1stPhase - 1)*100) << "\t" <<
				//~ isInfeasible << "\t" <<
				//~ stopsByGap << "\t" <<
				//~ stopsByTime << "\t" <<
				//~ endl;
		///For iRace tests: <obj,time>
		cout << obj1stPhase << " " << time1stPhase/1000 << endl;
		///Local search
		//~ cout << obj2ndPhase << " " << time2ndPhase/1000 << endl;
		
        #endif
        
        #ifndef NRelaxation        
        /// Data
        cout << file << "\t" <<
				opt_time/1000 << "\t" <<
				obj1stPhase << "\t" << 				
				endl;
        #endif
	}catch (IloException& e) {		
		///For iRace tests: <obj, time>
		cout << 99999999 << " " << time1stPhase/1000 << endl;		
		cerr << "Concert exception caught: " << e << endl;		
		e.end();
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}
	env.end();
}
