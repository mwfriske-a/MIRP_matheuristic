#pragma once
#include <vector>
#include <utility>
#include <ilcplex/ilocplex.h>

#include "util.h"

namespace mirp{		
	typedef IloArray<IloNumVarArray> NumVarMatrix;
	typedef IloArray<IloIntVarArray> IntVarMatrix;
	typedef IloArray<IloNumArray> NumNumMatrix;
	typedef IloArray<IloIntArray> IntMatrix;
	typedef IloArray<IloBoolVarArray> BoolMatrix; 
	
	void milp(std::string file, const double& timeLimit, std::string optStr);
	
	void fixAndRelax(std::string file, std::string optStr, const double& nIntervals, const double& gapFirst, const int& f, const double& overLap, const int& endBlockFirst,
const int& timePerIter, const double& mIntervals, const int& timePerIter2, const double& gapSecond, const double& overlap2,const double& timeLimit);
	
	void fixAndRelaxH(std:: string file, std::string optStr, const double& gapFirst, const int& outVessels, const int& timeLimitFirst, 
	const double& mIntervals, const int& timeLimitSecond, const double& gapSecond,const double& overlap2);
	
	void fixAndRelaxHV(std::string file, const double& nIntervals, const int& f, const double& overLap, const int& endBlock,
	const double& gapFirst, const int& outVessels, const int& timeLimitFirst, 
	const double& mIntervals, const int& timeLimitSecond, const double& gapSecond);
	
	struct Model{
		IloObjective obj;
		IloModel model;
		IloCplex cplex;	
			
		/*Variables
		 * OBS: Arc variables has indexes for port 0 and port J+1, corresponding to the source and sink 'ports'.
		 * This not happen with 
		 */
		NumVarMatrix alpha; 								//amount of product purchasses from or sells to spot market in time time period	
		NumVarMatrix beta;	 								//Slack for amount of product purchasses from or sells to spot market in time time period	
		NumVarMatrix theta;	 								//Slack for amount of product purchasses from or sells to spot market in time time period	
		NumVarMatrix sP;									//Number of units of inventory at port j available at the end of time period t
		IloArray<NumVarMatrix> f;		 					//amount loaded/discharged by a vessel in a port-time node
		
		IloArray<IloArray<NumVarMatrix> >fX;				//Load on board vessel when travelling from port i to port j, leaving at time period t
		IloArray<NumVarMatrix> fOA;							//Load on board vessel when starting to operate at port i in time period t
		IloArray<NumVarMatrix> fOB;							//Load on board vessel before continuing to operate at port i in time period t
		IloArray<NumVarMatrix> fW;							//Load on board vessel while waiting at port i in time period t
		IloArray<NumVarMatrix> fWB;							//Load on board vessel while waiting at port i in time period t after operated in time t-1
		
		IloArray<IloArray<IloArray<IloBoolVarArray> > > x;	//Takes value 1 if vessel v travesses an arc departing from i and arriving at j, starting in time period t. Port 0 corresponds to the source node and port T+1 corresponds to the sink node
		IloArray<IloArray<IloBoolVarArray> > z;				//Takes value 1 if vessel v attmpts to load/discharge at port j in time period t
		IloArray<IloArray<IloBoolVarArray> > w;				//Takes value 1 if vessel v waits at port j in time period t
		IloArray<IloArray<IloBoolVarArray> > wB;			//Takes value 1 if vessel v waits at port j in time period t after operated in time period t-1
		IloArray<IloArray<IloBoolVarArray> > oA;			//Takes value 1 if vessel v starts to operate at port i in time t
		IloArray<IloArray<IloBoolVarArray> > oB;			//Indicates the succeding operation at port i vessel v starts to operate at port i in time t
		
		//Branching variables 
		IloArray<IloIntVarArray> sX;						//Number of times a port i is visited by vessel v during all planing horizon
		IloArray<IloIntVarArray> sOA;						//Number of times a vessel v starts to operate (or operate) at port i
		//~ IloArray<IloBoolVarArray> w;			//Used to force a vessel to depart from a region after it is empty(or full)
	
		//Information
		IloArray<IloArray<IloArray<IloIntArray> > > hasArc; 	//Takes value 1 if there is a travelling arc for vessel v travelling from i to j, starting at time t
		IloArray<IloArray<IloArray<IloNumArray> > > arcCost;	//Cost of travelling arc if hasArc == 1
		IloArray<IloArray<IloIntArray> > hasEnteringArc1st;		//Takes value 1 if node v,i,t has an entering arc in the first level(source,traveling,or waiting). i \in J
		
		//Conversion - For relax-and-fix 
		IloArray<IloArray<IloArray<IloArray<IloConversion> > > > convertX; //Convert x variables according [v][i][j][t]
		IloArray<IloArray<IloArray<IloConversion> > > convertZ; //Convert z variables according to time ([v][i][t])
		IloArray<IloArray<IloArray<IloConversion> > > convertW; //Convert w variables according to time ([v][i][t])
		IloArray<IloArray<IloArray<IloConversion> > > convertWB; //Convert wB variables according to time ([v][i][t])
		IloArray<IloArray<IloArray<IloConversion> > > convertOA; //Convert oA variables according to time ([v][i][t])
		IloArray<IloArray<IloArray<IloConversion> > > convertOB; //Convert oB variables according to time ([v][i][t])
		
		//Store variables values 
		IloArray<IloArray<IloArray<IloNumArray> > > xValue;
		IloArray<IloArray<IloNumArray> > zValue;
		IloArray<IloArray<IloNumArray> > wValue;
		IloArray<IloArray<IloNumArray> > wBValue;
		IloArray<IloArray<IloNumArray> > oAValue;
		IloArray<IloArray<IloNumArray> > oBValue;
		
		IloArray<IloArray<IloNumArray> > fValue;
		IloArray<IloArray<IloArray<IloNumArray> > >  fXValue;
		IloArray<IloArray<IloNumArray> > foAValue;
		IloArray<IloArray<IloNumArray> > foBValue;
		IloArray<IloArray<IloNumArray> > fWValue;
		IloArray<IloArray<IloNumArray> > fWBValue;
		IloArray<IloNumArray> sPValue;
		IloArray<IloNumArray> alphaValue;

		///Constraints
		//Balance between nodes
		IloRangeArray sinkNodeBalance;
		IloRangeArray sourceNodeBalance;
		IloArray<IloArray<IloRangeArray> > firstLevelBalance;
		IloArray<IloArray<IloRangeArray> > secondLevelBalance;
		IloArray<IloArray<IloRangeArray> > linkBalance;
		
		//Flow vessels
		IloArray<IloArray<IloRangeArray> > firstLevelFlow;
		IloArray<IloArray<IloRangeArray> > secondLevelFlow;

		IloArray<IloRangeArray> berthLimit;
		IloArray<IloArray<IloArray<IloRangeArray> > >  travelAtCapacity;
		IloArray<IloArray<IloArray<IloRangeArray> > >  travelEmpty;
		IloArray<IloRangeArray> portInventory;
		IloRangeArray cumSlack;
		IloArray<IloArray<IloRangeArray> > operationLowerLimit;
		IloArray<IloArray<IloRangeArray> > operationUpperLimit;
		
		IloArray<IloArray<IloArray<IloRangeArray> > >  flowCapacityX;
		IloArray<IloArray<IloRangeArray> >   flowCapacityOA;
		IloArray<IloArray<IloRangeArray> >   flowCapacityOB;
		IloArray<IloArray<IloRangeArray> >   flowCapacityW;
		IloArray<IloArray<IloRangeArray> >   flowCapacityWB;
		IloArray<IloArray<IloArray<IloRangeArray> > >  flowMinCapacityX;
		IloArray<IloArray<IloRangeArray> >   flowMinCapacityOA;
		IloArray<IloArray<IloRangeArray> >   flowMinCapacityOB;
		IloArray<IloArray<IloRangeArray> >   flowMinCapacityW;
		IloArray<IloArray<IloRangeArray> >   flowMinCapacityWB;
		
		///Valid inequalities
		//For multiples denominators
		//~ //Loading ports
		//~ IloArray<IloArray<IloRangeArray> >knapsack_P_1; //Inequalities for the case where T_v = T (for each Q_v, j and it)
		//~ IloArray<IloRangeArray> knapsack_P_2; //Inequalities for the case where T_v = \emptyset
		//~ //Dischargin ports
		//~ IloArray<IloArray<IloRangeArray> > knapsack_D_1; //Case R²_v = T 
		//~ #ifndef WaitAfterOperate
		//~ IloArray<IloArray<IloRangeArray> > knapsack_D_2; //Case R¹_v = T
		//~ IloArray<IloRangeArray> knapsack_D_3; //Case R⁰_v = T
		//~ #endif
		//~ #ifdef WaitAfterOperate
		//~ IloArray<IloRangeArray>  knapsack_D_2; //Case R¹_v = T
		//~ #endif
		
		//For one denominator
		//Loading ports
		IloArray<IloRangeArray> knapsack_P_1; //Inequalities for the case where T_v = T (for each Q_v, j and it)
		IloArray<IloRangeArray> knapsack_P_2; //Inequalities for the case where T_v = \emptyset
		//Dischargin ports
		IloArray<IloRangeArray> knapsack_D_1; //Case R²_v = T
		IloArray<IloRangeArray> knapsack_D_2; //Case R¹_v = T
		IloArray<IloRangeArray> knapsack_D_3; //Case R⁰_v = T

		//Additional constraints;
		IloArray<IloRangeArray> branchingX;		//Branching rule where sX = sum{j\inJ}{t\inT} x[v][i][j][t]
		IloArray<IloRangeArray> branchingOA;	//Branching rule where sOA = sum{j\inJ}{t\inT} oA[v][i][t]
		IloRangeArray y_Sum; 					//For branching (total number visits in port j during all planing horizon)
		IloArray<IloRangeArray> minVisits;		//Valid inequality - minimum number of operations that must be done in each port from t=0 until time t' 
		IloArray<IloArray<IloRangeArray> > operateAndDepart;	//Ensure that a vessel must exit a region after operate when the vessel capacity is lesser than the maximum operation at port

		//Wagner-whitin e Lot-sizing with constant capacity (1 level)- variables and constraints
		IloArray<IloRangeArray> wwcc_relaxation;
		
		//Lot-sizing with start up relaxation
		IloArray<IloRangeArray> startup_sumStartIfOperate;  //(i,t)
		IloArray<IloRangeArray> startup_sumForceStartup; 	//(i,t)
		IloArray<IloRangeArray> startup_dlsccs;				//(i,t)
		IloArray<IloRangeArray> startup_validInequality;	//(i,x), where x is a combination in funcion of [k,l] and p
		
	
		///Branching rules
		IntVarMatrix sumX;
		IloArray<IloRangeArray> priorityX;
		IntVarMatrix sumOA;
		IloArray<IloRangeArray> priorityOA;
		
		//Others
		std::vector<std::pair<int,int> > ordV;
			
		Model(IloEnv& env) : obj(env), model(env), cplex(model){
			model.add(obj);		
		};	
		void buildModel(IloEnv& env, Instance inst); 
		void buildFixAndRelaxModel(IloEnv& env, Instance inst, const double& nIntervals, const int& endBlock); 
		void buildFixAndRelaxHorizontalModel(IloEnv& env, Instance inst, const int& outVessels); 
		void buildFixAndRelaxHVModel(IloEnv& env,Instance inst, const double& nIntervals, const int& endBlock, const int& outVessels);
		void setParameters(IloEnv& env, const double& timeLimit, const double& gap);
		void fixSolution(IloEnv& env, Instance inst, const int& t3S, const int& t3F, const int& p, const bool& fixSinkArc);
		void modifyModel(IloEnv& env, Instance inst, const int& nIntervals, const int& tS_fix, const int& tF_fix, const int& tS_add, const int& tF_add, const int& tS_int, const int& tF_int);
		void fixAllSolution(IloEnv& env, const Instance& inst);
		void decreaseEndBlock (IloEnv& env, Instance inst, const double& nIntervals, const int& t2S, const int& t2F); //Add to the model FO and constranints the variables between t2S and t2F
		void reIntegralize(IloEnv& env, Instance inst, const int& t1S, const int& t1F);
		void improvementPhase_timeIntervals(IloEnv& env, Instance inst, const double& mIntervals, const double& timePerIter, const double& gap, const double& overlap, Timer<std::chrono::milliseconds>& timer_cplex,float& opt_time,const double& timeLimit, float& elapsed_time, double& incumbent);
		void improvementPhaseVND_timeIntervals(IloEnv& env, Instance inst, const double& mIntervals, const double& timePerIter, const double& gap, const double& overlap, Timer<std::chrono::milliseconds>& timer_cplex,float& opt_time,const double& timeLimit, float& elapsed_time, double& incumbent);
		void unFixInterval(Instance inst, const int& tS, const int& tF);
		void improvementPhase_vessels(IloEnv& env, Instance inst, const double& timePerIter, const double& gap, double& obj1stPhase, Timer<std::chrono::milliseconds>& timer_cplex,float& opt_time,const double& timeLimit, float& elapsed_time);
		void fixVessel(IloEnv env,Instance inst, const int& v);
		void fixVesselPair(IloEnv env, Instance inst, const int& v,const int& v1);
		void fixVesselInterval(IloEnv env,Instance inst, const int& v,const int& tS, const int& tF);
		void unFixVessel(Instance inst, const int& v);
		void polish(IloEnv& env, Instance inst,const double& timeLimit, const double& gap);
		void printSolution(IloEnv env, Instance inst, const int& tF);
		void getSolVals(IloEnv& env, const Instance& inst);
		void getSolValsW(IloEnv& env, Instance inst, const int& tS, const int& tF, const bool& fixSinkArc);
		void addVesselToModel(IloEnv& env,Instance inst, const int& vAdd);
		void addVesselIntervalToModel(IloEnv& env,Instance inst, const int& vAdd, const double& nIntervals, const int& tS, const int& tF, const int& outVessels, const int& itInt);
		void integralizeVessel(const Instance& inst, const int& vInt);
		void integralizeVesselInterval(Instance inst, const int& vInt, const int& tS, const int& tF);
		void resetObjFunction(IloEnv& env, Instance inst);
		void fixVesselLessInterval(IloEnv env, Instance inst, const int& v, const int& tS, const int& tF);
		void improvementPhase_intervalVessel(IloEnv& env, Instance inst, const double& mIntervals, const double& timePerIter, 
			const double& gap, const double& overlap, Timer<std::chrono::milliseconds>& timer_cplex,float& opt_time, 
			const double& timeLimit, float& elapsed_time, double& incumbent);
		
		void improvementPhase_typePortsLS(IloEnv env, Instance inst, const double& timePerIter, const int& gap,Timer<std::chrono::milliseconds>& timer_cplex, float& opt_time, const double& timeLimit, float& elapsed_time, double& incumbent);
		void warmStart(IloEnv env, Instance inst, const double& timePerIter);
		void addInventoryConstraints(Instance inst, IloEnv& env, bool& feasible);
	};
}
