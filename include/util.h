#pragma once
#include <chrono>
#include <unordered_map>
#include <ilcplex/ilocplex.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <random>

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloNumArray> NumNumMatrix;
typedef IloArray<IloIntArray> IntMatrix;

template <typename T>
double truncateBetween(T x, T xmin, T xmax);

template<typename T>
class Timer{
  typedef typename T::rep R;
public:
	Timer();
  Timer(R time_limit);
	void start();
	void stop();
  bool isRunning() const;
	R total() const;
  bool reachedTimeLimit() const;
private:
  std::chrono::high_resolution_clock::time_point start_, finish_;
  bool running_;
  R time_limit_;
};

bool DBL_EQL(double, double);
bool DBL_GE (double, double);
bool DBL_G (double, double);
bool DBL_L (double, double);

double fRand(double fMin, double fMax,const int& seed);
int iRand(int fMin, int fMax);

extern int seed1;
extern std::default_random_engine engine;

namespace mirp{
	
struct Instance{
//Random 
int iRand2(int iMin, int iMax);	
float fRand2(float fMin, float fMax);
//~ std::default_random_engine generator;
///Parameters
//Metadata
IloInt t; 					//# Time periods
IloInt hoursPerPeriod;		//# hors per period
IloInt vC; 					//# Vessel classes
IloIntArray v;		 		//# vessels per class
IloInt loadReg; 			//# loading regions
IloInt discReg;				//# unloading regions 
IloIntArray loadPorts; 		//# loading ports of each region
IloIntArray discPorts;	 	//# unloading ports of each region
IloNum spotMarketPricePerUnit;	
IloNum spotMarketDiscountFactor;
IloNum attemptCost;
IloNum perPeriodRewardForFinishingEarly;
IloNum constantForSinglePeriodAlphaSlack;
IloNum constantForCumulativeAlphaSlack;     
int numTotalPorts;

IloNumArray capacity;

//Ports (loading and discharging and regions merged)
IloNumArray alp_max_j; 		//upper bound on the cummulative amount of product that can be sold to/bought from the spot market at port j
IloNumArray b_j;			//# of berths
IloNumArray s_j0;			// initial inventory
IloIntArray delta;			// +1 if loading port, -1 if discharging port
IloNumArray portFee;		// Port fee
IloIntArray x_coordinate;	// X coordinate
IloIntArray y_coordinate;	// Y coordinate

//Identifier - just for reading instance
IloArray<IntMatrix> identifyPort;	// Identify for a loading or unloading region, what is the index of port in the parameters above. [l(d)][r][p] where l(d) is the loading(dischargin), r is the region id and p is the port id
IloIntArray typePort; //Opposite of the identifyPort, for each port if id=0 is loading port,  and id=1 is unloading port

//Inter ports
NumNumMatrix distanceMatrix;		//Distance between port i and j;

//Regular nodes - port time pairs (j,t) 
NumNumMatrix alp_max_jt; 	//upper bound on the amount of product that can be sold to/bought fromthe spot market at port j in time period t
NumNumMatrix d_jt;			//# units produced/consumed at port j in time t
NumNumMatrix dM_jt;			//# units produced/consumed at port j in time t (depends of sMinM_jt and sMaxM_jt)
NumNumMatrix f_min_jt;		// minimum amount of product that can be loaded/discharged at port j in time period t;
NumNumMatrix f_max_jt;		// max amount of product that can be loaded/discharged at port j in timeperiod t;
NumNumMatrix p_jt; 			// penality parameter associated  with one unit of lost production or stockout at port j in time period t;
NumNumMatrix r_jt;			// the unit sales revenue for product discharged at por j in time period t (equals 0 in the loading ports)
NumNumMatrix sMin_jt;		// lower bound of port j in time period t
NumNumMatrix sMinM_jt;		// modified lower bound of port j in time period t - Only for discharging ports
NumNumMatrix sMax_jt;		// capacity of port j in time period t
NumNumMatrix sMaxM_jt;		// modified upper bound of port j in time period t - Only for loading ports 

//For disxrete lot-sizing relaxation with constrait capacity and startup
NumNumMatrix lb_oper_jt;	// Lower bound on the number of operating periods needed during the first t periods (depends of dM_jt) - = \tilde{O}_{it} defined in Agra et al 2013
NumNumMatrix delta_it;		// Vector of diference of lb_oper_jt and lb_oper_j,t-1
IloIntArray p_delta_j;		//Storage for each i \in J the index of delta_it in which delta_it = 1 (we will use the size of p_delta_j)

//Arcs
IloArray< std::vector<std::string> >arcs;				// For each vessel, a vector keys that pointing for register in a hash table
IloArray <IloArray<IntMatrix> >inArcs;					// For each vessel and each node n=(j,t), subset of entering arcs in the node (index of arcs)
IloArray <IloArray<IntMatrix> >inRegionArcs;			// For each vessel and each node n=(j,t), subset of entering arcs in the node which depart from another region(index of arcs)
IloArray <IloArray<IntMatrix> >outArcs;					// For each vessel and each node n=(j,t), subset of exiting arcs in the node  (index of arcs)
IloArray <IloArray<IntMatrix> >outRegionArcs;			// For each vessel and each node n=(j,t), subset of exiting arcs in the node which arrive in another region (index of arcs)
IloArray<std::unordered_map<std::string,double> > c_va; //Cost for vessel v travesse an arc a

//Additional information
IloNumArray f_min_r; //Min f_min_jt per region R (assuming that f_min_r is fixed for all t \in T)

//Vessel - classes are merged
IloNumArray s_v0;			// Initial inventory at vessel v of class c
IloIntArray initialPort;	// Initial port of vessel v
IloIntArray firstTimeAv;	// First time available vessel
//Equal for each vessel in the same class
IloNumArray q_v; 			// Capacity of vessel v 
IloNumArray speed;			// Speed of vessel v
IloNumArray costKm;			// Cost per km of vessel v
IloNumArray trav_empt;		// Discont factor for traveling empty
int maxVesselCapacity;			// Maximum capacity between vessel classes;
IloNumArray max_travelTime;	// Maximum travel time between ports of one vessel
int maxTravelTimeInstance;	//Maximum travel time of an instance
//Identifiers
IntMatrix identifyVessel;
IloIntArray idRegion; //Storage the region id of each vessel j
//Vessel x interpoints
IloArray<IntMatrix> travelTime;	//Travel time between  port i and j for vessel v 

//Additional data
IloArray<IloIntArray> maxTimeIntraReg; 		//Maximum travel time between port j and any other port of same region for vessel v


Instance(IloEnv& env){};
void readInstance(IloEnv& env, const std::string& name);
std::string sourceArc(const int& vesselId, const int& t1);
std::string travelArc(const int& j1, const int& t1, const int& j2, const int& t2);
std::string sinkArc(const int& j1, const int& t1);
void addArc(const int& vesselId, const int& type, const int& j1=0, const int& t1=0, const int& j2=0, const int& t2=0);
int travelArcType(const std::string& arcName, int& timeJ1, int& timeJ2);
int getArcType(const std::string& arcName, int& timeJ1, int& timeJ2, int& j1, int& j2);
bool isInModel(const int& tWEB, const int& v, const int& a);	//FOR RELAX-AND-FIX : Return true if an arc is in the model (and out of the end block), false otherwise
bool shouldRelaxArc(const int& tInteger, const int& v, const int& a);  //FOR RELAX-AND-FIX: Return true if an arc must be relaxed given tInteger as the last t of integer block
bool shouldAddArc(const int& tS, const int& tF,const int& v, const int& a); //For Relax-and-fix: Retunr true if arc should be inserted on the model while end block is reduced 
bool shouldIntegralizeArc(const int& v, const int& a, const int& t1S, const int& t1F);
void perturb(IloEnv& env, const std::string name, const int& numInstances);
};
}
