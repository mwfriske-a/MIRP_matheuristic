//~ #define NDEBUG
#include <assert.h>

#include "../include/util.h"
//~ #define zero_cost 1.0e-6
#define zero_cost 0
#define DecPrecision 1.0e+14
#define NPrintPortPoints

typedef IloArray<IloNumArray> NumNumMatrix;
typedef IloArray<IloIntArray> IntMatrix;

using namespace std;
using mirp::Instance;
using std::chrono::high_resolution_clock;

template <typename T>
double truncateBetween(T x, T xmin, T xmax) {
  return max(min(x, xmax), xmin);
}

template <typename T>
Timer<T>::Timer() : start_(high_resolution_clock::now()),
    running_(true)  { }

template <typename T>
Timer<T>::Timer(Timer<T>::R time_limit) : start_(high_resolution_clock::now()),
    running_(true), time_limit_(time_limit) { }

template <typename T>
void Timer<T>::start() { start_ = high_resolution_clock::now(); running_ = true; }

template <typename T>
void Timer<T>::stop() { finish_ = high_resolution_clock::now(); running_ = false; }

template <typename T>
typename Timer<T>::R Timer<T>::total() const {
  if(running_)
    return std::chrono::duration_cast<T>(high_resolution_clock::now() - start_).count();
  else
    return std::chrono::duration_cast<T>(finish_ - start_).count();
}

template <typename T>
bool Timer<T>::reachedTimeLimit() const {
  return total() > time_limit_;
}

template <typename T>
bool Timer<T>::isRunning() const {
  return running_;
}

template class Timer<std::chrono::milliseconds>;
template class Timer<std::chrono::seconds>;

void Instance::readInstance(IloEnv& env, const string& name){
	/// Reading Metadata
	ifstream meta;
	string fname = name + "metadata.txt";
	meta.open(fname);
	if (meta.fail()) {
      cerr << "Unable to open metada"<< endl;
      exit(1);
	}
	string s, s1;
	int i,num;
	
	getline(meta,s);
	stringstream(s) >> s1 >> t; 
	getline(meta,s); //jump num commodities
	getline(meta,s);
	stringstream(s) >> s1 >> loadReg;
	
	getline(meta,s);
	stringstream(s) >> s1 >> discReg;
	
	loadPorts = IloIntArray(env, loadReg);
	getline(meta,s);
	stringstream line(s);
	line >> s1;
	for (i=0;i<loadReg;i++){
		line >> num;
		loadPorts[i] = num;	
	}
	
	discPorts = IloIntArray(env,discReg);
	getline(meta,s);
	line.clear();
	line.str(s);
	line >> s1;
	for(i=0;i<discReg;i++){
		line >> num;
		discPorts[i]=num;		
	}
	
 	getline(meta,s);
	stringstream(s) >> s1 >> vC;
		
	v = IloIntArray(env, vC);
	getline(meta,s);
	line.clear();
	line.str(s);
	line >> s1;
	for(i=0;i<vC;i++){
		line >> num;
		v[i]=num;			
	}
	
	getline(meta,s);
	stringstream(s) >> s1 >> hoursPerPeriod;
	
	getline(meta,s);
	stringstream(s) >> s1 >> spotMarketPricePerUnit;
	
	getline(meta,s);
	stringstream(s) >> s1 >> spotMarketDiscountFactor;
	
	getline(meta,s);
	stringstream(s) >> s1 >> attemptCost;
	
	getline(meta,s);
	stringstream(s) >> s1 >> perPeriodRewardForFinishingEarly;
	
	getline(meta,s);
	stringstream(s) >> s1 >> constantForSinglePeriodAlphaSlack; 
	
	getline(meta,s);
	stringstream(s) >> s1 >> constantForCumulativeAlphaSlack;
	meta.close();
	///Ports 
	numTotalPorts = IloSum(loadPorts)+IloSum(discPorts);
	typePort = IloIntArray(env, numTotalPorts);
	idRegion = IloIntArray(env,numTotalPorts);
	identifyPort = IloArray<IntMatrix>(env, 2); // two types = loading(0) and discharging(1)
	identifyPort[0] = IntMatrix(env, loadReg);
	identifyPort[1] = IntMatrix(env, discReg);
	
	int j=0;
	// Loading region first
	for(i=0;i<loadReg;i++){ 
		identifyPort[0][i] = IloIntArray(env, loadPorts[i]);
		for(int k=0;k<loadPorts[i];k++){
			identifyPort[0][i][k] = j;
			typePort[j] = 0;
			idRegion[j] = i;
			j++;			
		}		
	}
	// Discharging region second
	identifyPort[1] = IntMatrix(env, discReg);
	for(i=0;i<discReg;i++){
		identifyPort[1][i] = IloIntArray(env, discPorts[i]);
		for(int k=0;k<discPorts[i];k++){
			identifyPort[1][i][k] = j;
			typePort[j] = 1;
			idRegion[j] = loadReg + i;
			j++;
		}					
	}
	
	//Init ports vectors	
	alp_max_j = IloNumArray(env, numTotalPorts);
	b_j = IloNumArray(env, numTotalPorts);
	s_j0 = IloNumArray(env, numTotalPorts);
	delta = IloIntArray(env, numTotalPorts);
	portFee = IloNumArray(env, numTotalPorts);
	x_coordinate = IloIntArray(env,numTotalPorts);
	y_coordinate = IloIntArray(env,numTotalPorts);
	
	//Inter ports
	distanceMatrix = NumNumMatrix(env,numTotalPorts);
	for (i=0;i<numTotalPorts;i++) distanceMatrix[i] = IloNumArray(env,numTotalPorts);
	
	//Init port-time vectors
	alp_max_jt = NumNumMatrix(env, numTotalPorts);
	d_jt = NumNumMatrix(env, numTotalPorts);
	dM_jt = NumNumMatrix(env, numTotalPorts);
	f_min_jt = NumNumMatrix(env, numTotalPorts);
	f_max_jt = NumNumMatrix(env, numTotalPorts);
	p_jt = NumNumMatrix(env, numTotalPorts);
	r_jt = NumNumMatrix(env, numTotalPorts);
	sMin_jt = NumNumMatrix(env, numTotalPorts);
	sMinM_jt = NumNumMatrix(env, numTotalPorts);
	sMax_jt = NumNumMatrix(env, numTotalPorts);
	sMaxM_jt = NumNumMatrix(env, numTotalPorts);
	int sumLoadPorts = IloSum(loadPorts);
	for(i=0;i<numTotalPorts;i++){
		alp_max_jt[i] = IloNumArray(env, t);
		d_jt[i] = IloNumArray(env, t);
		dM_jt[i] = IloNumArray(env, t+1); //CAUTION: index 0 is desconsidered
		f_min_jt[i] = IloNumArray(env, t);
		f_max_jt[i]  = IloNumArray(env, t);
		p_jt[i]  = IloNumArray(env, t);
		for(int time=0;time<t;time++){
			p_jt[i][time] = pow(spotMarketDiscountFactor,time)*spotMarketPricePerUnit;
		}
		r_jt[i] = IloNumArray(env, t);
		sMin_jt[i] = IloNumArray(env, t);
		sMinM_jt[i] = IloNumArray(env, t+1); //CAUTION: Considering 0 as the initial inventory.
		sMax_jt[i]  = IloNumArray(env, t);
		sMaxM_jt[i]  = IloNumArray(env, t+1); //CAUTION: Considering 0 as the initial inventory.
		if (i<sumLoadPorts)
			delta[i] = 1;
		else
			delta[i] = -1;
		//~ cout << i << " " << delta[i] << endl;
	}	
	
	///Loading ports 
	ifstream lPorts;
	fname.clear();
	fname = name + "loading_port_data.txt";
	lPorts.open(fname);
	if (lPorts.fail()) {
      cerr << "Unable to open loading ports data"<< endl;
      exit(1);
	}
	string s2,line1;
	//For each loading port
	for(j=0;j<IloSum(loadPorts);j++){
		getline(lPorts,s);
		stringstream(s) >> s1 >> s2; 
		//Partitioning the name of port
		vector<std::string>   namePort;
		stringstream ss(s2);
		while(getline(ss,line1,'_')){
			namePort.push_back(line1);
		}
		int idR = stoi(namePort[1]), idP = stoi(namePort[3]);
		int idPort = identifyPort[0][idR][idP];
		//~ cout << "Region " << idR << " Port " << idP << endl;
		for(i=0;i<4;i++) getline(lPorts,s); //Jump some info
		
		stringstream(s) >> s1 >> x_coordinate[idPort];
		getline(lPorts,s);
		stringstream(s) >> s1 >> y_coordinate[idPort];
		getline(lPorts,s);
		
		//~ cout << "Port " << identifyPort[1][idR][idP] << " x: " << x_coordinate[idPort] 
		//~ << " y: " << y_coordinate[idPort] << endl;
		
		stringstream(s) >> s1 >> portFee[idPort];
		//~ cout << "Fee " << portFee[idPort] << endl;
		getline(lPorts,s);
		stringstream(s) >> s1 >> b_j[idPort];
		//~ cout << "Berths " << b_j[idPort];
		//max and min amount per period -> considered for port-time nodes 
		getline(lPorts,s); 
		IloNum maxApp, minApp;
		stringstream(s) >> s1 >> maxApp;
		
		getline(lPorts,s); 	
		stringstream(s) >> s1 >> minApp;
		//~ cout << "Max/Min amount per period" << maxApp << "/"<<minApp <<endl;
		
		//Capacity
		getline(lPorts,s); 
		IloNum capacity;
		stringstream(s) >> s1 >> capacity;
		//~ cout << "Capacity " << capacity << endl;
		
		getline(lPorts,s); 
		stringstream(s) >> s1 >> s_j0[idPort];
		//~ cout << "Initial Inventory " << s_j0[idPort] << endl;
			
		getline(lPorts,s); // line of production
		stringstream(s) >> s1;
		ss.clear();
		ss.str(s);
		ss >> s1;
			
		//Storaging in the port-time pairs
		for (i=0;i<=t;i++){
			if(i<t){
				f_min_jt[idPort][i] = minApp;
				f_max_jt[idPort][i] = maxApp;		
				sMin_jt[idPort][i] = 0;
				sMax_jt[idPort][i] = capacity;		
				ss >> d_jt[idPort][i];
				alp_max_jt[idPort][i] = round(constantForSinglePeriodAlphaSlack*d_jt[idPort][i]);
				//~ r_jt[idPort][i] = 0;
				//~ cout << i << " " << d_jt[idPort][i] << endl;
				//~ cout << "Alpha j t " << alp_max_jt[idPort][i] << endl;
			}
			if(i==0){
				sMaxM_jt[idPort][i] = s_j0[idPort];
			}else{
				sMaxM_jt[idPort][i] = min(sMax_jt[idPort][i-1],  sMaxM_jt[idPort][i-1]+d_jt[idPort][i-1]);
				dM_jt[idPort][i]  = d_jt[idPort][i-1] + sMaxM_jt[idPort][i-1] - sMaxM_jt[idPort][i];
			}
			
		}
		
		alp_max_j[idPort] = constantForCumulativeAlphaSlack * d_jt[idPort][0];
		//~ cout << "Alpha J " << alp_max_j[idPort] << endl;
	}	
	
	///Discharging Ports
	ifstream dPorts;
	fname.clear();
	fname = name + "discharge_port_data.txt";
	dPorts.open(fname);
	if (dPorts.fail()) {
      cerr << "Unable to open discharging ports data"<< endl;
      exit(1);
	}
	//For each discharging port
	for(j=0;j<IloSum(discPorts);j++){
		getline(dPorts,s);
		stringstream(s) >> s1 >> s2; 
		//Partitioning the name of port
		vector<std::string>   namePort;
		stringstream ss(s2);
		while(getline(ss,line1,'_')){
			namePort.push_back(line1);
		}
		int idR = stoi(namePort[1]), idP = stoi(namePort[3]);
		int idPort = identifyPort[1][idR][idP];
		for(i=0;i<4;i++) getline(dPorts,s); //Jump some info
		
		stringstream(s) >> s1 >> x_coordinate[idPort];
		getline(dPorts,s);
		
		stringstream(s) >> s1 >> y_coordinate[idPort];		
		getline(dPorts,s);
		//~ cout << "Port " << idPort << " x: " << x_coordinate[idPort] 
		//~ << " y: " << y_coordinate[idPort] << endl;
		stringstream(s) >> s1 >> portFee[idPort];		
		
		getline(dPorts,s);
		stringstream(s) >> s1 >> b_j[idPort];
		
		getline(dPorts,s);
		IloNum maxApp, minApp; //storage the min and max amount to discharge
		stringstream(s) >> s1 >> maxApp;
		
		getline(dPorts,s);		 	
		stringstream(s) >> s1 >> minApp;
		
		getline(dPorts,s);
		IloNum capacity;
		stringstream(s) >> s1 >> capacity;
		
		getline(dPorts,s);
		stringstream(s) >> s1 >> s_j0[idPort];
		
		getline(dPorts,s); // line of consumption
		stringstream(s) >> s1;
		ss.clear();
		ss.str(s);
		ss >> s1;
		
		getline(dPorts,s); // line of revenue
		stringstream(s) >> s1;
		stringstream ss1(s);
		ss1 >> s1;
			
		//Storaging in the port-time pairs
		for (i=0;i<=t;i++){
			if (i<t){
				f_min_jt[idPort][i] = minApp;
				f_max_jt[idPort][i] = maxApp;		
				sMin_jt[idPort][i] = 0;			
				sMax_jt[idPort][i] = capacity;		
				ss >> d_jt[idPort][i];

				ss1 >> r_jt[idPort][i];
				alp_max_jt[idPort][i] = round(constantForSinglePeriodAlphaSlack*d_jt[idPort][i]);	
				//~ cout << round(constantForSinglePeriodAlphaSlack*d_jt[idPort][i]) << endl;
				//~ cout << "Alpha j t " << alp_max_jt[idPort][i] << endl;
			}
			if(i==0){
				sMinM_jt[idPort][i] = s_j0[idPort];	
			}else{
				sMinM_jt[idPort][i] = max(sMin_jt[idPort][i-1],  sMinM_jt[idPort][i-1]-d_jt[idPort][i-1]);
				dM_jt[idPort][i]  = d_jt[idPort][i-1] - sMinM_jt[idPort][i-1] + sMinM_jt[idPort][i];
			}
			//~ cout << "sMinM_jt " << idPort << " - " << i << " " << sMinM_jt[idPort][i] << "/ dMjt " << dM_jt[idPort][i] << endl;
		}
		alp_max_j[idPort] = constantForCumulativeAlphaSlack * d_jt[idPort][0];
		//~ cout << "Alpha J " << alp_max_j[idPort] << endl;
	}
	dPorts.close();
	
	#ifndef NPrintPortPoints	
	//Random distribute the ports between plan
	int minX,minY,maxX,maxY;
	minY = 100000;
	minX = 100000;
	maxX = -100000;
	maxY = -100000;
	for (j=0;j<numTotalPorts;j++){		
		if(y_coordinate[j] < minY)
			minY = y_coordinate[j];
		if(x_coordinate[j] < minX)
			minX = x_coordinate[j];
		if(y_coordinate[j] > maxY)
			maxY = y_coordinate[j];
		if(x_coordinate[j] > maxX)
			maxX = x_coordinate[j];
		
		//~ if(typePort[j]==0)
			//~ cout << "production\t\t";
		//~ else
			//~ cout << "consumption\t\t";
		//~ for(int t1=0;t1<t;t1++){
			//~ d_jt[j][t1] = d_jt[j][t1]*0.90; //percent of capacity
			//~ cout << d_jt[j][t1] << " ";
		//~ }
		//~ cout << endl;
			
		
	}
	//~ cout << "x: "<<minX << " " << maxX << endl;
	//~ cout << "y: "<<minY << " " << maxY << endl;
	srand(maxX); //seed value
	for (j=0;j<numTotalPorts;j++){
		x_coordinate[j] = iRand(minX,maxX);
		y_coordinate[j] = iRand(minY,maxY);
		cout << "Port " << j << endl <<
		"x_coordinate\t" << x_coordinate[j] << endl <<
		"y_coordinate\t" << y_coordinate[j]<< endl << endl;
	}
	
	//Print
	ofstream logPlot("ports_xy_New");	
	logPlot << "#Port\tx\ty\n";
	for(j=0;j<numTotalPorts;j++){
		logPlot << j << "\t" << x_coordinate[j] << "\t" << y_coordinate[j] << "\n";	
	}
	logPlot.close();	
	//~ exit(1);
	#endif
	///Distancematrix
	//~ ifstream distanceFile;
	//~ fname.clear();
	//~ fname = name + "distances.txt";
	//~ distanceFile.open(fname);
	//~ if (distanceFile.fail()) {
      //~ cerr << "Unable to open distances data"<< endl;
      //~ exit(1);
	//~ }
	int idI;
	//~ for (i=0;i<3;i++) getline(distanceFile,s);
	//~ for(j=0;j<numTotalPorts;j++){
		//~ stringstream ss(s);		
		//~ ss >> idI;		
		//~ for (i=0;i<numTotalPorts;i++){
			//~ ss >> distanceMatrix[idI][i];			
		//~ }
		//~ getline(distanceFile,s);
	//~ }
	//~ distanceFile.close();
	
	///Calculating the Euclidean distance between ports
	for (i=0;i<numTotalPorts;i++){
		for (j=0;j<numTotalPorts;j++){
			distanceMatrix[i][j] = sqrt(pow(x_coordinate[i]-x_coordinate[j],2) + pow(y_coordinate[i] - y_coordinate[j],2));
			//cout.precision(14);
			//cout << "D_" <<i << "," << j << " = " << distanceMatrix[i][j] << endl;
		}
	}
	
	///Additional calculus
	//Minimum amount to load/discharge in each Region 
	int sumLoadRegions = identifyPort[0].getSize();
	int sumDischRegions = identifyPort[1].getSize();
	f_min_r = IloNumArray(env,sumLoadRegions+sumDischRegions);
	int region=-1;	
	for(i=0;i<2;i++){ //Each type		
		for(int r=0;r<identifyPort[i].getSize();r++){ //Each region
			region++;			
			double min_f = DecPrecision;
			for(j=0;j<identifyPort[i][r].getSize();j++){ //Each port in region r
				if (f_min_jt[identifyPort[i][r][j]][0] < min_f)
					min_f = f_min_jt[identifyPort[i][r][j]][0];
			}
			f_min_r[region] = min_f;
		}
	}	
	
	///Vessel class
	ifstream vClass;
	fname.clear();
	fname = name + "vessel_class_data.txt";
	vClass.open(fname);
	if (vClass.fail()) {
      cerr << "Unable to open vessel class data"<< endl;
      exit(1);
	}
	//Auxiliar data for using after on the vessel data
	capacity = IloNumArray(env,vC);   
	IloNumArray avgSpeedInKnots(env,vC);   
	IloNumArray travelCostAsTermPerKm(env,vC);   
	IloNumArray discountTravelingEmpty(env,vC);   
	maxCapacity = 0;
	for (i=0;i<vC;i++){
		for (j=0;j<3;j++) getline(vClass,s);
		stringstream(s) >> s1 >> capacity[i];
		if(capacity[i] > maxCapacity)
			maxCapacity = capacity[i];
		
		getline(vClass,s);
		stringstream(s) >> s1 >> avgSpeedInKnots[i];
		
		getline(vClass,s);
		stringstream(s) >> s1 >> travelCostAsTermPerKm[i];
		
		getline(vClass,s);
		stringstream(s) >> s1 >> discountTravelingEmpty[i];		
	}
	vClass.close();
	
	///Vessel data
	ifstream vessel;
	fname.clear();
	fname = name + "vessel_data.txt";
	vessel.open(fname);
	if (vessel.fail()) {
      cerr << "Unable to open vessel data"<< endl;
      exit(1);
	}
	//Init vessels arrays
	int totalVessels = IloSum(v);
	s_v0 = IloNumArray(env, totalVessels);
	initialPort = IloIntArray(env, totalVessels);
	firstTimeAv = IloIntArray(env, totalVessels);
	q_v = IloNumArray(env, totalVessels);
	speed = IloNumArray(env, totalVessels);
	costKm = IloNumArray(env, totalVessels);
	trav_empt = IloNumArray(env, totalVessels);
	identifyVessel = IntMatrix(env,vC);
	idI = 0;
	for(i=0;i<vC;i++){
		identifyVessel[i] = IloIntArray(env,v[i]);
		for(j=0;j<v[i];j++){
			identifyVessel[i][j] = idI;			
			idI++;
		}		
	}
	
	for(i=0;i<totalVessels;i++){		
		for (j=0;j<3;j++) getline(vessel,s);
		int idClass;
		stringstream(s) >> s1 >> idClass;
		
		//Uptade vessel data from vessel class data 
		q_v[i] = capacity[idClass];
		speed[i] = avgSpeedInKnots[idClass];
		costKm[i] = travelCostAsTermPerKm[idClass];
		trav_empt[i] = discountTravelingEmpty[idClass];
		getline(vessel,s);
		getline(vessel,s);
		stringstream(s) >> s1 >> s_v0[i];

		getline(vessel,s);
		//Partitioning the name of initial port
		vector<std::string>   namePort;
		stringstream(s) >> s1 >> s2;
		stringstream ss(s2);
		while(getline(ss,line1,'_')){
			namePort.push_back(line1);
		}
		int typeReg;
		if (namePort[0].compare("LoadingRegion")==0) typeReg = 0;
		else typeReg = 1;
		int idR = stoi(namePort[1]), idP = stoi(namePort[3]);		
		initialPort[i] = identifyPort[typeReg][idR][idP];

		getline(vessel,s);
		stringstream(s) >> s1 >> firstTimeAv[i];
		//~ cout << "Vessel " << i << " class " << idClass << endl;
		//~ cout << "Capacity " << q_v[i] << endl;
		//~ cout << "Speed " << speed[i] << endl;
		//~ cout << "Km/cost " << costKm[i] << endl;
		//~ cout << "Trav. empty " << trav_empt[i] << endl;
		//~ cout << "Initial inventory " << s_v0[i] << endl;
		//~ cout << "Initial port " << initialPort[i] << endl;
		//~ cout << "First time avaliable " << firstTimeAv[i] << endl << endl;		
	}
	vessel.close();
	
	///Travel times for each vessel class	
	travelTime = IloArray<IntMatrix>(env,totalVessels);
	max_travelTime = IloNumArray(env,totalVessels);
	maxTravelTimeInstance = 0;
	for(i=0;i<totalVessels;i++){
		travelTime[i] = IntMatrix(env,numTotalPorts);
		for(j=0;j<numTotalPorts;j++) travelTime[i][j] = IloIntArray(env, numTotalPorts);
	}
	
	#ifdef NPrintPortPoints
	ifstream travelTimes;
	fname.clear();
	fname = name + "travelTimes.txt";
	travelTimes.open(fname);
	if (travelTimes.fail()) {
      cerr << "Unable to open travel times data"<< endl;
      exit(1);
	}
	getline(travelTimes,s);
	
	for(i=0;i<vC;i++){
		int maxTT = 0;
		for(j=0;j<2;j++) getline(travelTimes,s);
		int idPort;			
		for(int l=0;l<numTotalPorts;l++){			
			getline(travelTimes,s);			
			stringstream ss(s);
			ss >> idPort;			
			for(j=0;j<numTotalPorts;j++){
				int tTime;
				ss >> tTime;
				if (tTime > maxTT)
					maxTT = tTime;
				for (int k=0;k<v[i];k++){ 
					travelTime[identifyVessel[i][k]][idPort][j] = tTime;					
				}
			}				
		}
		for (int k=0;k<v[i];k++){
			max_travelTime[identifyVessel[i][k]] = maxTT;
			if(maxTT > maxTravelTimeInstance)
				maxTravelTimeInstance = maxTT;
		}
			
	}
	//~ for (i=0;i<totalVessels;i++){
		//~ cout << "Vessel " << i << endl;
		//~ for (j=0;j<numTotalPorts;j++){
			//~ cout << j << " ";
			//~ for (int k=0;k<numTotalPorts;k++){
				//~ cout << travelTime[i][j][k] << " " ;
			//~ }
			//~ cout << endl;
		//~ }
		//~ cout << endl;
	//~ }
	travelTimes.close();
	#endif
	//Calculating from the new coordinates
	#ifndef NPrintPortPoints
	for(i=0;i<vC;i++){
		for(int l=0;l<numTotalPorts;l++){		
			for(j=l;j<numTotalPorts;j++){
				for (int k=0;k<v[i];k++){ //For each vessel k on class i
					travelTime[identifyVessel[i][k]][l][j] = ceil(distanceMatrix[l][j]/(1.852*avgSpeedInKnots[i]*hoursPerPeriod));
					travelTime[identifyVessel[i][k]][j][l] = travelTime[identifyVessel[i][k]][l][j];
					//~ if (k==0)
						//~ cout << "Travel time between " << l << " and " << j << " Vessel class " << i << " " << travelTime[identifyVessel[i][k]][l][j] << endl;
				}				
			}			
		}
		//~ cout << endl;
	}	
	#endif
	
	///Additional information (for valid inequalities)
	lb_oper_jt = NumNumMatrix(env,numTotalPorts);
	delta_it = NumNumMatrix(env,numTotalPorts);
	for(int i=0;i<numTotalPorts;i++){
		lb_oper_jt[i] = IloNumArray(env,t);
		delta_it[i] = IloNumArray(env,t);
		p_delta_j = IloIntArray(env);
		for(int ti=0;ti<t;ti++){			
			double value=0;
			for(int u=0;u<=ti;u++){
				value += dM_jt[i][u];
			}
			lb_oper_jt[i][ti] = ceil(value/maxCapacity);
			//~ cout << i << " " << ti << ": "  << lb_oper_jt[i][ti] << endl;
			if(ti == 0)
				delta_it[i][ti] = lb_oper_jt[i][ti];
			else
				delta_it[i][ti] = lb_oper_jt[i][ti] - lb_oper_jt[i][ti-1];
			if(delta_it[i][ti] == 1)
				p_delta_j.add(ti);
		}		
	}
	
	
	///Building the arcs between port-time nodes 
	//For each vessel: 1 hash table for storage arc cost and 1 vector of strings to identify each arc
	arcs = IloArray< vector<string> >(env,totalVessels);
	c_va = IloArray<unordered_map<string,double> >(env, totalVessels);
	//For each vessel, node port-time
	inArcs = IloArray < IloArray<IntMatrix> > (env, totalVessels);
	inRegionArcs = IloArray < IloArray<IntMatrix> > (env, totalVessels);
	outArcs = IloArray < IloArray<IntMatrix> > (env, totalVessels);
	outRegionArcs = IloArray < IloArray<IntMatrix> > (env, totalVessels);

	//For each vessel and port
	maxTimeIntraReg = IloArray<IloIntArray>(env, totalVessels);
	
	for (i=0;i<totalVessels;i++){
		c_va[i] = unordered_map<string,double>();
		inArcs[i] = IloArray <IntMatrix>(env,numTotalPorts+1);	// Last index is used for the sink and source node
		outArcs[i] = IloArray <IntMatrix>(env,numTotalPorts+1); // Last index is used for the sink and source node
		outRegionArcs[i] = IloArray <IntMatrix>(env,numTotalPorts); // Only for regular nodes
		inRegionArcs[i] = IloArray <IntMatrix>(env,numTotalPorts); // Only for regular nodes
		maxTimeIntraReg[i] = IloIntArray(env, numTotalPorts);
		for(j=0;j<numTotalPorts;j++){
			inArcs[i][j] = IntMatrix(env,t);
			outArcs[i][j] = IntMatrix(env,t);
			outRegionArcs[i][j] = IntMatrix(env,t);
			inRegionArcs[i][j] = IntMatrix(env,t);
			for (int k=0;k<t;k++){
				inArcs[i][j][k] = IloIntArray(env);
				outArcs[i][j][k] = IloIntArray(env);
				outRegionArcs[i][j][k] = IloIntArray(env);
				inRegionArcs[i][j][k] = IloIntArray(env);
			}
			int maxTime = 0;
			for(int p=0;p<numTotalPorts;p++){
				if(idRegion[j] == idRegion[p]){
					if (travelTime[i][j][p] > maxTime)
						maxTime = travelTime[i][j][p];
				}
			}
		}
		//Special treatment for source and sink node
		inArcs[i][numTotalPorts] = IntMatrix(env,1);
		outArcs[i][numTotalPorts] = IntMatrix(env,1);
		inArcs[i][numTotalPorts][0] = IloIntArray(env);
		outArcs[i][numTotalPorts][0] = IloIntArray(env);		
				
		//~ addArc(i,0); //not used vessel arc		
		addArc(i,1); //Source arc(s)
		
		//Arcs from the initial Port and initial time available
		int j1 = initialPort[i];
		for (int t1=firstTimeAv[i];t1<t;t1++){
			//Travel arcs from initial port
			for(int j2=0;j2<numTotalPorts;j2++){
				if (j2 != j1){
					int t2 = t1 + travelTime[i][j1][j2];
					if (t2<t){
						addArc(i,2,j1,t1,j2,t2);
						//Create waiting arc for j2
						if (t2+1 < t) addArc(i,2,j2,t2,j2,t2+1);
						//Sink arc from port j2
						addArc(i,3, j2, t2);					
						//when the time t1 reach a time that can be considered as t2 for j2
						if (t1>=firstTimeAv[i]+travelTime[i][j1][j2]){
							addArc(i,2,j2,t1,j1,t2); //mirror arc to the initial port
							//Create arc from j2,t1 to others ports in time t3
							//~ for(int j3=0;j3<numTotalPorts;j3++){
								//~ if(j3 != j2 && j3 != j1){
									//~ int t3 = t1+travelTime[i][j2][j3];
									//~ if(t3<t) 
										//~ addArc(i,2,j2,t1,j3,t3);
								//~ }
							//~ }							
						}
						//Create arc from j2,t2 to others ports in time t3
						for(int j3=0;j3<numTotalPorts;j3++){
							if(j3 != j2 && j3 != j1){
								int t3 = t2+travelTime[i][j2][j3];
								if(t3<t) 
									addArc(i,2,j2,t2,j3,t3);
							}
						}						
					}
				}
			}
			//Waiting arcs of the initial port (j1)
			if (t1+1 < t) addArc(i,2,j1,t1,j1,t1+1);
			//Sink arc from port j1
			addArc(i,3, j1, t1);
		}	
	}
	
}
//type: 0 not used arc, 1 source arc, 2 regular arc(travel/waiting), 3 sink arc
void Instance::addArc(const int& vesselId, const int& type, const int& j1, const int& t1, const int& j2, const int& t2){
	string arc;	
	switch(type){	
	case 0:		
		arcs[vesselId].push_back("01");		
		c_va[vesselId].insert({"01",zero_cost});
		outArcs[vesselId][numTotalPorts][0].add(arcs[vesselId].size()-1); //exiting source
		inArcs[vesselId][numTotalPorts][0].add(arcs[vesselId].size()-1); //entering sink
		break;	
	case 1:
		//~ for (int t1=firstTimeAv[vesselId];t1<t;t1++){ //for each time after the time available of the vessel
			//~ arc.clear();
			//~ arc = sourceArc(vesselId,t1);
			//~ arcs[vesselId].push_back(arc);
			//~ c_va[vesselId].insert({arc,portFee[initialPort[vesselId]]});
			//~ inArcs[vesselId][initialPort[vesselId]][t1].add(arcs[vesselId].size()-1);
			//~ outArcs[vesselId][numTotalPorts][0].add(arcs[vesselId].size()-1);
		//~ }		
		//Only one source arc		
		arc = sourceArc(vesselId, firstTimeAv[vesselId]);
		arcs[vesselId].push_back(arc);
		c_va[vesselId].insert({arc,portFee[initialPort[vesselId]]});
		inArcs[vesselId][initialPort[vesselId]][firstTimeAv[vesselId]].add(arcs[vesselId].size()-1);
		inRegionArcs[vesselId][initialPort[vesselId]][firstTimeAv[vesselId]].add(arcs[vesselId].size()-1);
		outArcs[vesselId][numTotalPorts][0].add(arcs[vesselId].size()-1);
		break;
	case 2:
		//Waiting arc
		if (j1==j2){
			arc = travelArc(j1,t1,j2,t2);
			arcs[vesselId].push_back(arc);
			c_va[vesselId].insert({arc,zero_cost});
			outArcs[vesselId][j1][t1].add(arcs[vesselId].size()-1);
			inArcs[vesselId][j2][t2].add(arcs[vesselId].size()-1);			
		}else{ //Travel arc
			// j1 -> j2
			arc = travelArc(j1,t1,j2,t2);
			arcs[vesselId].push_back(arc);			
			double arc_cost;
			if (typePort[j1]==1 && typePort[j2]==0){ //If port j1 is consuming and j2 is a producer port
				arc_cost = costKm[vesselId]*distanceMatrix[j1][j2]*(1-trav_empt[vesselId]) + portFee[j2];		
			}else{
				arc_cost = costKm[vesselId]*distanceMatrix[j1][j2] + portFee[j2];
			}							
			c_va[vesselId].insert({arc,arc_cost});
			outArcs[vesselId][j1][t1].add(arcs[vesselId].size()-1);
			inArcs[vesselId][j2][t2].add(arcs[vesselId].size()-1);	
			//~ if (idRegion[j1] != idRegion[j2]) 
			if (typePort[j1] != typePort[j2])							
			//~ if (j1 != j2)							
				outRegionArcs[vesselId][j1][t1].add(arcs[vesselId].size()-1);
				inRegionArcs[vesselId][j2][t2].add(arcs[vesselId].size()-1);
				
			
			//~ //Mirror j2 -> j1
			//~ arc = travelArc(j2,t1,j1,t2);
			//~ arcs[vesselId].push_back(arc);			
			//~ c_va[vesselId].insert({arc,arc_cost});
			//~ outArcs[vesselId][j2][t1].add(arcs[vesselId].size()-1);
			//~ inArcs[vesselId][j1][t2].add(arcs[vesselId].size()-1);
		}
		break;
	case 3:
		arc = sinkArc(j1,t1);
		arcs[vesselId].push_back(arc);
		c_va[vesselId].insert({arc,-(t-t1-1)*perPeriodRewardForFinishingEarly});
		outArcs[vesselId][j1][t1].add(arcs[vesselId].size()-1);
		outRegionArcs[vesselId][j1][t1].add(arcs[vesselId].size()-1);
		inArcs[vesselId][numTotalPorts][0].add(arcs[vesselId].size()-1);
		break;
	default: 
		cout << "Erro: no arc type found for type = " << type << endl;
		exit(1);
		
	}
}

string portSegment(const int& j){	
	if (j > 9) 
		return to_string(j);
	else
		return "0" + to_string(j);		
}
string timeSegment(const int& t){	
	if (t > 99) 
		return to_string(t);
	else if(t>9)
		return "0" + to_string(t);		
	else if (t<=9)
		return "00"+to_string(t);
}

 string Instance::sourceArc(const int& vesselId, const int& t1){ //return a key for source arc 
	#ifndef NDEBUG
	assert(t1 >= firstTimeAv[vesselId]);
	#endif
	return "0" + portSegment(initialPort[vesselId]) + timeSegment(t1);	
}

 string Instance::sinkArc(const int& j1, const int& t1){ //return a key for source arc 
	return "1" + portSegment(j1) + timeSegment(t1);	
}

string Instance::travelArc(const int& j1, const int& t1, const int& j2, const int& t2){ //return a key for travel arc 
	//~ cout << j1 << "," << t1 << " -> " << j2 << "," << t2 << endl;
	#ifndef NDEBUG
	assert(t1<t2);
	#endif
	string arc;
	arc =  portSegment(j1) + timeSegment (t1) + portSegment(j2) + timeSegment(t2);
	//~ cout << "Travel arc " << arc << endl; 
	#ifndef NDEBUG
	assert(arc.size()==10);
	#endif
	
	return arc;
}
/*Returns the type of arc:
 0-Source arc; 
 1-Travel arc - travel from producer to consumer port
 2-Travel arc - travel from consumer to producer port
 3-Sink ARC - travel from producer to sink node
 4-Sink ARC - travel from consumer to sink node
 -1 - Waiting arc/travelling same region type
 -2 Otherwise (default)
 - otherwise
  */
int Instance::travelArcType(const std::string& arcName, int& timeJ1, int& timeJ2){
	int type = -2;
	timeJ1 = 999;
	timeJ2 = 999;
	if (arcName.size()==10){ //Regular arc
		int j1,j2;
		j1 = stoi(arcName.substr(0,2));
		j2 = stoi(arcName.substr(5,2));
		timeJ1 = stoi(arcName.substr(2,3)); //Getting the time from j1	
		timeJ2 = stoi(arcName.substr(7,3));
		assert(timeJ1 < t+1);
		if(typePort[j1] == 0 && typePort[j2] == 1)			
			return 1;
		else if(typePort[j1] == 1 && typePort[j2] == 0)
			return 2;
		else
			return -1; //waiting arc/travel same region type
	}
	else if(arcName.size()==6){ //Entering/exiting system
		if (stoi(arcName.substr(0,1)) == 1){ //Exiting system
			int j1 = stoi(arcName.substr(1,2));
			timeJ1 = stoi(arcName.substr(3,3)); //Getting the time from j1			
			assert(timeJ1 < t+1);
			if (typePort[j1] == 0)
				return 3;
			else
				return 4;
		
		}else{ //Entering system
			return 0;
		}
	}
	return type;	
}
int Instance::getArcType(const std::string& arcName, int& timeJ1, int& timeJ2, int& j1, int& j2){
	int type = -2;
	j1 = 999;
	j2 = 999;
	timeJ1 = 999;
	timeJ2 = 999;
	if (arcName.size()==10){ //Regular arc		
		j1 = stoi(arcName.substr(0,2));
		j2 = stoi(arcName.substr(5,2));
		timeJ1 = stoi(arcName.substr(2,3)); //Getting the time from j1	
		timeJ2 = stoi(arcName.substr(7,3));
		assert(timeJ1 < t+1);
		if(typePort[j1] == 0 && typePort[j2] == 1)			
			return 1;
		else if(typePort[j1] == 1 && typePort[j2] == 0)
			return 2;
		else
			return -1; //waiting arc
	}
	else if(arcName.size()==6){ //Entering/exiting system
		if (stoi(arcName.substr(0,1)) == 1){ //Exiting system
			j1 = stoi(arcName.substr(1,2));
			timeJ1 = stoi(arcName.substr(3,3)); //Getting the time from j1			
			assert(timeJ1 < t+1);
			if (typePort[j1] == 0)
				return 3;
			else
				return 4;
		
		}else{ //Entering system
			return 0;
		}
	}
	return type;	
}

///Instances arc verifications
bool Instance::shouldAddArc(const int& tS, const int& tF, const int& v, const int& a){
	int t1,t2, arcType;		
	arcType = travelArcType(arcs[v][a],t1,t2);
	//Verify the type of arc: if is sink arc add, travel and waiting arc must be into the interval tS and tF	
	switch (arcType)
	{	
		case 1 :
			if (t2 >= tS && t2 < tF) //Both sides of arc are in the interval of the appended block AND Arc between float and end block in the previous iteration
				return true;
			else
				return false;								
			break;
		case 2 :
			if (t2 >= tS && t2 < tF) //Both sides of arc are in the interval of the appended block AND Arc between float and end block in the previous iteration
			return true;
		else
			return false;								
			break;
		case -1 :
			if (t2 >= tS && t2 < tF) //Both sides of arc are in the interval of the appended block AND Arc between float and end block in the previous iteration
				return true;
			else
				return false;								
			break;
		case 3 :
		if (t1 >= tS && t1 < tF) // Depart node is in the model			
			return true;
		else
			return false;								
		break;
		case 4 :
		if (t1 >= tS && t1 < tF) // Depart node is in the model			
			return true;
		else
			return false;								
		break;
		default:			
			return false; //Source arcs should have already been added
			break;
	}		
}

bool Instance::isInModel(const int& tOEB, const int& v, const int& a){
	int t1,t2, arcType;		
	arcType = travelArcType(arcs[v][a],t1,t2);		
	if(arcType == 0 ||  //Source arcs - always add, assuming that no end block will contain a source arc
	  (t1 < tOEB && t2 < tOEB) ||   //Both sides of arc are out of the end block (considering travel and waiting arcs)
	  ((arcType == 3 || arcType == 4) && t1 < tOEB) )  //Sink arc for nodes in the MIP block
			return true;
	else
		return false;							
}

/*Return true if an arc should be relaxed when the model is constructed */
bool Instance::shouldRelaxArc(const int& tInteger, const int& v, const int& a){
	int t1,t2, arcType;		
	arcType = travelArcType(arcs[v][a],t1,t2);		
	assert(arcType !=-2);
	switch (arcType)
	{
		case 0: //Source arcs - never relax (they are already fixed)
			return false;
			break;
		case 1:
			if ((t1 >= tInteger ) ||  //Both sides of arc are out the integer block (TODO t2 > t1 - remove 2nd term)
				(t1 < tInteger && t2 >= tInteger))  //Arc travessing between integer and relaxed block is relaxed
				return true;
			else
				return false;								
			break;
		case 2:
			if ((t1 >= tInteger ) ||  //Both sides of arc are out the integer block (TODO t2 > t1 - remove 2nd term)
				(t1 < tInteger && t2 >= tInteger))  //Arc travessing between integer and relaxed block is relaxedock is relaxed
				return true;
			else
				return false;				
			break;
		case 3: //sink arc 
			if(t1 >= tInteger) //- If (t1 > tInteger) then the last port-Time will have a sink arc integer (but the exiting arc will be float)
				return true;
			else
				return false;
			break;
		case 4: //sink arc
			if(t1 >= tInteger)
				return true;
			else
				return false;
			break;
		case -1:  //Waiting arc
			if ((t1 >= tInteger) ||  //Both sides of arc are out the integer block (TODO t2 > t1 - remove 2nd term)
				(t1 < tInteger && t2 >= tInteger))  //Arc travessing between integer and relaxed block is relaxed
				return true;
			else
				return false;				
			break;
		default:
			cout << "ERROR! Impossible evaluate arc for relax or not" << endl;
			exit(1);
			break;			
	}

}
bool Instance::shouldIntegralizeArc(const int& v, const int& a, const int& tS, const int& tF){
	int t1, t2;				
	int arcType = travelArcType(arcs[v][a], t1, t2);			
		switch (arcType)
		{
			case 1: //Integralize travel and waiting arcs between tS and tF and before tS that reach until tF
				if ((t1 >= tS && t2 < tF) ||  //Both sides of arc are in the interval 
				(t1 < tS && (t2 >= tS && t2 < tF)) )//Arc between float and end block in the previous iteration
					return true;
				else
					return false;								
			break;
			case 2: //Integralize travel and waiting arcs between tS and tF and before tS that reach until tF
				if ((t1 >= tS && t2 < tF) ||  //Both sides of arc are in the interval 
				(t1 < tS && (t2 >= tS && t2 < tF)) )//Arc between float and end block in the previous iteration
					return true;
				else
					return false;
				break;
			case -1: //Integralize travel and waiting arcs between tS and tF and before tS that reach until tF
				if ((t1 >= tS && t2 < tF) ||  //Both sides of arc are in the interval 
				(t1 < tS && (t2 >= tS && t2 < tF)) )//Arc between float and end block in the previous iteration
					return true;
				else
					return false;
				break;
			case 3:  //Sink arc
				if(t1 >= tS && t1 < tF)
					return true;
				else
					return false;
				break;
			case 4:  //Sink arc
				if(t1 >= tS && t1 < tF)
					return true;
				else
					return false;
				break;
			default:
				return false;
				break;
		}
}



double fRand(double fMin, double fMax,const int& seed = -1){
	if (seed != -1)
		srand(seed);
	else
		srand(time(NULL));
		
	double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
int iRand(int iMin, int iMax){
	//~ srand(seed);
	return iMin + (rand() % (int)(iMax - iMin + 1));
	
}
