#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include "../include/util.h"
#include "../include/model.h"
// #define NTestInstanceSet
using namespace std;

int main(const int argc, char** argv) {
 if (argc < 7) {
	cout << "Error! At least 4 parameters must be defined" << endl;	
 }
stringstream ss;
ss << "usuage is \n -v : Version (1 - for MILP; 2- for fix-and-relax vertically; 3-fix-and-relax horizontally) \n " 
	<< "-p : instance path \n " 
	<< "-t : time limit in seconds - for -v = 2, corresponds to timelimit of local searches \n " 
	
	<< "-q : 1 if use valid inequalities,0 otherwise \n"
	<< "-c : 1 if use additional constraints, 0 otherwise \n"
	<< "-d : 1 if use tightening of inventory contraints, 0 otherwise (only when e > 0) \n"
	<< "-k : 1 if use of proportional alphaMaxJ, 0 otherwise (only when e > 0) \n"
	<< "-z : 1 if simplify instance, 0 otherwise \n"
	<< "-b : 1 if bounds on flow fOB(fWB) and fX are tighten\n"
	
	<< "For -v = {2,3} (Fix and relax algorithm):\n " 
	<< "-n : Number of intervals in the first phase (-v =2) \n "
	<< "-g : Initial GAP for solving the first phase (-v = 2,3) \n"
	<< "-f : Number of intervals that leaves the end block in each iteration (-v =2) \n "
	<< "-o : percentage of overlapping between the fixed and the integer interval (-v =2) \n "
	<< "-e : Size of the end block/vessels out of model (must be at most n-2) \n "
	<< "-s : A string defining the fix-and-optimize procedures and their order to be used in the second phase \n" 
	<< "-r : Vessels out of model (must be at most v-2) \n "
	<< "-l : Time limit for solving each interval \n "
	<< "-m : Number of intervals in the improvement phase (-v = 2) - use m=0 for just produce a  R&F solution \n "
	<< "-i : Time limit for one loop in the improvement phase \n "
	<< "-h : Fixed GAP in the imnprovement phase \n "
	<< "-u : Overlap in the imnprovement phase \n ";
int modelId;
stringstream file;
string fixOptStr;
double timeLimit, gapFirst, gapSecond, overLap, mIntervals, nIntervals, overlap2;
int opt, endBlock, outVessel, timePerIterFirst, timePerIterSecond, f ;
bool validIneq, addConstr, tightenInvConstr, proportionalAlpha, reduceArcs, tightenFlow;
vector <string> instances;
instances.push_back("../Group1_data_format_only_files/LR1_1_DR1_3_VC1_V7a/");
instances.push_back("../Group1_data_format_only_files/LR1_1_DR1_4_VC3_V11a/");
instances.push_back("../Group1_data_format_only_files/LR1_1_DR1_4_VC3_V12a/");
instances.push_back("../Group1_data_format_only_files/LR1_1_DR1_4_VC3_V12b/");
instances.push_back("../Group1_data_format_only_files/LR1_1_DR1_4_VC3_V8a/");
instances.push_back("../Group1_data_format_only_files/LR1_1_DR1_4_VC3_V9a/");
instances.push_back("../Group1_data_format_only_files/LR1_2_DR1_3_VC2_V6a/");
instances.push_back("../Group1_data_format_only_files/LR1_2_DR1_3_VC3_V8a/");
instances.push_back("../Group1_data_format_only_files/LR2_11_DR2_22_VC3_V6a/");
instances.push_back("../Group1_data_format_only_files/LR2_11_DR2_33_VC4_V11a/");
instances.push_back("../Group1_data_format_only_files/LR2_11_DR2_33_VC5_V12a/");
instances.push_back("../Group1_data_format_only_files/LR2_22_DR2_22_VC3_V10a/");
instances.push_back("../Group1_data_format_only_files/LR2_22_DR3_333_VC4_V14a/");
instances.push_back("../Group1_data_format_only_files/LR2_22_DR3_333_VC4_V17a/");

instances.push_back("../Group1_data_format_only_files/t60/LR1_1_DR1_3_VC1_V7a/");
instances.push_back("../Group1_data_format_only_files/t60/LR1_1_DR1_4_VC3_V11a/");
instances.push_back("../Group1_data_format_only_files/t60/LR1_1_DR1_4_VC3_V12a/");
instances.push_back("../Group1_data_format_only_files/t60/LR1_1_DR1_4_VC3_V12b/");
instances.push_back("../Group1_data_format_only_files/t60/LR1_1_DR1_4_VC3_V8a/");
instances.push_back("../Group1_data_format_only_files/t60/LR1_1_DR1_4_VC3_V9a/");
instances.push_back("../Group1_data_format_only_files/t60/LR1_2_DR1_3_VC2_V6a/");
instances.push_back("../Group1_data_format_only_files/t60/LR1_2_DR1_3_VC3_V8a/");
instances.push_back("../Group1_data_format_only_files/t60/LR2_11_DR2_22_VC3_V6a/");
instances.push_back("../Group1_data_format_only_files/t60/LR2_11_DR2_33_VC4_V11a/");
instances.push_back("../Group1_data_format_only_files/t60/LR2_11_DR2_33_VC5_V12a/");
instances.push_back("../Group1_data_format_only_files/t60/LR2_22_DR2_22_VC3_V10a/");
instances.push_back("../Group1_data_format_only_files/t60/LR2_22_DR3_333_VC4_V14a/");
instances.push_back("../Group1_data_format_only_files/t60/LR2_22_DR3_333_VC4_V17a/");

while ((opt = getopt(argc,argv,"v:p:t:s:q:c:d:k:z:b:n:g:f:o:e:r:l:m:i:h:u:")) != EOF)
	switch(opt)
	{
		case 'v': 
			modelId = atoi(optarg);
			break;
		case 'p':
			file << optarg;
			break;
		case 't': 
			 timeLimit = stod(optarg);
			break;		
		case 's':
			fixOptStr = optarg;
			break;
		case 'q':
			validIneq = atoi(optarg);
			break;			
		case 'c':
			addConstr = atoi(optarg);
			break;
		case 'd':
			tightenInvConstr = atoi(optarg);
			break;
		case 'k':
			proportionalAlpha = atoi(optarg);
			break;
		case 'z':
			reduceArcs = atoi(optarg);
			break;
		case 'b':
			tightenFlow = atoi(optarg);
			break;
		case 'n':
			nIntervals = stod(optarg);
			break;		
		case 'g':
			gapFirst = stod(optarg);
			break;
		case 'f':
			f = atoi(optarg);
			break;	
		case 'o':
			overLap = stod(optarg);
			break;
		case 'e':
			endBlock = atoi(optarg);
			break;
		case 'r':
			outVessel = atoi(optarg);
			break;
		case 'l':
			timePerIterFirst = atoi(optarg);
			break;
		case 'm':
			mIntervals = stod(optarg);
			break;
		case 'i':
			timePerIterSecond = atoi(optarg);
			break;
		case 'h':
			gapSecond = stod(optarg);
			break;
		case 'u':
			overlap2 = stod(optarg);
			break;
		case '?': 
			fprintf(stderr, ss.str().c_str());
			break;
		default: cout<<endl; abort();
	}

 switch (modelId){
	case 1: 
		#ifdef NTestInstanceSet
		mirp::milp(file.str(), timeLimit, fixOptStr);
		break;
		#endif
		#ifndef NTestInstanceSet
        for(int i=0;i<instances.size();i++){
			mirp::milp(instances[i], timeLimit, fixOptStr);
		}
		#endif
				
	case 2: 
		//~ if (endBlock + 2 > ceil(nIntervals)) {
			//~ cout << "Error! Size of end block must be at leat 2 units less than the number of intervals)" << endl;
		//~ }
		//Header
		cout << "Instance\t"  <<
				"Iter(RF)\t"  <<
				"Time(RF)\t"  <<
				"Obj(RF)\t"   <<
				"Time(FO)\t"  <<
				"Obj(FO)\t"   << 
				"CPX time\t"  << 
				"Otr time \t" << 
				"Impr \%\t"   << 
				"Inf\t"		  <<				
				"Stop by gap\t" <<
				"Stop by time\t" <<
				endl;
        #ifndef NTestInstanceSet
        for(int i=0;i<instances.size();i++){
            mirp::fixAndRelax(instances[i], fixOptStr, nIntervals, gapFirst, f, overLap, endBlock, timePerIterFirst, 
				mIntervals, timePerIterSecond, gapSecond, overlap2, timeLimit, validIneq, addConstr, tightenInvConstr,
				proportionalAlpha, reduceArcs,tightenFlow);
        }
        #endif
        #ifdef NTestInstanceSet
        mirp::fixAndRelax(file.str(), fixOptStr, nIntervals, gapFirst, f, overLap, endBlock, timePerIterFirst, 
			mIntervals, timePerIterSecond, gapSecond, overlap2, timeLimit, validIneq, addConstr, tightenInvConstr,
			proportionalAlpha, reduceArcs,tightenFlow);
        #endif
		break;
	case 3:
		//~ mirp::fixAndRelaxH(file.str(), fixOptStr, gapFirst, outVessel, timePerIterFirst, mIntervals, timePerIterSecond, gapSecond,overlap2);	
		break;
	case 4:				
		//~ mirp::fixAndRelaxHV(file.str(), nIntervals, f , overLap, endBlock, gapFirst, outVessel, timePerIterFirst, mIntervals, timePerIterSecond, gapSecond);
		break;
	default: 
		cout << "Error, no Model with index " << modelId << " was not found!"<< endl;
		break;
 } 
 
return 0;
}
