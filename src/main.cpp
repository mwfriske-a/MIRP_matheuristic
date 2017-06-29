#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include "../include/util.h"
#include "../include/model.h"

using namespace std;

int main(const int argc, char** argv) {
 if (argc < 7) {
	cout << "Error! At least 4 parameters must be defined" << endl;	
 }
stringstream ss;
ss << "usuage is \n -v : Version (1 - for MILP; 2- for fix-and-relax vertically; 3-fix-and-relax horizontally) \n " 
	<< "-p : instance path \n " 
	<< "-t : time limit in seconds - for -v = 2, corresponds to timelimit of local searches \n " 
	<< "-s : A string to put together in the logFile name \n" 
	<< "For -v = {2,3} (Fix and relax algorithm):\n " 
	<< "-n : Number of intervals in the first phase (-v =2) \n "
	<< "-g : Initial GAP for solving the first phase (-v = 2,3) \n"
	<< "-f : Number of intervals that leaves the end block in each iteration (-v =2) \n "
	<< "-o : percentage of overlapping between the fixed and the integer interval (-v =2) \n "
	<< "-e : Size of the end block/vessels out of model (must be at most n-2) \n "
	<< "-r : Vessels out of model (must be at most v-2) \n "
	<< "-l : Time limit for solving each interval \n "
	<< "-m : Number of intervals in the improvement phase (-v =2)\n "
	<< "-i : Time limit for one loop in the improvement phase \n "
	<< "-h : Fixed GAP in the imnprovement phase \n "
	<< "-u : Overlap in the imnprovement phase \n ";
int modelId;
stringstream file;
string optstr;
double timeLimit, gapFirst, gapSecond, overLap, mIntervals, nIntervals, overlap2;
int opt, endBlock, outVessel, timePerIterFirst, timePerIterSecond, f ;
while ((opt = getopt(argc,argv,"v:p:f:t:y:s:n:g:o:e:r:l:m:i:h:u:")) != EOF)
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
			optstr = optarg;
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
		case '?': fprintf(stderr, ss.str().c_str());
		
		default: cout<<endl; abort();
	}

 switch (modelId){
	case 1: 
		mirp::milp(file.str(), timeLimit, optstr);
		break;
	case 2: 
		if (endBlock + 2 > ceil(nIntervals)) {
			cout << "Error! Size of end block must be at leat 2 units less than the number of intervals)" << endl;
		}
		//~ mirp::fixAndRelax(file.str(), optstr, nIntervals, gapFirst, f, overLap, endBlock, timePerIterFirst, mIntervals, timePerIterSecond, gapSecond, overlap2, timeLimit);
		break;
	case 3:
		//~ mirp::fixAndRelaxH(file.str(), optstr, gapFirst, outVessel, timePerIterFirst, mIntervals, timePerIterSecond, gapSecond,overlap2);	
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
