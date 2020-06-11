#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <boost/program_options.hpp>
#include "../include/util.h"
#include "../include/model.h"

#define NTestInstanceSet
using namespace std;
namespace po = boost::program_options;

int main(const int argc, char** argv) {
	 try {	 	
	 	string file;
		string fixOptStr;
		int modelId;
		double timeLimit, gapFirst, gapSecond, overLap, intervalsA, intervalsB, nIntervals, overlapA, overlapB;
		int opt, endBlock, outVessel, timePerIterFirst, timePerIterSecond;
		bool validIneq, addConstr, tightenInvConstr, proportionalAlpha, reduceArcs, tightenFlow;
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "Produce help message")
            ("version,V", po::value<int>(&modelId), "Relax-and-fix (R&F) version (1 - for MILP; 2- for fix-and-relax vertically; 3-fix-and-relax horizontally)")
            ("path,P", po::value<string>(&file), "Path to the instance file")
            
            ("sizeInterval,N", po::value<double>(&nIntervals), "Number of time periods per interval")
            ("endBlock,E", po::value<int>(&endBlock)->default_value(0), "Number of intervals initially at the endblock (default 0)")
            ("overlapRf,O", po::value<double>(&overLap)->default_value(0), "Overlap between intervals in the R&F iterations (default 0)")
            ("gap,G", po::value<double>(&gapFirst)->default_value(0.0001), "Initial GAP at the first iteration (default 0.0001)")
            ("timeItRf,L", po::value<int>(&timePerIterFirst), "Time limit per iteration of R&F")

            ("addConst,C", po::value<bool>(&addConstr)->default_value(0), "Use additional constraints (1) (default 0 - do not use)")
            ("validIneq,Q", po::value<bool>(&validIneq)->default_value(0), "Use valid inequalities (1) (default 0 - do not use)")
            ("tightInventory,D", po::value<bool>(&tightenInvConstr)->default_value(0), "(For endblock > 0) Use tight inventory constraints (1) (default 0 - do not use)")
            ("propAlpha", po::value<bool>(&proportionalAlpha)->default_value(0), "(For endblock > 0) Use proportional alpha limits (1) (default 0 - do not use)")
            ("reducedArc", po::value<bool>(&reduceArcs)->default_value(1), "Ignore inter-regional arcs between ports of the same type (default 1) ")
            ("tightFlow", po::value<bool>(&tightenFlow)->default_value(0), "Use thight a flow by adding new constraints (default 0 - do not use)")

			("fixOpt", po::value<string>(&fixOptStr), "Order and fix-and-optmize procedures used: A - Vessel time intervals; B - Time Intervals; C - Vessels Pairs; D - Port Regions")
			("numintervalsA", po::value<double>(&intervalsA)->default_value(1), "Number of intervals in the Fix-and-optimize A: Interval Vessel (default 1)")
			("numintervalsB", po::value<double>(&intervalsB)->default_value(1), "Number of intervals in the Fix-and-optimize B: Time intervals (default 1)")
			("overlapA", po::value<double>(&overlapA)->default_value(0), "Overlap between intervals in the Fix-and-optimize A: Interval Vessel (default 0)")
			("overlapB", po::value<double>(&overlapB)->default_value(0), "Overlap between intervals in the Fix-and-optimize B: Time intervals (default 0)")
			("timeItFo",po::value<int>(&timePerIterSecond)->default_value(0), "Time limit per Fix-and-optimize iteration (default 0)")
			("timeFo",po::value<double>(&timeLimit)->default_value(0), "Time limit per Fix-and-optimize procedure (default 0)")
			("gaoFo",po::value<double>(&gapSecond)->default_value(0.0001), "Gap considered for fix-and-optimize (default 0.0001)")

        ;
        po::variables_map vm;        
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;
        }      	    

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

	// instances.push_back("../Group1_data_format_only_files/t60/LR1_1_DR1_3_VC1_V7a/");
	// instances.push_back("../Group1_data_format_only_files/t60/LR1_1_DR1_4_VC3_V11a/");
	// instances.push_back("../Group1_data_format_only_files/t60/LR1_1_DR1_4_VC3_V12a/");
	// instances.push_back("../Group1_data_format_only_files/t60/LR1_1_DR1_4_VC3_V12b/");
	// instances.push_back("../Group1_data_format_only_files/t60/LR1_1_DR1_4_VC3_V8a/");
	// instances.push_back("../Group1_data_format_only_files/t60/LR1_1_DR1_4_VC3_V9a/");
	// instances.push_back("../Group1_data_format_only_files/t60/LR1_2_DR1_3_VC2_V6a/");
	// instances.push_back("../Group1_data_format_only_files/t60/LR1_2_DR1_3_VC3_V8a/");
	// instances.push_back("../Group1_data_format_only_files/t60/LR2_11_DR2_22_VC3_V6a/");
	// instances.push_back("../Group1_data_format_only_files/t60/LR2_11_DR2_33_VC4_V11a/");
	// instances.push_back("../Group1_data_format_only_files/t60/LR2_11_DR2_33_VC5_V12a/");
	// instances.push_back("../Group1_data_format_only_files/t60/LR2_22_DR2_22_VC3_V10a/");
	// instances.push_back("../Group1_data_format_only_files/t60/LR2_22_DR3_333_VC4_V14a/");
	// instances.push_back("../Group1_data_format_only_files/t60/LR2_22_DR3_333_VC4_V17a/");
	
	 switch (modelId){
		case 1: 
			#ifdef NTestInstanceSet
			mirp::milp(file, timeLimit, fixOptStr);			
			#endif
			#ifndef NTestInstanceSet
	        for(int i=0;i<instances.size();i++){
				mirp::milp(instances[i], timeLimit, fixOptStr);
			}			
			#endif
			break;					
		case 2: 
			//Header
			// cout << "Instance\t"  <<
					// "Iter(RF)\t"  <<
					// "Time(RF)\t"  <<
					// "Obj(RF)\t"   <<
					// "Time(FO)\t"  <<
					// "Obj(FO)\t"   << 
					// "CPX time\t"  << 
					// "Otr time \t" << 
					// "Impr \%\t"   << 
					// "Inf\t"		  <<				
					// "Stop by gap\t" <<
					// "Stop by time\t" <<
					// endl;
	        #ifndef NTestInstanceSet
	        for(int i=0;i<instances.size();i++){
	            mirp::fixAndRelax(instances[i], nIntervals, endBlock, overLap, gapFirst, timePerIterFirst, 
					validIneq, addConstr, tightenInvConstr,	proportionalAlpha, reduceArcs, tightenFlow, 
					fixOptStr, intervalsA, intervalsB, overlapA, overlapB, timePerIterSecond, timeLimit, gapSecond);
	        }	        
	        #endif
	        #ifdef NTestInstanceSet
	        mirp::fixAndRelax(file, nIntervals, endBlock, overLap, gapFirst, timePerIterFirst, 
					validIneq, addConstr, tightenInvConstr,	proportionalAlpha, reduceArcs, tightenFlow, 
					fixOptStr, intervalsA, intervalsB, overlapA, overlapB, timePerIterSecond, timeLimit, gapSecond);
	        #endif
			break;
		case 3:
			//~ mirp::fixAndRelaxH(file.str(), fixOptStr, gapFirst, outVessel, timePerIterFirst, intervalsA, timePerIterSecond, gapSecond,overlapA);	
			break;
		case 4:				
			//~ mirp::fixAndRelaxHV(file.str(), nIntervals, f , overLap, endBlock, gapFirst, outVessel, timePerIterFirst, intervalsA, timePerIterSecond, gapSecond);
			break;
		default: 
			cout << "Error, no Model with index " << modelId << " was not found!"<< endl;
			break;
 	} 
 	}	
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }
 
return 0;
}
