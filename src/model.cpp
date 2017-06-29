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
#define UselessVariables
//~ #define NSpotMarket
//~ #define NTravelAtCapacity
//~ #define NTravelEmpty
//~ #define NBerthLimit
//~ #define NBranching
#define NValidInequalities
//~ #define NOperateAndGo
#define NOperateAndGo2
#define NBetas
#define NNoRevisits
#define NAddAtLeftBlock

ILOSTLBEGIN

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloNumArray> NumNumMatrix;

using namespace std;
using mirp::Model;
using mirp::Instance;

#define zero -0.0001

