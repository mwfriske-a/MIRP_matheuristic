cmake_minimum_required(VERSION 2.7)
project(csp)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#set(BOOST_DIR /usr/local/boost_1_54_0/)
#set(Boost_INCLUDE_DIR /usr/local/boost_1_54_0)
#set(Boost_LIBRARY_DIR /usr/local/boost_1_54_0/libs)
#~ FIND_PACKAGE(Boost 1.54.0 REQUIRED)
FIND_PACKAGE(cplex REQUIRED)

#find_package(Boost COMPONENTS program_options timer system chrono)
include_directories(include ${Boost_INCLUDE_DIR} ${CPLEX_INCLUDE_DIR} ${CONCERT_INCLUDE_DIR})
#link_directories(${Boost_LIBRARY_DIR} ${CPLEX_LIBRARIES})

# set(cplex_DIR /opt/cplex)

set(CMAKE_CXX_FLAGS "-Wno-unknown-pragmas -g -pedantic -lm -std=c++0x -O3 -fPIC -DIL_STD -fexceptions -m64")

set(SRC src/mirp_lp.cpp src/util.cpp src/main.cpp src/fix_and_relax.cpp)# src/fix_and_relaxH.cpp  src/fix_and_relaxHV.cpp 
ADD_EXECUTABLE(mirp.out ${SRC})
TARGET_LINK_LIBRARIES(mirp.out ${CPLEX_LIBRARIES} ${Boost_LIBRARIES})
#message(${BOOST_library_dirs})

# tests
#add_executable(test.out test/main.cpp ${SRC})
#TARGET_LINK_LIBRARIES(test.out ${CPLEX_LIBRARIES} ${Boost_LIBRARIES})

#enable_testing()

#add_test(backtrackXcplex test.out ../data/simple/sub1test.gtt ../data/simple/vir3test.gtt)
