cmake_minimum_required(VERSION 2.7)
project(MIRP)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
# set(BOOST_ROOT /usr/local/boost_1_73_0/)
# set(BOOST_INCLUDEDIR /usr/local/boost_1_73_0/)
# set(BOOST_LIBRARYDIR /usr/local/boost_1_73_0/libs/)
# find_package(Boost 1.73.0 COMPONENTS program_options REQUIRED)
find_package(cplex REQUIRED)

include_directories(include ${Boost_INCLUDEDIR} ${CPLEX_INCLUDE_DIR} ${CONCERT_INCLUDE_DIR})
#link_directories(${Boost_LIBRARY_DIR} ${CPLEX_LIBRARIES})

# set(cplex_DIR /opt/cplex)

set(CMAKE_CXX_FLAGS "-Wno-unknown-pragmas -g -pedantic -lm -std=c++0x -O0 -O3 -fPIC -DIL_STD -fexceptions -m64")

set(SRC src/mirp_lp.cpp src/util.cpp src/main.cpp src/fix_and_relax.cpp)# src/fix_and_relaxH.cpp  src/fix_and_relaxHV.cpp 
ADD_EXECUTABLE(mirp.out ${SRC})
TARGET_LINK_LIBRARIES(mirp.out ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} -ldl)
#message(${BOOST_library_dirs})

# tests
#add_executable(test.out test/main.cpp ${SRC})
#TARGET_LINK_LIBRARIES(test.out ${CPLEX_LIBRARIES} ${Boost_LIBRARIES})

#enable_testing()
