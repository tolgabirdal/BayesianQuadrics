CXX:=clang++

LFLAGS:="`mkoctfile -p LFLAGS` -O2 -std=c++11"
CXXFLAGS:="`mkoctfile -p CXXFLAGS` -O2 -std=c++11"

INCLUDES:=/usr/local/Cellar/eigen/3.2.8/include/eigen3/

all: bindir Quadric.oct GenerateUnitSphere.oct mvbinaries
	
GenerateUnitSphere.oct: src/GenerateUnitSphere.cpp
	DL_LD=${CXX} CXX=${CXX} LFLAGS=${LFLAGS} CXXFLAGS=${CXXFLAGS} mkoctfile $< -o $@ -v
	
Quadric.oct : src/QuadricOct.cpp src/Quadric.cpp src/Matrix.cpp 
	DL_LD=${CXX} CXX=${CXX} LFLAGS=${LFLAGS} CXXFLAGS=${CXXFLAGS} mkoctfile -I${INCLUDES} $^ -o $@ -v
	
mex: src/QuadricOct.cpp src/Quadric.cpp src/Matrix.cpp 
	DL_LD=${CXX} CXX=${CXX} LFLAGS=${LFLAGS} CXXFLAGS=${CXXFLAGS} mkoctfile --mex -I${INCLUDES} $^ -o Quadric.mex -v

bindir:
	-mkdir bin

mvbinaries:
	mv *.oct bin
	mv *.o bin

clean:
	-rm -rf bin
