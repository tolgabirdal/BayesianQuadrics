CXX:=clang++

LFLAGS:="`mkoctfile -p LFLAGS` -O2 -std=c++11"
CXXFLAGS:="`mkoctfile -p CXXFLAGS` -O2 -std=c++11"

INCLUDES:=/usr/local/Cellar/eigen/3.2.6/include/eigen3/

all: Quadric.oct GenerateUnitSphere.oct
	
GenerateUnitSphere.oct: GenerateUnitSphere.cpp
	DL_LD=${CXX} CXX=${CXX} LFLAGS=${LFLAGS} CXXFLAGS=${CXXFLAGS} mkoctfile $< -o $@ -v
	
Quadric.oct : QuadricOct.cpp Quadric.cpp Matrix.cpp 
	DL_LD=${CXX} CXX=${CXX} LFLAGS=${LFLAGS} CXXFLAGS=${CXXFLAGS} mkoctfile -I${INCLUDES} $^ -o $@ -v
	
mex: QuadricOct.cpp Quadric.cpp Matrix.cpp 
	DL_LD=${CXX} CXX=${CXX} LFLAGS=${LFLAGS} CXXFLAGS=${CXXFLAGS} mkoctfile --mex -I${INCLUDES} $^ -o Quadric.mex -v

clean:
	-rm *.o
	-rm *.oct
	-rm *.mex
