CXX:=clang++

LFLAGS:="`mkoctfile -p LFLAGS` -O2 -std=c++11"
CXXFLAGS:="`mkoctfile -p CXXFLAGS` -O2 -std=c++11"

INCLUDES:=libs/eigen/

octexistence:=$(shell which mkoctfil | wc -c | sed "s/ //g")

targets= bindir Quadric.oct GenerateUnitSphere.oct Quadric.a mvbinaries
ifeq ("$(octexistence)","0")
	targets=bindir Quadric.a mvbinaries err
endif

all: $(targets)

err:
	@echo ""
	@echo "Octave is not installed on this machine, or it is not on the system path."
	@echo ""
	@echo "To install using Homebrew:"
	@echo "   brew install octave"
	@echo ""
	@echo "Library files are installed in ./bin"
	@echo ""
	
GenerateUnitSphere.oct: src/GenerateUnitSphere.cpp
	DL_LD=${CXX} CXX=${CXX} LFLAGS=${LFLAGS} CXXFLAGS=${CXXFLAGS} mkoctfile $< -o $@ -v
	
Quadric.oct : src/QuadricOct.cpp src/Quadric.cpp src/Matrix.cpp 
	DL_LD=${CXX} CXX=${CXX} LFLAGS=${LFLAGS} CXXFLAGS=${CXXFLAGS} mkoctfile -I${INCLUDES} $^ -o $@ -v
	
mex: src/QuadricOct.cpp src/Quadric.cpp src/Matrix.cpp 
	DL_LD=${CXX} CXX=${CXX} LFLAGS=${LFLAGS} CXXFLAGS=${CXXFLAGS} mkoctfile --mex -I${INCLUDES} $^ -o Quadric.mex -v

Quadric.a: src/Matrix.cpp src/Quadric.cpp
	${CXX} -c src/Matrix.cpp src/Quadric.cpp -Isrc -I${INCLUDES} -std=c++11
	ar rvs $@ *.o 

bindir:
	-mkdir bin

mvbinaries:
	-mv *.oct bin
	mv *.o bin
	mv *.a bin

clean:
	-rm -rf bin
