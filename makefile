CXX:=clang++
LFLAGS:="`mkoctfile -p LFLAGS` -O2 -std=c++11"
CXXFLAGS:="`mkoctfile -p CXXFLAGS` -O2 -std=c++11"

all: GenerateUnitSphere.oct
	
%.oct: %.cpp
	DL_LD=${CXX} CXX=${CXX} LFLAGS=${LFLAGS} CXXFLAGS=${CXXFLAGS} mkoctfile $< -o $@ -v

clean:
	rm *.o
	rm *.oct
