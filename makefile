CPPFLAGS=-std=c++11 -O3 

all: a.out

a.out: main.cpp
	$(CXX) $(CPPFLAGS) $<

clean:
	$(RM) *.o a.out *.lammps foam_judge.o* foam_judge.e* *.log *.dat
