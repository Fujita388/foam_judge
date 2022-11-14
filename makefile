CPPFLAGS=-std=c++11
CPPFILE=main.cpp split.cpp
OFILE=main.o split.o

all: o

o:
	$(CXX) -c $(CPPFLAGS) $(CPPFILE)
	$(CXX) $(CPPFLAGS) $(OFILE)

clean:
	$(RM) *.o a.out *.lammps foam_judge.o* foam_judge.e* *.log check.dat
