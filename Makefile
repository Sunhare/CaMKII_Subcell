DEF= -DCPU -D___PERIODIC -D___SIGMOID -D___NCX

CXX=g++
CXXFLAGS=-O3 -fopenmp


uc: uc.o cell.o ap.o recsubcell.o subcell.o log.o
	$(CXX) -o $@ $(CXXFLAGS) uc.o cell.o ap.o recsubcell.o subcell.o log.o

.cc.o:
	$(CXX) -c -MMD -MP $< $(CXXFLAGS) $(DEF)

clean:
	rm *.o uc
