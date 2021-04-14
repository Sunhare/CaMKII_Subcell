DEF= -DCPU -D___PERIODIC -D___SIGMOID -D___NCX

CXX=g++
CXXFLAGS=-O3 -fopenmp -std=c++11
# CXXFLAGS=-O3


SRC_DIR := ./lib
OBJ_DIR := ./obj
SRC_FILES := $(wildcard $(SRC_DIR)/*.cc)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cc,$(OBJ_DIR)/%.o,$(SRC_FILES))

uc: uc.cc $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(DEF) -o $@ $^ 

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	$(CXX) $(CXXFLAGS) $(DEF) -c -o $@ $< 


clean: 
	rm $(OBJ_FILES) uc
	mv slurm*.out ./slurms