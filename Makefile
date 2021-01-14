DEF= -DCPU -D___PERIODIC -D___SIGMOID -D___NCX -D___KOSRCA

CXX=g++
CXXFLAGS=-O3 -fopenmp


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