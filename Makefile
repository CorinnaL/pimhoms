
SRC = src
INCL = include
OBJ_PATH = obj
BIN_PATH = bin
DEBUG_BIN_PATH = bin/Debug

# defines search paths for certain types of files. syntax: vpath pattern path path path....
vpath %.o $(OBJ_PATH) 
vpath %.cpp $(SRC) 
vpath %.hpp $(INCL) 

OBJ = CLinComb.o CSymm.o CTensor.o functions.o main.o CWreathModuleElement.o calculationstex.o calculations.o homomorphism.o checkhoms.o main.o

COMPILER_FLAGS = -Wall -std=c++11 -I$(INCL)

all:  $(BIN_PATH)/out
	@echo "success!"

debug: COMPILER_FLAGS += -g
debug: LINKER_FLAGS += -g
debug: $(DEBUG_BIN_PATH)/out
	@echo "success!"

.PHONY: all debug clean

# $@ = target $^ = all prerequisites
bin/out: $(OBJ) | $(BIN_PATH)
	g++ -o $@ $^ $(LINKER_FLAGS)

bin/Debug/out: $(OBJ) | $(DEBUG_BIN_PATH) 
	g++ -o $@ $^ $(LINKER_FLAGS)

# $< = first prerequisite
$(OBJ_PATH)/main.o: main.cpp  | $(OBJ_PATH)
	g++ $(COMPILER_FLAGS) -c $< -o $@ 

$(OBJ_PATH)/calculations.o: calculations.cpp  | $(OBJ_PATH)
	g++ $(COMPILER_FLAGS) -c $< -o $@ 

$(OBJ_PATH)/calculationstex.o: calculationstex.cpp  | $(OBJ_PATH)
	g++ $(COMPILER_FLAGS) -c $< -o $@ 

$(OBJ_PATH)/CWreathModuleElement.o: CWreathModuleElement.cpp CWreathModuleElement.hpp | $(OBJ_PATH)
	g++ $(COMPILER_FLAGS) -c $< -o $@ 

$(OBJ_PATH)/CLinComb.o: CLinComb.cpp CLinComb.hpp | $(OBJ_PATH)
	g++ $(COMPILER_FLAGS) -c $< -o $@ 

$(OBJ_PATH)/CSymm.o: CSymm.cpp CSymm.hpp | $(OBJ_PATH)
	g++ $(COMPILER_FLAGS) -c $< -o $@ 

$(OBJ_PATH)/CTensor.o: CTensor.cpp CTensor.hpp | $(OBJ_PATH)
	g++ $(COMPILER_FLAGS) -c $< -o $@ 

$(OBJ_PATH)/functions.o: functions.cpp  | $(OBJ_PATH)
	g++ $(COMPILER_FLAGS) -c $< -o $@ 

$(OBJ_PATH)/homomorphism.o: homomorphisms.cpp homomorphism.hpp | $(OBJ_PATH)
	g++ $(COMPILER_FLAGS) -c $< -o $@ 

$(OBJ_PATH)/checkhoms.o: checkhoms.cpp | $(OBJ_PATH)
	g++ $(COMPILER_FLAGS) -c $< -o $@ 

$(OBJ_PATH):
	mkdir $(OBJ_PATH)
$(BIN_PATH):
	mkdir $(BIN_PATH)
$(DEBUG_BIN_PATH): |$(BIN_PATH)
	mkdir $(DEBUG_BIN_PATH)

clean:
	rm obj/*.o 
