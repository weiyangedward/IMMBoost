CC = g++ 
CCFLAGS = -Iinclude -g -Wno-deprecated -O3
OBJ = .obj
SRC = src
BIN = bin
PREFIX = ../..

install: $(BIN)/imm_build $(BIN)/imm_score


#cp $(BIN)/* $(PREFIX)/bin

all: $(BIN)/imm_build $(BIN)/imm_score
	chmod 755 $(BIN)/*

$(OBJ)/%.o: $(SRC)/%.cpp
	$(CC) $(CCFLAGS) -c $< -o $@  

$(BIN)/imm_score : $(OBJ)/imm_score.o ${OBJ}/sequence.o ${OBJ}/tree_node.o ${OBJ}/tree.o
	$(CXX) $^ -o $@

$(BIN)/imm_build : $(OBJ)/imm_build.o ${OBJ}/sequence.o ${OBJ}/tree_node.o ${OBJ}/tree.o
	$(CXX) $^ -o $@

clean:
	rm -rf ${OBJ}/*.o ${BIN}/*

