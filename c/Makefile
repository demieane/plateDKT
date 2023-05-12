# Name your final target, that is, your executable:
# It compiles with 
# >> gcc -Wall -g3 -fsanitize=address ./src/mainDKT.cpp ./src/funcFEM.cpp ./src/funcMat.cpp ./src/funcFSI.cpp -o mainDKT_CPP -lm -lblas -llapack
TARGET_EXEC := mainDKT_CPP
CC = gcc
#CC = nvcc

# List your source and build output directories:
SRC_DIR := src
OBJ_DIR := obj

# Find all the C and C++ files we want to compile
# * searches your filesystem for matching filenames
SRC := $(wildcard $(SRC_DIR)/*.cpp)

BIN_DIR := $(wildcard *.bin)

# From the source files, list the object files:
OBJ := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC))

.PHONY: clean 

# Now let's handle the flags
# If you want to print every warning so that you can examine 
# your code and make it better than ever, use the -Wall flag.
CFLAGS   := -Wall -g3             # some warnings about bad code
#CFLAGS   :=               # some warnings about bad code
LDFLAGS  :=  -lm -lblas -llapack               # Left empty if no libs are needed

# this tells the make file which targets to execute based on their names
all: mydir printInfo file1 file2 file3 file4 $(TARGET_EXEC) 

#make obj directory only if it doesnt exist
mydir: 
	mkdir -p obj 

printInfo: 
	echo $(SRC)
	echo $(OBJ)
#echo $(BIN_DIR)

# major compilation into .o object file
file1: $(SRC_DIR)/mainDKT.cpp
	$(CC) $(CFLAGS) -c -o $(OBJ_DIR)/mainDKT.o $(SRC_DIR)/mainDKT.cpp $(LDFLAGS) 

file2: $(SRC_DIR)/funcFEM.cpp
	$(CC) $(CFLAGS) -c -o $(OBJ_DIR)/funcFEM.o $(SRC_DIR)/funcFEM.cpp $(LDFLAGS) 

file3: $(SRC_DIR)/funcFSI.cpp
	$(CC) $(CFLAGS) -c -o $(OBJ_DIR)/funcFSI.o $(SRC_DIR)/funcFSI.cpp $(LDFLAGS) 

file4: $(SRC_DIR)/funcMat.cpp
	$(CC) $(CFLAGS) -c -o $(OBJ_DIR)/funcMat.o $(SRC_DIR)/funcMat.cpp $(LDFLAGS) 


$(TARGET_EXEC): $(OBJ) 
	$(CC) -o $@ $^ $(LDFLAGS)

clean: 
	rm -f $(TARGET_EXEC) 
	rm -f $(OBJ)
#	rm -f $(BIN_DIR)

# The file name of the target of the rule: #@
# The names of all the prerequisites, with spaces between them. $^
# The name of the first prerequisite. #<