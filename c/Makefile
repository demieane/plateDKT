# Name your final target, that is, your executable:
TARGET_EXEC := mainDKT
CC = gcc

# List your source and build output directories:
SRC_DIR := src
OBJ_DIR := obj

# Find all the C and C++ files we want to compile
# * searches your filesystem for matching filenames
SRC := $(wildcard $(SRC_DIR)/*.c)

# From the source files, list the object files:
OBJ := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC))

.PHONY: clean 

# Now let's handle the flags
# If you want to print every warning so that you can examine 
# your code and make it better than ever, use the -Wall flag.
CFLAGS   := -Wall              # some warnings about bad code
LDLIBS   := -lm -lblas -llapack                # Left empty if no libs are needed

# this tells the make file which targets to execute based on their names
all: mydir printInfo file1 file2 $(TARGET_EXEC) 

#make obj directory only if it doesnt exist
mydir: 
	mkdir -p obj 

printInfo: 
	echo $(SRC)
	echo $(OBJ)

# major compilation into .o object file
file1: $(SRC_DIR)/mainDKT.c
	$(CC) $(CFLAGS) -c -o $(OBJ_DIR)/mainDKT.o $(SRC_DIR)/mainDKT.c $(LDLIBS)

file2: $(SRC_DIR)/funcBLAS.c
	$(CC) $(CFLAGS) -c -o $(OBJ_DIR)/funcBLAS.o $(SRC_DIR)/funcBLAS.c $(LDLIBS)

$(TARGET_EXEC): $(OBJ) 
	$(CC) -o $@ $^ $(LDLIBS)

clean: 
	rm -f $(TARGET_EXEC) 
	rm -f $(OBJ)

# The file name of the target of the rule: #@
# The names of all the prerequisites, with spaces between them. $^
# The name of the first prerequisite. #<