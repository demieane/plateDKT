# Name your final target, that is, your executable:
TARGET_EXEC := final_program
CC = gcc

# List your source and build output directories:
SRC_DIR := src
OBJ_DIR := obj

# Find all the C and C++ files we want to compile
# * searches your filesystem for matching filenames
SRC := $(wildcard $(SRC_DIR)/*.c)

# From the source files, list the object files:
OBJ := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC))

.PHONY: clean #myprint

# Now let's handle the flags
# If you want to print every warning so that you can examine 
# your code and make it better than ever, use the -Wall flag.
CFLAGS   := -Wall              # some warnings about bad code
LDFLAGS  := #-Llib              # -L is a linker flag
LDLIBS   := -lm -lblas -llapack                # Left empty if no libs are needed

# this tells the make file which targets to execute based on their names
all: compileCode printInfo
#compileCode greetings

printInfo:
	echo $(SRC)
	echo $(OBJ)

compileCode:
	$(CC) $(CFLAGS) -o $(TARGET_EXEC) $(SRC) $(LDLIBS)
# $(CC) -o $(TARGET_EXEC) mainFem2.c femLibrary.c $(LDLIBS)

clean: 
	rm -f $(TARGET_EXEC)