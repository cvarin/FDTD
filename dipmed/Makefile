################################################################################
BIN  = dipmed1D.exe
BUILD_DIR = ./build
DATA_DIR  = ./output
SRC_DIR   = ./src
SRC_EXT = .cpp
HEAD_EXT = .hpp

################### Fetch names ################################################
VPATH = $(SRC_DIR)
SRC   = $(foreach DIR, $(SRC_DIR), $(wildcard $(DIR)/*$(SRC_EXT)))
HEADS = $(foreach DIR, $(SRC_DIR), $(wildcard $(DIR)/*$(HEAD_EXT)))
INCLUDES = $(addprefix -I./, $(SRC_DIR) .)
xnames = $(notdir $(SRC))
NAMES = $(subst $(SRC_EXT),,$(xnames))
OBJS  = $(addsuffix .o, $(NAMES))
LINKS = $(addprefix $(BUILD_DIR)/,$(OBJS))
DATA  = $(wildcard $(DATA_DIR)/*)

################### Commands ###################################################
CPP   = g++
FLAGS = -fPIC -O3 -Wall
# FLAGS = -O3
RM    = rm -fr

################### Compiling and linking ######################################
all: $(BIN)

$(BIN): $(LINKS) $(HEADS)
	@echo "LINKING: $(CPP) $(OBJS) -o $(BIN)"
	@$(CPP) $(LINKS) $(INCLUDES) -o $(BIN)

$(BUILD_DIR)/%.o: %$(SRC_EXT) $(HEADS)
	@echo "COMPILING: $(CPP) $(FLAGS) -c $< -o $@"
	@$(CPP) $(FLAGS) $(INCLUDES) -c $< -o $@

################### Run ########################################################
.PHONY: r run
r: run
run: $(BIN)
	./$(BIN)

################### Clean ######################################################
.PHONY: clean c
c: clean
clean:
	${RM} $(LINKS) $(BIN)

.PHONY: clean_data cd
cd: clean_data
clean_data:
	$(RM) $(DATA)

################### END of FILE ################################################

