
CC = g++
BIN = fdtd_omp.exe
SRC = fdtd_omp.cpp
CFLAGS = -O0 -g -DYDEBUG -Wall -fopenmp

###############################################################################
# Compilation
$(BIN) : $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(BIN)

###############################################################################
# For optimized compilation
opt:$(BIN)
	@echo "Optimization turned on."
ifneq ($(filter opt, $(MAKECMDGOALS) ),)
    CFLAGS      := -O3 -DNDEBUG -fopenmp
endif

###############################################################################
# Clean
.PHONY:c
c:clean
clean:
	rm $(BIN)
	
###############################################################################