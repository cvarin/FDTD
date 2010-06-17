#====================================================================================================
# Makefile test
#====================================================================================================
cc	:=  g++
cflags	:=  -Wall
lib	:=  -lm
opt	:=  -O3

#####################################################################################################
# Build bin.
all: 1d 2d
1d:main_1d.exe
main_1d.exe: main_1d.cpp
	@echo "COMPILING: $(cc) $(cflags) $(opt) $<"
	@$(cc) $(cflags) $(opt) $< -o $@
2d:main_2d.exe
main_2d.exe: main_2d.cpp
	@echo "COMPILING: $(cc) $(cflags) $(opt) $<"
	@$(cc) $(cflags) $(opt) $< -o $@

#####################################################################################################
# Use the .PHONY declaration went you want to prevent the target ("debug" here) from being
# interpreted as the name of a file. Here, if there was a file named "debug" in the directory,
# make would not do anything because there are no dependencies (at the right of the colon), so
# the target is always up to date!
.PHONY:	debug
debug:
	make -e opt="-g"

#####################################################################################################
.PHONY:	prof
prof:
	make -e opt="-pg"

#####################################################################################################
.PHONY: clean c
c:clean
clean:
	@echo "Removing executable and data files"
	@( /bin/rm -f *.exe *.dat)

#####################################################################################################
