#====================================================================================================
# Makefile test
#====================================================================================================
bin	= em1d.exe
sources	= $(shell ls *.cpp)
headers = $(shell ls *.hpp)
objects	= $(sources:.cpp=.o)

# other variables
cc	:=  g++
cflags	:=  -Wall
lib	:=  -lm
opt	:=  -O3

#####################################################################################################
# By default, build bin.
# If any object is more recent than the program, link.
$(bin):	$(objects)
	@echo "LINKING:  $(bin)"
	@$(cc) $(cflags) $(lib) $(opt) $(objects) -o $(bin)

#####################################################################################################
# If any source file is more recent than its object file, or if any header has changed,
# it is compiled
%.o: %.cpp $(headers) 
	@echo "COMPILING: $(cc) $(cflags) $(opt) $<"
	@$(cc) -c $(cflags) $(opt) $< -o $@

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
	@echo "Removing object files"
	@( /bin/rm -f *.o )
	
.PHONY: cleanbin cb
cb:cleanbin
cleanbin:
	@echo "Removing executables"
	@( /bin/rm -f *.exe)
	
.PHONY: cleandat cd
cd:cleandat
cleandat:
	@echo "Removing data files"
	@( /bin/rm -f *.dat)
	
.PHONY: cleanall ca
ca:cleanall
cleanall: clean cleanbin cleandat

####################################################################################################
.PHONY: doc
doc: $(sources)
	@echo "Building documentation"
	doxygen doc/Doxyfile

#####################################################################################################
