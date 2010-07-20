

#
# Compiler specification
#
CC = cc
#
# C optimization, only one CFLAGS line should be uncommented
# The -Xcpluscomm flags is required for the SGI cc compiler to recognize the C++
# style comments ( // ) used in the code.
#
#CFLAGS = -Xcpluscomm 
#CFLAGS = -Xcpluscomm -O
#CFLAGS = -Xcpluscomm -O2
#CFLAGS = -Xcpluscomm -O3
#CFLAGS = -Xcpluscomm -n32
#CFLAGS = -Xcpluscomm -O -n32 
#CFLAGS = -Xcpluscomm -O2 -n32
#CFLAGS = -Xcpluscomm -O3 -n32
CFLAGS = -Xcpluscomm -g

#
# By default, create ToyBoxFDTDbezhig.
#
ToyBoxFDTDbezhig: ToyBoxFDTDbezhig.o Makefile
	${CC} ${CFLAGS} -o ToyBoxFDTDbezhig ToyBoxFDTDbezhig.o -lm


ToyBoxFDTDbezhig.o: ToyBoxFDTDbezhig.c Makefile
	${CC} ${CFLAGS} -c ToyBoxFDTDbezhig.c


#
# Type: make clean
# to remove executables, core files, object files, et cetera
#
clean:
	rm -f ToyBoxFDTDbezhig *.o *~ core
#
# Type: make cleandata
# to remove only data files created by ToyBoxFDTDbezhig
#
cleandata:
	rm -f *.bob
	rm -f runlog
	rm -f bezhigStimulus.dat
	rm -f bezhigCenter.dat
	rm -f ToyBoxFDTDbezhig.viz










