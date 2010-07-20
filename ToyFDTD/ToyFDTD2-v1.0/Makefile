

#
# C compiler specification
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

all:
	@echo "To make an executable, type the following:"
	@echo " "
	@echo "	make ToyFDTD2c"
	@echo " "
	@echo "To test the executable, type the following:"
	@echo " "
	@echo "	make testc"
	@echo " "


#
# Creation of ToyFDTD2 from C source code
#
ToyFDTD2c: ToyFDTD2c.o Makefile
	${CC} ${CFLAGS} -o ToyFDTD2c ToyFDTD2c.o -lm


ToyFDTD2c.o: ToyFDTD2.c Makefile
	${CC} ${CFLAGS} -c -o ToyFDTD2c.o ToyFDTD2.c


testc: ToyFDTD2c
	rm -f c_runLog
	time ToyFDTD2c > c_runLog


#
# Type: make clean
# to remove executables, core files, object files, et cetera
#
clean:
	rm -f ToyFDTD2c *.o *~ core
#
# Type: make cleandata
# to remove only data files created by ToyFDTD2
#
cleandata:
	rm -f c_*0.bob
	rm -f c_*1.bob
	rm -f c_*2.bob
	rm -f c_*3.bob
	rm -f c_*4.bob
	rm -f c_*5.bob
	rm -f c_*6.bob
	rm -f c_*7.bob
	rm -f c_*8.bob
	rm -f c_*9.bob
	rm -f c_runLog
	rm -f ToyFDTD2c.viz








