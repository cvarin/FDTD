

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
#
# F90 compiler specification
#
F90 = f90
#
# F90 optimization, only one FLAGS line should be uncommented
#
#F90FLAGS = -cpp
#F90FLAGS = -cpp -O
#F90FLAGS = -cpp -O2
#F90FLAGS = -cpp -O3
#F90FLAGS = -cpp -n32
#F90FLAGS = -cpp -O -n32 
#F90FLAGS = -cpp -O2 -n32
#F90FLAGS = -cpp -O3 -n32
F90FLAGS = -cpp -g


all:
	@echo "To make an executable, type one of the following:"
	@echo " "
	@echo "	make ToyFDTD1c"
	@echo "	make ToyFDTD1f90"
	@echo " "
	@echo "To test one of the executables, type one of the following:"
	@echo " "
	@echo "	make testc"
	@echo "	make testf90"
	@echo " "


#
# Creation of ToyFDTD1 from C source code
#
ToyFDTD1c: ToyFDTD1c.o Makefile
	${CC} ${CFLAGS} -o ToyFDTD1c ToyFDTD1c.o -lm


ToyFDTD1c.o: ToyFDTD1.c Makefile
	${CC} ${CFLAGS} -c -o ToyFDTD1c.o ToyFDTD1.c
#
# Creation of ToyFDTD1 from F90 source code
#
ToyFDTD1f90: ToyFDTD1f90.o Makefile
	${F90} ${F90FLAGS} -o ToyFDTD1f90 ToyFDTD1f90.o -lm


ToyFDTD1f90.o: ToyFDTD1.f90 Makefile
	${F90} ${F90FLAGS} -c -o ToyFDTD1f90.o ToyFDTD1.f90




testf90: ToyFDTD1f90
	rm -f f90_runLog
	time ToyFDTD1f90 > f90_runLog


testc: ToyFDTD1c
	rm -f c_runLog
	time ToyFDTD1c > c_runLog


#
# Type: make clean
# to remove executables, core files, object files, et cetera
#
clean:
	rm -f ToyFDTD1c ToyFDTD1f90 *.o *~ core
#
# Type: make cleandata
# to remove only data files created by ToyFDTD1
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
	rm -f f90_*0.bob
	rm -f f90_*1.bob
	rm -f f90_*2.bob
	rm -f f90_*3.bob
	rm -f f90_*4.bob
	rm -f f90_*5.bob
	rm -f f90_*6.bob
	rm -f f90_*7.bob
	rm -f f90_*8.bob
	rm -f f90_*9.bob
	rm -f c_runLog f90_runLog
	rm -f ToyFDTD1c.viz ToyFDTD1f90.viz








