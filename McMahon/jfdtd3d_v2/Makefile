# SET THE C++ COMPILER
#CC = mpic++
# CC = CC
CC = mpiCC
#CC = mpicxx
#CCJazz = mpiCC

# SET THE C++ LINKER
# i. this may be the same as the C++ compiler
#CCL = mpic++
# CCL = CC
CCL = mpiCC
#CCL = mpicxx
#CCJazz = mpiCC

# SET THE C++ COMPILER FLAGS
CFLAGS = -c -O3 -DMPICH_IGNORE_CXX_SEEK

# SET THE PATH TO THE JSCIENCE LIBRARY
JSCIPATH=../jscience/

# SOURCE CODE PATH FOR JFEM3D
SRCPATH=./src/

# OUTPUT PATH
OUTPATH=./output/

# OBJECT FILES FOR THIS PROGRAM
OBJS = $(SRCPATH)jfdtd3d.o $(SRCPATH)global.o $(SRCPATH)structure.o $(JSCIPATH)jio.o $(JSCIPATH)jmaterials.o $(JSCIPATH)jconsts.o $(JSCIPATH)jmath_mesh.o $(JSCIPATH)jmath_vector.o 

all: jfdtd3d filecombine 3dto2d

jfdtd3d: jfdtd3d.o global.o structure.o
	$(CCL) $(OBJS) -o jfdtd3d

jfdtd3d.o: 
	$(CC) -I$(JSCIPATH) -I$(SRCPATH) $(CFLAGS) $(SRCPATH)jfdtd3d.cpp
	mv jfdtd3d.o $(SRCPATH)jfdtd3d.o

global.o: 
	$(CC) -I$(JSCIPATH) -I$(SRCPATH) $(CFLAGS) $(SRCPATH)global.cpp
	mv global.o $(SRCPATH)global.o

structure.o: 
	$(CC) -I$(JSCIPATH) $(CFLAGS) structure.cpp
	mv structure.o $(SRCPATH)structure.o

filecombine: filecombine.o
	$(CCL) $(SRCPATH)filecombine.o -o filecombine

filecombine.o:
	$(CC) $(CFLAGS) $(SRCPATH)filecombine.cpp
	mv filecombine.o $(SRCPATH)filecombine.o

3dto2d: 3dto2d.o
	$(CCL) $(SRCPATH)3dto2d.o -o 3dto2d

3dto2d.o:
	$(CC) $(CFLAGS) $(SRCPATH)3dto2d.cpp
	mv 3dto2d.o $(SRCPATH)3dto2d.o

clean:
	\rm $(SRCPATH)*.o *~ jfdtd3d filecombine 3dto2d

clean_results:
	\rm $(OUTPATH)*
