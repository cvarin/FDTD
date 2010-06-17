# SET THE C++ COMPILER
# CC = icc
CC = g++
#CC = CC
#CC = mpicxx

# SET THE C++ COMPILER FLAGS
CFLAGS = -c -O3

# OBJECT FILES FOR THIS PROGRAM
OBJS = jmath_basic.o jmath_geometry.o jmath_vector.o jmath_matrix.o jmath_special.o jmath_mesh.o jconsts.o jmaterials.o jplanewave.o jio.o

all: jscience

jscience: jmath_basic.o jmath_geometry.o jmath_vector.o jmath_matrix.o jmath_special.o jmath_mesh.o jconsts.o jmaterials.o jplanewave.o jio.o

jmath_basic.o: jmath_basic.h jmath_basic.cpp
	$(CC) $(CFLAGS) jmath_basic.cpp

jmath_geometry.o: jmath_geometry.h jmath_geometry.cpp
	$(CC) $(CFLAGS) jmath_geometry.cpp

jmath_vector.o: jmath_vector.h jmath_vector.cpp
	$(CC) $(CFLAGS) jmath_vector.cpp

jmath_matrix.o: jmath_matrix.h jmath_matrix.cpp
	$(CC) $(CFLAGS) jmath_matrix.cpp

jmath_special.o: jmath_special.h jmath_special.cpp
	$(CC) $(CFLAGS) jmath_special.cpp

jmath_mesh.o: jmath_mesh.h jmath_mesh.cpp
	$(CC) $(CFLAGS) jmath_mesh.cpp

jconsts.o: jconsts.h jconsts.cpp 
	$(CC) $(CFLAGS) jconsts.cpp

jmaterials.o: jmaterials.h jmaterials.cpp
	$(CC) $(CFLAGS) jmaterials.cpp

jplanewave.o: jplanewave.h jplanewave.cpp
	$(CC) $(CFLAGS) jplanewave.cpp

jio.o: jio.h jio.cpp
	$(CC) $(CFLAGS) jio.cpp

clean:
	rm *.o *~
