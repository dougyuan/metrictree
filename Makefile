NUM_PROC = 16
CCFLAGS = -O3 -fopenmp -Wall -o
LIBS = -lm -lpthread
RFLAGS = -n $(NUM_PROC)
GCC  = gcc

TARGET = serial pthreads openmp

all: s p o

s: 
	$(GCC) $(CCFLAGS) serial metrictree_serial.c $(LIBS)

p:
	$(GCC) $(CCFLAGS) pthreads metrictree_pthreads.c $(LIBS)

o:
	$(GCC) $(CCFLAGS) openmp metrictree_openmp.c $(LIBS)

clean:
	rm  -rf *.o $(TARGET)
