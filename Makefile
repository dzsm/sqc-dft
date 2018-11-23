CC      = g++
CFLAGS  = -O3 -I/usr/include/eigen3  --std=c++11
LDFLAGS = -fopenmp  -I/usr/include/eigen3

all: sqc-dft clean

sqc-dft: dft.o compute_K_Un_S.o gausslegendre128_quadrature.o jsoncpp.o lebedev59_quadrature.o scf.o site.o site_type.o slater_poisson_maxn7.o solid_real_sphericalharmonics_maxl5.o sto_params.o diis.o
	$(CC) -o $@ $^ $(LDFLAGS)

dft.o : main.cpp
	$(CC) -o dft.o -c $(CFLAGS) $<


compute_K_Un_S.o : compute_K_Un_S.cpp computations.h
	$(CC) -c $(CFLAGS) $<


%.o: %.cpp
	$(CC) -c $(CFLAGS) $<

.PHONY: clean cleanest

clean:
	rm *.o

cleanest: clean
	rm sqc-dft


