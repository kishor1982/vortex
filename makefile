# BEC-GP-ROT-OMP programs are developed by:
#
# R. Kishor Kumar
# (Instituto de Fisica, Universidade de Sao Paulo, Sao Paulo, Brazil)
#
# Vladimir Loncar, Antun Balaz
# (Scientific Computing Laboratory, Center for the Study of Complex Systems, Institute of Physics Belgrade, Serbia)
#
# Paulsamy Muruganandam
# (Department of Physics, Bharathidasan University, Tiruchirappalli, Tamil Nadu, India)
#
# Sadhan K. Adhikari
# (Instituto de Fisica Teorica, UNESP - Sao Paulo State University, Brazil)
#
# Public use and modification of these codes are allowed provided that the following papers are cited:
# [1] L. E. Young-S. et al., Comput. Phys. Commun. 220 (2017) 503.
# [2] P. Muruganandam and S. K. Adhikari, Comput. Phys. Commun. 180 (2009) 1888.
# [3] D. Vudragovic et al., Comput. Phys. Commun. 183 (2012) 2021.
# [4] R. Kishor Kumar et al., Comput. Phys. Commun. 195 (2015) 117.
# [5] B. Sataric et al., Comput. Phys. Commun. 200 (2016) 411.
# [6] V. Loncar et al., Comput. Phys. Commun. 200 (2016) 406.
# [7] L. E. Young-S. et al., Comput. Phys. Commun. 204 (2016) 209.
# [8] V. Loncar et al., Comput. Phys. Commun. 209 (2016) 190.
#
# The authors would be grateful for all information and/or comments
# regarding the use of the programs.

ifndef $(compiler)
    compiler = gnu
endif

ifeq ($(compiler), gnu)
   FC = gfortran
   OMPFLAGS = -O3 -fopenmp
endif

ifeq ($(compiler), intel)
   FC = ifort
   OMPFLAGS = -O3 -qopenmp -w -mcmodel medium -shared-intel -warn unused
#If older version of Intel C compiler is present, the option -qopenmp in the
#line above should be replaced by -openmp
endif

ifeq ($(compiler), pgi)
   FC = pgfortran
   OMPFLAGS = -O3 -fast -mp=allcores
#Explicit control of the number of OpenMP threads can be achieved by removing
#'=allcores' option and using OMP_NUM_THREADS environment variable
endif

ifeq ($(compiler), oracle)
   FC = sunf90
   OMPFLAGS = -fopenmp -fast
endif

all: bec-gp-rot-2d-th bec-gp-rot-3d-th

help:	README.md
	less $^

bec-gp-rot-2d-th:
	$(FC) $(OMPFLAGS) -c src/bec-gp-rot-2d-th.f90 -o bec-gp-rot-2d-th.o
	$(FC) $(OMPFLAGS) bec-gp-rot-2d-th.o -o bec-gp-rot-2d-th
	rm *.o *.mod

bec-gp-rot-3d-th:
	$(FC) $(OMPFLAGS) -c src/bec-gp-rot-3d-th.f90 -o bec-gp-rot-3d-th.o
	$(FC) $(OMPFLAGS) bec-gp-rot-3d-th.o -o bec-gp-rot-3d-th
	rm *.o *.mod

full-clean: clean # Cleans all the output files also
	rm -rf imag*.txt real*.txt

clean:
	rm -rf *.mod *~ *.o bec-gp-rot-2d-th  bec-gp-rot-3d-th
