# For compilation details see:
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#

F95ROOT=/opt/intel/compilers_and_libraries/linux/mkl/gfortran
MKLROOT=/opt/intel/compilers_and_libraries/linux/mkl

COREdir		= /home/darek/Documents/REPOSITORIES/spinor-phase-transition/core
ARPACKdir	= /home/darek/src/ARPACK
BLASdir		= ${F95ROOT}/lib/intel64
LAPACKdir	= ${F95ROOT}/lib/intel64
INCLUDEdir	= $(COREdir)/include
LIBdir		= $(COREdir)/lib
SRCdir		= $(COREdir)/src


ARPACKLIB	= $(ARPACKdir)/libarpack.a
LAPACKLIB	= $(LAPACKdir)/libmkl_lapack95_lp64.a
BLASLIB		= $(BLASdir)/libmkl_blas95_lp64.a
SPINORLIB	= $(LIBdir)/libspinor.a

ALIBS		= $(SPINORLIB) $(BLASLIB) $(LAPACKLIB) $(ARPACKLIB)

.SUFFIXES: .o .f

.f90.o: 
	@$(FC) -c $(FCFLAGS) -$<

FC		      = gfortran
#FCFLAGS		= -I${INCLUDEdir} -I${F95ROOT}/include/intel64/lp64 -m64 -I${MKLROOT}/include -ffree-line-length-none -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=nan
FCFLAGS		= -I${INCLUDEdir} -I${F95ROOT}/include/intel64/lp64 -m64 -I${MKLROOT}/include -O3 -ffree-line-length-none
LDFLAGS		= $(ALIBS) -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl

ECHO		= echo
RM      	= rm
RMFLAGS 	= -f

AR = ar 
ARFLAGS = rv
RANLIB   = ranlib
