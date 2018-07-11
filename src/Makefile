include ../Make.inc

OBJS	= data_struct.o sparse.o hamiltonian.o density_matrix.o measures.o

all: lib

lib: $(OBJS)
	mv $(SRCdir)/*.mod $(INCLUDEdir)
	$(AR) $(ARFLAGS) $(SPINORLIB) $(OBJS)
	$(RANLIB) $(SPINORLIB)

%.o: %.f90
	$(FC) -c $(FCFLAGS) $<

clean:
	( cd $(SRCdir); $(RM) $(RMFLAGS) *.o )
	( cd $(INCLUDEdir); $(RM) $(RMFLAGS) *.mod )
	( cd $(LIBdir); $(RM) $(RMFLAGS) *.a )
