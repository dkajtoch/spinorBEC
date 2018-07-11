include Make.inc

all: spinorlib

spinorlib:
	mkdir include
	mkdir lib
	make -C $(SRCdir)

clean: 
	( cd $(SRCdir); make clean )
	( rm -r ./include ./lib )


