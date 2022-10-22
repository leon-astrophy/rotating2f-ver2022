include ./make.inc

SOURCES=Definition.f90 \
	RotatingStar.f90 \
	Analysis.f90 \
	EOSTable.f90 \
	FindDensity.f90 \
	FindRotation.f90 \
	FindPotential.f90 \
	Initial.f90 \
	Openfile.f90 \
       	Output.f90 \
	SCFConsistent.f90 \
	TidalLove.f90

MODULES=Definition.mod
OBJECTS=$(SOURCES:.f90=.o )

ROTATING: Definition.o $(OBJECTS)  
	$(F90) $(LDFLAGS) -o ./ROTATING $(OBJECTS) 

$(OBJECTS): %.o: %.f90 
	$(F90) $(F90FLAGS) -c $< -o $@

ROTATING.o: Definition.o

clean:
	rm -rf Definition
	rm -rf *.o 	
	rm -rf *.mod
	rm -rf ROTATING
cleanfile:
	rm -rf ./Profile/*.dat
	rm -rf ./Parameter/*.dat
	rm -rf tmp.txt
cleanpic:
	rm -rf ./Plot/*.png
