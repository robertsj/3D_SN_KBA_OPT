#===============================================================================
# User Options
#===============================================================================

PROGRAM  = 3DKBA
CC       = g++
DEBUG    = yes
OPTIMIZE = yes
PROFILE  = no
FLOAT    = float
CHUNK    = yes

#===============================================================================
# Object Files
#===============================================================================

OBJECTS = \
  main.o \
  auxiliary_function.o \
  miniapp.o \
  mesh.o \
  eas_mod.o \
  esa_mod.o \
  sae_mod.o \
  sea_mod.o \
  ase_mod.o \
  aes_mod.o \
  eas.o \
  esa.o \
  sae.o \
  sea.o \
  ase.o \
  aes.o \

#===============================================================================
# COMPILER FLAGS
#===============================================================================

ifeq ($(CC), icpc)
	CCFLAGS += -openmp -qopt-report=5 -mavx #-no-prec-div
else
	CCFLAGS += -fopenmp
endif

ifeq ($(OPTIMIZE),yes)
  CCFLAGS += -O3
else
  CCFLAGS += -O0
endif

ifeq ($(DEBUG),yes)
  CCFLAGS += -g
endif

ifeq ($(FLOAT),float)
  CCFLAGS += -DFLOAT
endif
#===============================================================================
# Targets
#===============================================================================

$(PROGRAM): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(CCFLAGS) $(LDFLAGS) $(LIB) $(INCLUDE)

div_vs_mult:
	$(CC) div_vs_mult.cc -o $@ $(CCFLAGS) $(LDFLAGS) $(LIB) $(INCLUDE)

clean:
	@rm -f *.o $(PROGRAM)

neat:
	@rm -f *.o

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .cc .o
.PHONY: clean neat

%.o: %.cc
	$(CC) $(CCFLAGS) $(INCLUDE) -c $<

#===============================================================================
# Dependencies
#===============================================================================

main.o: auxiliary_function.o
main.o: miniapp.o 

miniapp.o: auxiliary_function.o 
miniapp.o: eas_mod.o 
miniapp.o: esa_mod.o 
miniapp.o: sae_mod.o 
miniapp.o: sea_mod.o 
miniapp.o: ase_mod.o 
miniapp.o: aes_mod.o 
miniapp.o: eas.o 
miniapp.o: esa.o 
miniapp.o: sae.o 
miniapp.o: sea.o 
miniapp.o: ase.o 
miniapp.o: aes.o 

eas.o: auxiliary_function.o
esa.o: auxiliary_function.o
aes.o: auxiliary_function.o
ase.o: auxiliary_function.o
sae.o: auxiliary_function.o
sea.o: auxiliary_function.o

eas_mod.o: auxiliary_function.o
esa_mod.o: auxiliary_function.o
aes_mod.o: auxiliary_function.o
ase_mod.o: auxiliary_function.o
sae_mod.o: auxiliary_function.o
sea_mod.o: auxiliary_function.o
