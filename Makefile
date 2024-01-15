BINDIR := ../bin
SRCDIR := src
OBJDIR := obj

SRCEXT := .f90 
OBJEXT := .o

FC      := gfortran
FFLAGS  := -cpp
LDFLAGS := -cpp
LIB    :=  
.SUFFIXES: .f90 .o

#-----------------------------------------
OTHER :=  src/prec src/string src/iofile \
		src/clock \
		src/random \
		src/arrays_utility \
		src/HS_puregas \
		src/DMC_params \
		src/DMC_print \
		src/DMC \

PROGS := src/main

ALL_FILES  = $(OTHER) $(PROGS)
ALL_SOURCE = $(ALL_FILES:%=%$(SRCEXT))
ALL_OBJS   := $(subst $(SRCDIR)/,$(OBJDIR)/,$(ALL_FILES:%=%.o))



.PHONY: all
all: $(PROGS)

.PHONY: remake
remake: clean
remake: all

.PHONY: debug
debug: clean
debug: FFLAGS += -g -W -Wall -fcheck=all 
debug: LDFLAGS+= -g -W -Wall -fcheck=all
#decomment to turnoff the interaction between particles
debug: FFLAGS += #-DTURNOFF_INTERACTION
debug: LDFLAGS+= #-DTURNOFF_INTERACTION
debug: all

.PHONY:parallel
parallel: clean
parallel: FFLAGS += -fopenmp
parallel: LDFLAGS+= -fopenmp
parallel: all

.PHONY: fast
fast: clean
fast: FFLAGS += -O2
fast: all

.PHONY: profiling
profiling: clean
profiling: FFLAGS += -pg
profiling: LDFLAGS+= -pg
profiling: all


$(PROGS): $(ALL_OBJS)
	$(FC) $(LDFLAGS) $(ALL_OBJS) -o $(notdir $@).x 

$(OBJDIR)/%.o: $(SRCDIR)/%$(SRCEXT)
	$(FC) $(FFLAGS) -c -o $@ $<


.PHONY: clean
clean:
	-rm -f obj/*.o *.mod *.x