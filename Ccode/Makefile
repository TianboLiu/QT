# A Generic Makefile for compiling ROOT programs
# R. Michaels, rom@jlab.org, Aug 2001  See also README !!
# Version of this release
#
O = example
 
ROOTLIBS      = $(shell root-config --libs) -MathMore -lHistPainter
ROOTGLIBS     = $(shell root-config --glibs)
LHAPDFLIBS    = $(shell lhapdf-config --cflags --ldflags)
GSLLIBS       = -lgsl -lgslcblas -lm
CXX           = g++
CXXFLAGS      = -lgfortran -lgsl -lgslcblas -lm \
		-Wall -frtti -fexceptions -fPIC \
                   -DLINUXVERS -I$(ROOTSYS)/include -I/var/phy/project/mepg/tl190/LHAPDF/include -I/var/phy/project/mepg/tl190/LHAPDF/boost_1_59_0 -L/var/phy/project/mepg/tl190/LHAPDF/lib -lLHAPDF -O

FFLAGS	      = -lgfortran -O -W -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8

# Linux with g++
INCLUDES      = -I$(ROOTSYS)/include 
CXX           = g++
LD            = g++
LDFLAGS       = -O3

LIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS)

ALL_LIBS =  $(GLIBS) $(LIBS) $(LHAPDFLIBS) $(GSLLIBS)

# The following sources comprise the package of scaler classes by R. Michaels.
SRC = $(O).C

HEAD = $(SRC:.C=.h)
DEPS = $(SRC:.C=.d)
SCALER_OBJS = $(SRC:.C=.o)

F_OBJS =
cc_OBJS =
total_OBJS = $(F_OBJS) $(cc_OBJS)

# Test code executibles
PROGS = $(O)

%.o: %.f
	$(CXX) $(FFLAGS) -c $< -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

$(O): $(O).o $(O).C $(total_OBJS)
	rm -f $@
	$(CXX) $(CXXFLAGS) -o $@ $(O).o $(total_OBJS) $(ALL_LIBS)


clean:
	rm -f *.o core *~ *.d *.tar $(PROGS)

realclean:  clean
	rm -f *.d

###

.SUFFIXES:
.SUFFIXES: .c .cc .cpp .C .o .d

%.o:	%.C
	$(CXX) $(CXXFLAGS) -c $<

%.d:	%.C
	@echo Creating dependencies for $<
	@$(SHELL) -ec '$(CXX) -MM $(CXXFLAGS) -c $< \
		| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

-include $(DEPS)








