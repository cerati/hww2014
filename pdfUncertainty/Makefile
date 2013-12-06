CC = g++
ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
LHAPDF = /tas/dlevans/HWW2012/CMSSW_5_2_3/src/LHAPDF-5.8.92b/
INCLUDE = -I../ -I./ $(ROOFITINCLUDE) -I$(LHAPDF)/include
CFLAGS = -Wall -g -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB
LINKER = g++

LINKERFLAGS = $(shell root-config --ldflags)
ifeq ($(shell root-config --platform),macosx)
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif

SOURCES = MyScanChain.cc ../core/Selections.cc ../core/SmurfPlotUtilities.cc ../core/SmurfSample.cc ../SmurfScaleFactors.cc
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libSmurfPDFLooper.so

$(LIB):	$(OBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) -o $@ 

LIBS = $(LIB)

LinkDef_out.cxx: LinkDef.h MyScanChain.h ../core/SmurfPlotUtilities.h ../core/SmurfSample.h ../SmurfScaleFactors.h
	rootcint -f $@ -c $(INCLUDE)  MyScanChain.h  ../core/SmurfPlotUtilities.h ../core/SmurfSample.h ../SmurfScaleFactors.h $<

# General rule for making object files
%.d:	%.cc
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@
%.d:	%.cxx
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

%.o: 	%.cc 
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cxx
	$(CC) $(CFLAGS) $< -c -o $@

.PHONY: clean all dy
all:  $(LIBS)

pdf:  $(LIBS) 
	root -l -q doWW.C

clean:  
	rm -f *.d \
	rm -f *.o \
	rm -f *.so \
	rm -f *.cxx \
	rm -f ../core/*.d \
    rm -f ../core/*.o 


-include $(SOURCES:.cc=.d)
-include $(LIBDIR)/LinkDef_out.d

