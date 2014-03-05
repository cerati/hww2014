CC = g++
ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
INCLUDE = -I../ -I./ $(ROOFITINCLUDE)
CFLAGS = -Wall -g -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB
LINKER = g++

LINKERFLAGS = $(shell root-config --ldflags)
ifeq ($(shell root-config --platform),macosx)
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif

SOURCES = SmurfLooper.cc SmurfTable.cc SmurfTableWWXSec.cc SmurfScaleFactors.cc core/SmurfSample.cc core/TH1Keys.cc core/SmurfPlotUtilities.cc core/Selections.cc
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libSmurfLooper.so

$(LIB):	$(OBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) -o $@ 

LIBS = $(LIB) 

LinkDef_out.cxx: LinkDef.h SmurfLooper.h SmurfScaleFactors.h SmurfTable.h SmurfTableWWXSec.h core/SmurfSample.h core/TH1Keys.h core/SmurfPlotUtilities.h
	rootcint -f $@ -c $(INCLUDE)  SmurfLooper.h SmurfScaleFactors.h SmurfTable.h SmurfTableWWXSec.h core/SmurfSample.h core/TH1Keys.h core/SmurfPlotUtilities.h $<

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

.PHONY: clean all ww
all:  $(LIBS)

ww:  $(LIBS) 
	root -l -q doAllHWW.C

clean:  
	rm -f *.d \
	rm -f *.o \
	rm -f *.so \
	rm -f *.cxx \
	rm -f LinkDef_out.h \
	rm -f core/*.o \
	rm -f core/*.d 

-include $(SOURCES:.cc=.d)
-include $(LIBDIR)/LinkDef_out.d

