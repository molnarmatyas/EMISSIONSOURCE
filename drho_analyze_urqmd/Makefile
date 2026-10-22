#all codes ending on .cc will be attempted to compile with this Makefile
PROGRAMS = $(basename $(shell ls *.cc))

#At zeroth order, the user does not have to below this line.
WrkDir := $(shell pwd)
SrcDir = .
ExeDir = .
ObjDir = .
DepDir = .

CXX     = g++
LD      = g++
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit2
CFLAGS  = -O3 -Wall -fPIC -fno-inline $(ROOTCFLAGS) -I$(WrkDir)/$(RdrDir)
LDFLAGS = -O3

ifdef FFT
CPPFLAGS += -DEMISSIONSOURCE_USE_FFT=$(FFT)
endif
# Rebuild when switching FFT=0/1 on the command line as well as on header edits.
.PHONY: FORCE
pairsource_urqmd.o: FORCE
FORCE:

SOURCES = $(addprefix $(SrcDir)/,$(addsuffix .cc,$(PROGRAMS)))
ALL_SOURCES = $(sort $(SOURCES))

define COMPILE_TEMPLATE
-include $(DepDir)/$(notdir $(basename $(1))).d
$(ObjDir)/$(notdir $(basename $(1))).o:
	@echo ""
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c -MD -MP -MF $(DepDir)/$(notdir $(basename $(1))).d $(1) -o $$@
endef
$(foreach source, $(ALL_SOURCES), $(eval $(call COMPILE_TEMPLATE,$(source))))

define FINAL_LINK_TEMPLATE
$(ExeDir)/$(1).exe: $(2) $(ObjDir)/$(1).o
	@echo ""
	$(LD) $(LDFLAGS) $$^ $(ROOTLIBS) $(SYSLIBS) -o $$@
endef
$(foreach program, $(PROGRAMS), $(eval $(call FINAL_LINK_TEMPLATE,$(program),$(COMMON_OBJECTS))))

all: $(addprefix $(ExeDir)/,$(addsuffix .exe, $(PROGRAMS)))
	
clean:
	@rm -f $(ExeDir)/*.exe $(ObjDir)/*.o $(DepDir)/*.d
