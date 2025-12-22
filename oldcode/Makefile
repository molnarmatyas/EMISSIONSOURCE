ObjSuf        = o
SrcSuf        = cc
ExeSuf        = exe


ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit2
ROOTGLIBS    := $(shell root-config --glibs) -lMinuit2

CXX           = g++
CFLAGS        = -O3 -Wall -fPIC -fno-inline
RFLAGS        = -O3 -Wall -fPIC -fno-inline $(ROOTCFLAGS)
LD            = g++
LDFLAGS       = -O3 #-m32
RLIBS         = $(ROOTLIBS) $(SYSLIBS)
RGLIBS        = $(ROOTGLIBS) $(SYSLIBS)

SOURCES       = $(shell ls *.cc)
PROGRAMS      = $(basename $(SOURCES))

define MAKE_TEMPLATE
-include $(1).d
$(1).$(ObjSuf):
	@echo ""
	$(CXX) $(RFLAGS) -c -MD -MP -MF $(1).d $(1).$(SrcSuf) -o $$@

$(1).$(ExeSuf): $(1).$(ObjSuf)
	@echo ""
	$(LD) $(LDFLAGS) $$^ $(RGLIBS) -o $$@
endef
$(foreach program, $(PROGRAMS), $(eval $(call MAKE_TEMPLATE,$(program))))

all: $(addsuffix .exe, $(PROGRAMS))

clean:
	@rm -f *.$(ExeSuf) *.$(ObjSuf)