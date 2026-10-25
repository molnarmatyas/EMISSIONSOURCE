CXX ?= g++
ROOT_CONFIG ?= root-config
ROOT_CFLAGS := $(shell $(ROOT_CONFIG) --cflags)
ROOT_LIBS := $(shell $(ROOT_CONFIG) --libs)
CXXFLAGS ?= -O3 -Wall -Wextra
CPPFLAGS += $(ROOT_CFLAGS)
FFTW_CFLAGS ?= $(shell pkg-config --cflags fftw3)
FFTW_LIBS ?= $(shell pkg-config --libs fftw3)

# With FFT unset, the common header controls the method. FFT=1 overrides it for
# this build; FORCE prevents stale binaries when switching build options.
ifdef FFT
CPPFLAGS += -DEMISSIONSOURCE_USE_FFT=$(FFT)
endif

.PHONY: all check FORCE
all: transform_fft.exe CqFit.exe calc_and_plot_syserr_FFT.exe compare_methods.exe

transform_fft.exe: transform_fft.cc fft_common.h ../pair_source_tree.h ../header_for_all_emissionsource.h FORCE
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $(FFTW_CFLAGS) $< -o $@ $(ROOT_LIBS) $(FFTW_LIBS)

CqFit.exe: CqFit.cpp fft_common.h ../header_for_all_emissionsource.h FORCE
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $< -o $@ $(ROOT_LIBS) -lMinuit2

calc_and_plot_syserr_FFT.exe: calc_and_plot_syserr_FFT.cpp fft_common.h ../header_for_all_emissionsource.h FORCE
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $< -o $@ $(ROOT_LIBS)

compare_methods.exe: compare_methods.cpp fft_common.h ../header_for_all_emissionsource.h FORCE
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $< -o $@ $(ROOT_LIBS)

check:
	bash tests/check.sh

FORCE:
