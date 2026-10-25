#!/bin/bash
set -e
cd "$(dirname "$0")/.."

# These small test programs stay separate from the analysis/fit code.
make FFT=1 all
root_flags="$(root-config --cflags --libs)"
fftw_flags="$(pkg-config --cflags --libs fftw3)"

g++ -O2 -std=c++17 -DEMISSIONSOURCE_USE_FFT=1 ../drho_analyze_urqmd/pairsource_urqmd.cc -o /tmp/emissionsource_pairsource_fft.exe $root_flags
g++ -O2 -std=c++17 -DEMISSIONSOURCE_USE_FFT=0 ../drho_analyze_urqmd/pairsource_urqmd.cc -o /tmp/emissionsource_pairsource_legacy.exe $root_flags
g++ -O2 -std=c++17 -DEMISSIONSOURCE_USE_FFT=0 transform_fft.cc -o /tmp/emissionsource_transform_disabled.exe $root_flags $fftw_flags
g++ -O2 -std=c++17 -I../levyfit/classes tests/projection_cache.cc ../levyfit/classes/Levy_proj_reader.cc -o /tmp/emissionsource_projection_cache_tests.exe $root_flags -lMinuit2
g++ -O2 -std=c++17 tests/fft_pipeline.cc -o /tmp/emissionsource_fft_tests.exe $root_flags
/tmp/emissionsource_fft_tests.exe "$PWD"
g++ -O2 -std=c++17 tests/fft_systematics.cc -o /tmp/emissionsource_fft_systematics_tests.exe $root_flags
/tmp/emissionsource_fft_systematics_tests.exe "$PWD/calc_and_plot_syserr_FFT.exe"
g++ -O2 -std=c++17 tests/compare_methods.cc -o /tmp/emissionsource_compare_methods_tests.exe $root_flags
/tmp/emissionsource_compare_methods_tests.exe "$PWD/compare_methods.exe"
python3 tests/batch_pipeline.py
