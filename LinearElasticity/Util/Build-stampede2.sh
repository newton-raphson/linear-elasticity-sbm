#module load intel/17.0.4
#module load petsc/3.7-cxx
#module load cmake/3.10.2

#cmake_param="-DCMAKE_CXX_FLAGS=\"-Wno-unknown-pragmas\" -DCMAKE_CXX_COMPILER=mpiicpc -DCMAKE_C_COMPILER=mpiicc -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=\"-xCORE-AVX2 -axCORE-AVX512,MIC-AVX512 -O3 -DNDEBUG\""
cd ../cmake-build-release
cmake .. -DCMAKE_CXX_FLAGS="-Wno-unknown-pragmas" -DCMAKE_CXX_COMPILER=mpiicpc -DCMAKE_C_COMPILER=mpiicc -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-xCORE-AVX2 -axCORE-AVX512,MIC-AVX512 -O3 -DNDEBUG"
make -j 16 nshtibm-kt
cd ../cmake-build-release-2d
cmake .. -DCMAKE_CXX_FLAGS="-Wno-unknown-pragmas" -DCMAKE_CXX_COMPILER=mpiicpc -DCMAKE_C_COMPILER=mpiicc -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-xCORE-AVX2 -axCORE-AVX512,MIC-AVX512 -O3 -DNDEBUG" -DENABLE_2D=Yes
make -j 16 nshtibm-kt
