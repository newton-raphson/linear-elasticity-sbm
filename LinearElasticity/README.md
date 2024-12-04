
Compilation
============
```bash
mkdir build;
cd build;
```
### For 2D  Computation
```bash
cmake ..  -DENABLE_2D=Yes -DIBM=Yes -DCMAKE_BUILD_TYPE=Release
make -j8 sle-kt
```

### For 3D  Computation
```bash
cmake ..  -DENABLE_3D=Yes  -DIBM=Yes  -DCMAKE_BUILD_TYPE=Release
make -j8 sle-kt
```

HPC
============
### For stampede+frontera (frontera needs to compile in linux, and cannot use unix)
```bash
git submodule update --init --recursive
3D:cmake .. -DCMAKE_CXX_FLAGS="-Wno-unknown-pragmas" -DCMAKE_CXX_COMPILER=mpiicpc -DCMAKE_C_COMPILER=mpiicc -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-xCORE-AVX2 -axCORE-AVX512,MIC-AVX512 -O3 -DNDEBUG" -DENABLE_3D=Yes -DIBM=Yes
2D:cmake .. -DCMAKE_CXX_FLAGS="-Wno-unknown-pragmas" -DCMAKE_CXX_COMPILER=mpiicpc -DCMAKE_C_COMPILER=mpiicc -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-xCORE-AVX2 -axCORE-AVX512,MIC-AVX512 -O3 -DNDEBUG" -DENABLE_2D=Yes -DIBM=Yes
3D (with timer): cmake .. -DCMAKE_CXX_FLAGS="-Wno-unknown-pragmas" -DCMAKE_CXX_COMPILER=mpiicpc -DCMAKE_C_COMPILER=mpiicc -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-xCORE-AVX2 -axCORE-AVX512,MIC-AVX512 -O3 -DNDEBUG" -DENABLE_3D=Yes -DIBM=Yes -DPROFILING=Yes
make -j32 sle-kt
```
    `-xCORE-AVX2 -axCORE-AVX512,MIC-AVX512` are important flags to generate an executable cross compiled for both 
    `skx` and `knl` nodes.


### For nova
```bash
git submodule update --init --recursive
3D: cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_3D=Yes -DCMAKE_C_COMPILEER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DIBM=Yes
2D: cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_2D=Yes -DCMAKE_C_COMPILEER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DIBM=Yes
make -j32 sle-kt
```
