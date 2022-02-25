# poincare
Code for poincare surface creation using VTK-m

# Summit Build
```
module load git-lfs
module load gcc/9.1.0
module load cuda/11.0.3
module load cmake/3.21.3
module load adios2

mkdir build
cd build

cmake \
 -DCMAKE_CXX_COMPILER=/sw/summit/gcc/9.1.0-alpha+20190716/bin/g++ \
 -DCMAKE_C_COMPILER=/sw/summit/gcc/9.1.0-alpha+20190716/bin/gcc \
 -DCMAKE_BUILD_TYPE=Release \
 -DVTKm_NO_DEPRECATED_VIRTUAL=ON \
 -DVTKm_ENABLE_TESTING=OFF \
 -DBUILD_TESTING=OFF \
 -DVTKm_ENABLE_EXAMPLES=ON\
 -DVTKm_ENABLE_CUDA=ON \
 -DVTKm_ENABLE_OPENMP=ON \
 -DVTKm_ENABLE_MPI=OFF \
 -DADIOS2_DIR=/sw/summit/spack-envs/base/opt/linux-rhel8-ppc64le/gcc-7.5.0/adios2-2.7.1-qtvhwjxzdlbro3zv7zytmxm72fc3a5d5/lib64/cmake/adios2 \
 -DVTKm_USE_DOUBLE_PRECISION=ON \
../vtk-m

make

```

# Run on Summit
Below is an example to run on a cyclone dataset on Summit.

```
jsrun -n1 -a1 -c20 -g1 ./build/examples/poincare/Simple2.3  --vField B --dir /gpfs/alpine/proj-shared/csc143/pugmire/summit/vtkm/data/POINC --traces 0 --useHighOrder --turbulence 1 --range 2.5 3.5 --numSeeds 10000 --gpu  --worklet 2 --output OUTPUT --numPunc 200 --gpuParams 256 128 --stepSize 0.1

```