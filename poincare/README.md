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
Running cyclone case:

jsrun -n1 -a1 -c20 -g1 ./examples/poincare/Poincare  --vField B --dir ../data/run_1920 --traces 0 --useHighOrder --turbulence 1 --psiRange .05 .95 10 8 --gpu  --output OUTPUT --numPunc 200 --gpuParams 256 128 --stepSize 0.1 --test

Running Streaming case:
Reading panout.1

jsrun -n1 -a1 -c20 -g1 ./examples/poincare/Poincare  --dir <data directory> --gpu  --output OUT.1 --numPunc 100 --gpuParams 256 128 --stepSize 0.05 --streaming panout.1.bp  --psiRange 0.95 1.05 100 12 --xml ../poincare.xml


Note: this uses --useLinearB while I fix the issue with interpolation.

Running on whoopingcough:

./examples/poincare/Poincare --dir /media/dpn/disk2TB/proj/vtkm/XGCLocator/data/XGC_GB/su458_ITER_data --useHighOrder --turbulence 1 --openmp  --output OUT --numPunc 100 --gpuParams 256 128 --stepSize 0.01 --dumpSeeds --parseAdios xgc.particle.0000000.init.bp 1 |tee out



```