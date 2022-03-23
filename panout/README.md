# Panout
Read xgc.3d.bp and split the stream (or file) in a round-robin fashion.

# Run on Summit
Below is an example to run on Summit.

```
module load gcc/9.1.0

module use -a /gpfs/alpine/phy122/world-shared/lib/install/summit/modulefiles
module load adios2/devel
module load gptl4py

SRC=/gpfs/alpine/proj-shared/csc143/jyc/summit/test_GB_small_su455/xgc.3d.bp

LD_PRELOAD=/sw/summit/gcc/9.1.0-alpha+20190716/lib64/libstdc++.so:/sw/summit/gcc/9.1.0-alpha+20190716/lib64/libgomp.so \
jsrun -n 8 python adios2-panout.py --npanout=6 $SRC xgc.3d.panout
```
