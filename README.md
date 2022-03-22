# XGC_analysis
Code for analysis node from XGC particle data

# Build (Summit)
```
module load boost

ADIOS2_DIR=/gpfs/alpine/world-shared/phy122/lib/install/summit/adios2-nompi/devel/nvhpc \
CAMTIMER_DIR=/gpfs/alpine/world-shared/phy122/lib/install/summit/camtimers/gcc9.1.0 \
CC=gcc CXX=g++ cmake -DCMAKE_PREFIX_PATH=$CAMTIMER_DIR ..
```

# Additional Input files for Gordon Bell
located at
```
/gpfs/alpine/world-shared/phy122/sku/GB_XGC_input
```
sml_input_dir in 'input' should point this directory 
