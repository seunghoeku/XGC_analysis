# XGC_analysis
Code for analysis node from XGC particle data

# Build (Summit)
```
module load boost

ADIOS_DIR=/gpfs/alpine/world-shared/phy122/lib/install/summit/adios2/devel/gcc9.1.0
CAMTIMER_DIR=/gpfs/alpine/world-shared/phy122/lib/install/summit/camtimers/gcc9.1.0

CC=gcc CXX=g++ cmake -DCMAKE_PREFIX_PATH="$ADIOS_DIR;$CAMTIMER_DIR" ..

```

# Additional Input files for Gordon Bell
located at
```
/gpfs/alpine/world-shared/phy122/sku/GB_XGC_input
```
sml_input_dir in 'input' should point this directory 
