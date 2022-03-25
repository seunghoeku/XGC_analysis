# XGC_analysis
Code for analysis node from XGC particle data

# Build (Summit)
```
module purge
ml DefApps
ml gcc

module load cmake/3.20.2
module load boost
module load python/3.8-anaconda3
module load libfabric/1.12.1-sysrdma

ADIOS_DIR=/gpfs/alpine/world-shared/phy122/lib/install/summit/adios2/devel/gcc9.1.0
CAMTIMER_DIR=/gpfs/alpine/world-shared/phy122/lib/install/summit/camtimers/gcc9.1.0

CC=gcc CXX=g++ cmake -DCMAKE_PREFIX_PATH="$ADIOS_DIR;$CAMTIMER_DIR" ..

```

# Command line options
```
$ diffusion/diffusion -h
Allowed options:
  -h [ --help ]             produce help message
  -w [ --xgcdir ] arg       XGC directory
  -s [ --maxstep ] arg (=0) max steps

$ ./heatload/heatload -h
Allowed options:
  -h [ --help ]             produce help message
  -w [ --xgcdir ] arg       XGC directory
  -s [ --maxstep ] arg (=0) max steps
  -f [ --freshstart ]       fresh start (no restart)
  -i [ --ion_only ]         Ion only

$ python ../panout/adios2-panout.py -h
usage: adios2-panout.py [-h] [--outengine OUTENGINE] [--var VAR [VAR ...]] [--decomposition DECOMPOSITION [DECOMPOSITION ...]] \
                             [--append] [--npanout NPANOUT] [--start START] [-s NSTEP] infile outfile

positional arguments:
  infile                infile
  outfile               outfile

optional arguments:
  -h, --help            show this help message and exit
  --outengine OUTENGINE
                        output engine
  --var VAR [VAR ...]   var
  --decomposition DECOMPOSITION [DECOMPOSITION ...]
                        var
  --append              append mode
  --npanout NPANOUT     npanout
  --start START         start
  -s NSTEP, --nstep NSTEP
                        nstep

```

# Run (Summit)
```
#!/bin/bash
#BSUB -env "all"
#BSUB -P PHY122
#BSUB -W 0:30
#BSUB -nnodes 1
#BSUB -q debug
#BSUB -o run-%J.log
#BSUB -e run-%J.log

[ -z $JOBID ] && JOBID=$LSB_JOBID
[ -z $JOBSIZE ] && JOBSIZE=$(((LSB_DJOB_NUMPROC-1)/42))

export OMPI_MCA_coll_ibm_collselect_mode_barrier=failsafe

module purge
ml DefApps
ml gcc

module load cmake/3.20.2
module load forge
module load boost
module load python/3.8-anaconda3
module load libfabric/1.12.1-sysrdma

module use -a /gpfs/alpine/world-shared/phy122/lib/install/summit/modulefiles
module load adios2/devel
module load gptl4py/devel

NR=6
XGCDIR=/gpfs/alpine/world-shared/phy122/sku/GB_ITER_1024

echo "Run diffusion"
export OMP_NUM_THREADS=1
jsrun -n $((JOBSIZE*NR)) -a1 -c1 -g0 -r$NR -brs /usr/bin/stdbuf -oL -eL ./diffusion/diffusion -w $XGCDIR 2>&1 | tee run-diffusion-$JOBID.log
echo "Done."
sleep 3

echo "Run heatload"
export OMP_NUM_THREADS=14
jsrun -n $((JOBSIZE*NR)) -a1 -c7 -g0 -r$NR -brs /usr/bin/stdbuf -oL -eL ./heatload/heatload -w $XGCDIR -f 2>&1 | tee run-heatload-$JOBID.log
echo "Done."
sleep 3

echo "Run pan-out"
LD_PRELOAD=/sw/summit/gcc/9.1.0-alpha+20190716/lib64/libstdc++.so:/sw/summit/gcc/9.1.0-alpha+20190716/lib64/libgomp.so \
jsrun -n $((JOBSIZE*NR)) python adios2-panout.py --npanout=6 $XGCDIR/xgc.3d.bp xgc.3d.panout
echo "Done."

```

# Additional Input files for Gordon Bell
located at
```
/gpfs/alpine/world-shared/phy122/sku/GB_XGC_input
```
sml_input_dir in 'input' should point this directory 
