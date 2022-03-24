# XGC_analysis
Code for analysis node from XGC particle data

# Build (Summit)
```
module load boost

ADIOS_DIR=/gpfs/alpine/world-shared/phy122/lib/install/summit/adios2/devel/gcc9.1.0
CAMTIMER_DIR=/gpfs/alpine/world-shared/phy122/lib/install/summit/camtimers/gcc9.1.0

CC=gcc CXX=g++ cmake -DCMAKE_PREFIX_PATH="$ADIOS_DIR;$CAMTIMER_DIR" ..

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

date

NR=6
XGCDIR=/gpfs/alpine/world-shared/phy122/sku/GB_ITER_1024

echo "Run diffusion"
export OMP_NUM_THREADS=1
jsrun -n $((JOBSIZE*NR)) -a1 -c1 -g0 -r$NR -brs /usr/bin/stdbuf -oL -eL ./diffusion/diffusion -w $XGCDIR -s 5 2>&1 | tee run-diffusion-$JOBID.log
echo "Done."
sleep 3

echo "Run heatload"
export OMP_NUM_THREADS=14
jsrun -n $((JOBSIZE*NR)) -a1 -c7 -g0 -r$NR -brs /usr/bin/stdbuf -oL -eL ./heatload/heatload -w $XGCDIR -s 5 -f 2>&1 | tee run-heatload-$JOBID.log
echo "Done."
sleep 3

echo "Run pan-out"
LD_PRELOAD=/sw/summit/gcc/9.1.0-alpha+20190716/lib64/libstdc++.so:/sw/summit/gcc/9.1.0-alpha+20190716/lib64/libgomp.so \
jsrun -n $((JOBSIZE*NR)) python adios2-panout.py --npanout=6 -s 5 xgc.3d.bp xgc.3d.panout

```

# Additional Input files for Gordon Bell
located at
```
/gpfs/alpine/world-shared/phy122/sku/GB_XGC_input
```
sml_input_dir in 'input' should point this directory 
