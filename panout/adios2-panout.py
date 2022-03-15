#!/usr/bin/env python
from mpi4py import MPI
import numpy as np
import adios2 as ad2
import os
import argparse
import logging
import itertools


def split(nsize, nchunk, i):
    assert i < nchunk
    lsize = nsize // nchunk if nchunk > (i + 1) else nsize - (nsize // nchunk) * i
    offset = (nsize // nchunk) * i
    return (offset, lsize)


def adios2_get_shape(f, varname):
    nstep = int(f.available_variables()[varname]["AvailableStepsCount"])
    shape = f.available_variables()[varname]["Shape"]
    lshape = None
    if shape == "":
        ## Accessing Adios1 file
        ## Read data and figure out
        v = f.read(varname)
        lshape = v.shape
    else:
        lshape = tuple([int(x.strip(",")) for x in shape.strip().split()])
    return (nstep, lshape)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="infile")
    parser.add_argument("outfile", help="outfile")
    parser.add_argument("--outengine", help="output engine", type=str, default="BPFile")
    parser.add_argument("--var", help="var", nargs="+", type=str)
    parser.add_argument("--decomposition", help="var", nargs="+", type=int)
    parser.add_argument("--append", help="append mode", action="store_true")
    parser.add_argument("--npanout", help="npanout", type=int, default=1)
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    fmt = "[%d:%%(levelname)s] %%(message)s" % (rank)
    logging.basicConfig(level=logging.DEBUG, format=fmt)

    logging.info(f"rank,size: {rank} {size}")
    if rank == 0:
        logging.info(f"infile: {args.infile}")
        logging.info(f"outfile: {args.outfile}")
        logging.info(f"var: {args.var}")
        logging.info(f"decomposition: {args.decomposition}")

    decomposition = [
        1,
    ] * 10
    if args.decomposition is not None:
        for i in range(len(args.decomposition)):
            decomposition[i] = args.decomposition[i]
    else:
        decomposition[0] = size

    assert size == np.prod(decomposition)

    with ad2.open(
        args.infile,
        "r",
        comm=comm,
        config_file="adios2cfg.xml",
        io_in_config_file="field3D",
    ) as fh:
        for fstep in fh:
            istep = fstep.current_step()
            var_list = list()
            vars = list()
            if args.var is None:
                var_list.extend(fstep.available_variables())
            else:
                var_list.extend(args.var)

            for var in var_list:
                nstep, nsize = adios2_get_shape(fstep, var)
                ndim = len(nsize)
                start, count = (), ()
                # print (var, nstep, ndim, nsize)

                if ndim > 0:
                    x = list()
                    for i in range(ndim):
                        y = list()
                        for j in range(decomposition[i]):
                            s = split(nsize[i], decomposition[i], j)
                            y.append(s)
                        x.append(y)
                    # print(x)

                    z = list(itertools.product(*x))[rank]
                    start, count = list(zip(*z))

                logging.info((istep, var, nsize, start, count))
                val = fstep.read(var, start=start, count=count)
                vars.append((var, nsize, start, count, val))

            fname = "%s.%d.bp" % (args.outfile, istep % args.npanout)
            logging.info("Writing: %s" % fname)
            mode = "w" if istep < args.npanout else "a"
            with ad2.open(
                fname, mode=mode, comm=comm, engine_type=args.outengine
            ) as fw:
                for var, shape, start, count, val in vars:
                    fw.write(var, val, shape, start, count)
