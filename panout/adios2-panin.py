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


def define_variables(io, varinfo_list):
    for vname, shape, start, count, val in varinfo_list:
        io.DefineVariable(vname, val, shape, start, count, ad2.ConstantDims)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", nargs="+", help="infile")
    parser.add_argument("outfile", help="outfile")
    parser.add_argument("--outengine", help="output engine", type=str, default="BPFile")
    parser.add_argument("--var", help="var", nargs="+", type=str)
    parser.add_argument("--decomposition", help="var", nargs="+", type=int)
    parser.add_argument("--append", help="append mode", action="store_true")
    parser.add_argument("--npanout", help="npanout", type=int, default=1)
    parser.add_argument("--start", help="start", type=int, default=0)
    parser.add_argument("-s", "--nstep", help="nstep", type=int, default=1000)
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

    adios = ad2.ADIOS("adios2cfg.xml", ad2.DebugON)
    io = adios.DeclareIO("field3D")
    writer = io.Open(args.outfile, ad2.Mode.Write, comm)

    istep = args.start
    for infile in args.infile:
        with ad2.open(
            infile,
            "r",
            comm=comm,
            config_file="adios2cfg.xml",
            io_in_config_file="field3D",
        ) as fh:
            var_list = list()
            varinfo_list = list()
            if args.var is None:
                var_list.extend(fh.available_variables())
            else:
                var_list.extend(args.var)

            for vname in var_list:
                nstep, nsize = adios2_get_shape(fh, vname)
                ndim = len(nsize)
                start, count = (), ()
                # print (vname, nstep, ndim, nsize)

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

                logging.info((istep, vname, nsize, start, count))
                val = fh.read(vname, start=start, count=count)
                varinfo_list.append((vname, nsize, start, count, val))

            ## Define variables at a first time
            if (istep // args.npanout) == 0:
                define_variables(io, varinfo_list)

            ## Write a step
            writer.BeginStep()
            for vname, shape, start, count, val in varinfo_list:
                var = io.InquireVariable(vname)
                writer.Put(var, val)
            writer.EndStep()
        istep += 1

    ## Output close
    writer.Close()
