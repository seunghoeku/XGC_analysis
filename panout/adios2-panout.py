#!/usr/bin/env python
from mpi4py import MPI
import numpy as np
import adios2 as ad2
import os
import argparse
import logging
import itertools

try:
    import gptl4py as gp
except ImportError:

    class gp:
        def initialize():
            pass

        def finalize():
            pass

        def start(name):
            pass

        def stop(name):
            pass

        def pr_file(filename):
            pass

        def pr_summary_file(filename):
            pass


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
    parser.add_argument("infile", help="infile")
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

    fname_list = list()
    writer_list = list()
    for i in range(args.npanout):
        io = adios.DeclareIO("field3D.%d" % i)
        fname = "%s.%d.bp" % (args.outfile, i)
        fname_list.append(fname)
        writer = io.Open(fname, ad2.Mode.Write, comm)
        writer_list.append((io, writer))

    gp.initialize()
    gp.start("TOTAL")
    with ad2.open(
        args.infile,
        "r",
        comm=comm,
        config_file="adios2cfg.xml",
        io_in_config_file="field3D",
    ) as fh:
        for fstep in fh:
            istep = fstep.current_step()
            if (istep < args.start) or (istep >= args.start + args.nstep):
                continue

            gp.start("STEP")
            gp.start("ADIOS_STEP")
            var_list = list()
            varinfo_list = list()
            if args.var is None:
                var_list.extend(fstep.available_variables())
            else:
                var_list.extend(args.var)

            for vname in var_list:
                nstep, nsize = adios2_get_shape(fstep, vname)
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
                if np.prod(count) > 0:
                    gp.start("ADIOS_PERFORM_GETS")
                    val = fstep.read(vname, start=start, count=count)
                    gp.stop("ADIOS_PERFORM_GETS")
                    varinfo_list.append((vname, nsize, start, count, val))

            gp.stop("ADIOS_STEP")

            ## Output
            logging.info(f"Writing: %s" % (fname_list[istep % args.npanout]))
            io, writer = writer_list[istep % args.npanout]

            ## Define variables at a first time
            if (istep // args.npanout) == 0:
                define_variables(io, varinfo_list)

            ## Write a step
            gp.start("ADIOS_WRITE")
            writer.BeginStep()
            for vname, shape, start, count, val in varinfo_list:
                var = io.InquireVariable(vname)
                writer.Put(var, val)
            writer.EndStep()
            gp.stop("ADIOS_WRITE")
            gp.stop("STEP")

    ## Output close
    for io, writer in writer_list:
        writer.Close()

    gp.stop("TOTAL")
    gp.pr_file("adios2-panout-timing.%d" % rank)
    gp.pr_summary_file("adios2-panout-timing.summary")
    gp.finalize()
