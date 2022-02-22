from mpi4py import MPI
import numpy as np
import adios2 as ad2

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

"""
Block-by-block read demonstration:
Let assume there is a Adios file written by N writers (N blocks) and we have M readers (M < N).
Each reader will read about (N/M) number of blocks one by one.
"""


def range_split(n, size, rank):
    return [slice(x[0], x[-1] + 1) for x in np.array_split(range(n), size)][rank]


def adios2_block_read(IO, reader, varname, slice=None, dtype=np.double):
    if slice is None:
        slice = slice(None, None, None)
    istep = reader.CurrentStep()
    block_list = reader.BlocksInfo(varname, istep)
    arr_list = list()
    for block in block_list[slice]:
        shape = tuple([int(x) for x in block["Count"].strip().split(",")])
        arr = np.zeros(shape, dtype=dtype)
        arr_list.append(arr)

        var = IO.InquireVariable(varname)
        var.SetBlockSelection(int(block["BlockID"]))
        reader.Get(var, arr)

    return arr_list


def adios2_get_block_list(reader, varname, istep):
    block_list = reader.BlocksInfo(varname, istep)
    shape_list = list()
    for block in block_list:
        blockid = int(block["BlockID"])
        shape = block["Count"]
        lshape = tuple([int(x) for x in shape.strip().split(",")])
        shape_list.append({"id": blockid, "shape": lshape})

    return shape_list


ah = ad2.ADIOS(MPI.COMM_WORLD)
IO = ah.DeclareIO("diffusion_read")
IO.SetEngine("BP4")
reader = IO.Open("xgc.tracer_diag.su455_table.bp", ad2.Mode.Read)

while True:
    status = reader.BeginStep(ad2.StepMode.Read, timeoutSeconds=60)
    if status != ad2.StepStatus.OK:
        if rank == 0:
            print("No more data")
        break

    istep = reader.CurrentStep()
    shape_list = adios2_get_block_list(reader, "i_table", istep)
    my_block_list = np.array_split(shape_list, size)

    ## Read block by block
    for block in my_block_list[rank]:

        ## Prepare array
        i_ntLV = np.zeros(1, dtype=np.int64)
        i_ntriangles = np.zeros(1, dtype=np.int64)
        i_table = np.zeros(block["shape"], dtype=np.double)

        ## Inquire var info
        var_i_ntLV = IO.InquireVariable("i_ntLV")
        var_i_ntriangles = IO.InquireVariable("i_ntriangles")
        var_i_table = IO.InquireVariable("i_table")

        ## Set block info
        var_i_ntLV.SetBlockSelection(block["id"])
        var_i_ntriangles.SetBlockSelection(block["id"])
        var_i_table.SetBlockSelection(block["id"])

        ## Read
        reader.Get(var_i_ntLV, i_ntLV)
        reader.Get(var_i_ntriangles, i_ntriangles)
        reader.Get(var_i_table, i_table)
        reader.PerformGets()

        print(
            istep,
            rank,
            block["id"],
            block["shape"],
            np.min(i_table),
            np.max(i_table),
            i_ntLV.item(),
            i_ntriangles.item(),
        )

    reader.EndStep()

reader.Close()
