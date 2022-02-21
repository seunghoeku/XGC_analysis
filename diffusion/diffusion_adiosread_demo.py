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


def adios2_get_block_list(reader, varname, istep):
    block_list = reader.BlocksInfo(varname, istep)
    shape_list = list()
    for block in block_list:
        blockid = int(block["BlockID"])
        shape = block["Count"]
        lshape = tuple([int(x) for x in shape.strip().split(",")])
        shape_list.append({"id": blockid, "shape": lshape})

    return shape_list


def split_given_size(a, size):
    return np.split(a, np.arange(size, len(a), size))


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
    my_block_list = split_given_size(shape_list, int(np.ceil(len(shape_list) / size)))

    for block in my_block_list[rank]:

        ## Prepare data
        i_ntLV = np.zeros(1, dtype=np.int64)
        i_ntriangles = np.zeros(1, dtype=np.int64)
        i_table = np.zeros(block["shape"], dtype=np.double)

        var_i_ntLV = IO.InquireVariable("i_ntLV")
        var_i_ntriangles = IO.InquireVariable("i_ntriangles")
        var_i_table = IO.InquireVariable("i_table")
        

        var_i_ntLV.SetBlockSelection(block["id"])
        var_i_ntriangles.SetBlockSelection(block["id"])
        var_i_table.SetBlockSelection(block["id"])

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
