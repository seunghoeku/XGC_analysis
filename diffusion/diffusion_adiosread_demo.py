from mpi4py import MPI
import numpy as np
import adios2 as ad2

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


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

    """
    Block-by-block reading
    Assume there is a Adios data written by N writer (N blocks) and we have M readers (M < N)
    Each reader will read about (N/M) number of blocks one by one
    """
    istep = reader.CurrentStep()
    shape_list = adios2_get_block_list(reader, "i_table", istep)
    my_block_list = split_given_size(shape_list, int(np.ceil(len(shape_list) / size)))


    for block in my_block_list[rank]:

        ## Prepare data
        i_table = np.zeros(block["shape"], dtype=np.double)
        var = IO.InquireVariable("i_table")
        var.SetBlockSelection(block["id"])
        reader.Get(var, i_table)
        reader.PerformGets()

        print(
            istep, rank, block["id"], block["shape"], np.min(i_table), np.max(i_table)
        )

    reader.EndStep()

reader.Close()
