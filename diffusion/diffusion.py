import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib import colors
import adios2
from mpi4py import MPI
from math import sqrt


# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))


class Diffusion:
    def __init__(
        self,
        path="./",
        aggregate=True,
        engine="SST",
        channel_name="diffusion",
        comm=MPI.COMM_WORLD,
    ):
        self.comm = comm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()

        self.aggregate = aggregate

        # read in mesh information from XGC files (written at beginning of XGC sim)
        self.read_mesh(path + "xgc.mesh.bp", path + "xgc.equil.bp")

        # TODO read in timestep
        self.dt = 1e-6

        # create regions to aggregate the tracer data over
        # TODO add in optional region def. to __init__ call, or a config file
        if aggregate:
            self.regions, self.regpsin, self.regtheta = self.create_regions()

        self.setup_adios(engine, channel_name)

    def workflow(self):
        """Continually read in XGC data from adios stream, calc reduced stats, plot"""

        self.dr_stds = []
        self.dr_avgs = []
        self.En_dr_stds = []
        self.En_dr_avgs = []
        self.marker_dens = []

        tindex = 0
        while True:
            stepStatus = self.reader.BeginStep(adios2.StepMode.Read)
            if stepStatus == adios2.StepStatus.OK:
                tindex += 1
                data = self.read_stats()
                # TODO refactor possibly for dictonary
                self.dr_stds.append(data[0])
                self.En_dr_stds.append(data[1])
                self.dr_avgs.append(data[2])
                self.En_dr_avgs.append(data[3])
                self.marker_dens.append(data[4])

                self.reader.EndStep()

                if self.aggregate:
                    sigma_rs, sigma_Ens = self.aggregate_stats()
                    Deff = sigma_rs / (2 * self.dt * tindex)
                    chieff = sigma_Ens / (2 * self.dt * tindex)
                    try:
                        self.plot_regions(Deff, chieff)
                    except Exception as e:
                        print("ERROR: plotting error:", e)
                else:
                    sigma_rs, sigma_Ens = self.calc_sigma(
                        self.dr_stds[tindex],
                        self.En_dr_stds[tindex],
                        self.dr_avgs[tindex],
                        self.En_dr_avgs[tindex],
                        self.marker_dens[tindex],
                    )
                    Deff = sigma_rs / (2 * self.dt * tindex)
                    chieff = sigma_Ens / (2 * self.dt * tindex)
                    self.plot_tri2d(Deff, chieff)
            else:
                if self.rank == 0:
                    print("No more data")
                break

    def setup_adios(self, engine, channel_name):
        """setup adios2 reader"""
        self.adios = adios2.ADIOS("adios2cfg.xml")
        self.IO = self.adios.DeclareIO("tracer_diag")
        self.IO.SetEngine(engine)
        self.reader = self.IO.Open(channel_name, adios2.Mode.Read, self.comm)
        # self.IO.SetParameters()

    ######TODO: Update for ADIOS streaming#############
    def read_stats(self, inds=Ellipsis):
        # these come from XGC gathered locally on the rank, so this is an all-reduce performed by ADIOS2
        # these are not normalized to the 1/N factor, and so need to be done somewhere (N=marker_den)

        ntri = len(self.tri)
        dr_std, En_dr_std, dr_avg, En_dr_avg, marker_den = (
            np.zeros(ntri),
            np.zeros(ntri),
            np.zeros(ntri),
            np.zeros(ntri),
            np.zeros(ntri),
        )

        istep = self.reader.CurrentStep()
        shape_list = adios2_get_block_list(self.reader, "table", istep)
        my_block_list = np.array_split(shape_list, size)

        ## Read block by block
        for block in my_block_list[self.rank]:

            ## Prepare array
            ntLV = np.zeros(1, dtype=np.int64)
            ntriangles = np.zeros(1, dtype=np.int64)
            table = np.zeros(block["shape"], dtype=np.double)

            ## Inquire var info
            var_ntLV = self.IO.InquireVariable("ntLV")
            var_ntriangles = self.IO.InquireVariable("ntriangles")
            var_table = self.IO.InquireVariable("table")

            ## Set block info
            var_ntLV.SetBlockSelection(block["id"])
            var_ntriangles.SetBlockSelection(block["id"])
            var_table.SetBlockSelection(block["id"])

            ## Read
            self.reader.Get(var_ntLV, ntLV)
            self.reader.Get(var_ntriangles, ntriangles)
            self.reader.Get(var_table, table)
            self.reader.PerformGets()

            print(
                istep,
                rank,
                block["id"],
                block["shape"],
                np.min(table),
                np.max(table),
                ntLV.item(),
                ntriangles.item(),
            )

            ## Read each row and update dr/En_dr from table info
            for row in table:
                itri = int(row[0])
                i_dr_average = row[1]
                i_dr_squared_average = row[2]
                i_dE_average = row[3]
                i_dE_squared_average = row[4]
                i_marker_den = row[5]
                e_dr_average = row[6]
                e_dr_squared_average = row[7]
                e_dE_average = row[8]
                e_dE_squared_average = row[9]
                e_marker_den = row[10]

                i_dr_std = sqrt(i_dr_squared_average - i_dr_average**2)
                e_dr_std = sqrt(e_dr_squared_average - e_dr_average**2)
                i_dE_std = sqrt(i_dE_squared_average - i_dE_average**2)
                e_dE_std = sqrt(e_dE_squared_average - e_dE_average**2)

                dr_std[itri] += i_dr_std + e_dr_std
                En_dr_std[itri] += i_dE_std + e_dE_std
                dr_avg[itri] += i_dr_average + e_dr_average
                En_dr_avg[itri] += i_dE_average + e_dE_average
                marker_den[itri] += i_marker_den + e_marker_den

        return (
            dr_std[inds],
            En_dr_std[inds],
            dr_avg[inds],
            En_dr_avg[inds],
            marker_den[inds],
        )

    def read_mesh(self, filename_mesh, filename_eq):
        """Read in mesh info from xgc.mesh.bp"""
        fm = adios2.open(filename_mesh, "r")
        self.RZ = fm.read("/coordinates/values")
        self.tri = fm.read("/cell_set[0]/node_connect_list")
        psi = fm.read("psi")
        fm.close()
        self.triObj = Triangulation(self.RZ[:, 0], self.RZ[:, 1], self.tri)
        feq = adios2.open(filename_eq, "r")
        self.psi_x = feq.read("eq_x_psi")
        feq.close()

        self.psin = psi / self.psi_x
        R0, Z0 = self.RZ[0, :]
        self.theta = np.arctan2(self.RZ[:, 1] - Z0, self.RZ[:, 0] - R0)

        self.RZtri = np.mean(self.RZ[self.tri, :], axis=1)
        self.psintri = np.mean(self.psin[self.tri], axis=-1)
        self.thetatri = np.mean(self.theta[self.tri], axis=-1)

    def create_regions(
        self, theta0=-60, theta1=30, dtheta=10, psin0=1.0, psin1=1.01, dpsin=0.01
    ):
        """Create list of arrays of indices defining regions to aggregate over"""
        Nthetas = int((theta1 - theta0) / dtheta)
        thetas = theta0 + np.arange(Nthetas, dtype=np.double) * dtheta
        thetas *= np.pi / 180.0  # convert deg. to rad.
        Npsins = int((psin1 - psin0) / dpsin)
        psins = psin0 + np.arange(Npsins) * dpsin

        regions = []
        regpsin = []
        regtheta = []
        for t in range(Nthetas - 1):
            for p in range(Npsins - 1):
                inds = np.where(
                    (psins[p] >= self.psin)
                    & (psins[p + 1] < self.psin)(thetas[p] >= self.theta)
                    & (thetas[t + 1] < self.theta)
                )[0]
                regions.append(inds)
                regpsin.append((psins[p + 1] + psins[p]) / 2.0)
                regtheta.append((theta[p + 1] + theta[p]) / 2.0)

        return regions, regpsin, regtheta

    def aggregate_stats(self, tindex=-1):
        """Aggregate the sigmas over regions defined by create_regions"""
        sigma_rs = np.zeros((len(self.regions),))
        sigma_Ens = np.zeros((len(self.regions),))
        for i, region in enumerate(self.regions):
            sigma_rs[i], sigma_Ens[i] = self.calc_sigma(
                self.dr_stds[tindex][region],
                self.En_dr_stds[tindex][region],
                self.dr_avgs[tindex][region],
                self.En_dr_avgs[tindex][region],
                self.marker_dens[tindex][region],
            )
        return sigma_rs, sigma_Ens

    def calc_sigma(dr_std, En_dr_std, dr_avg, En_dr_avg, marker_den):
        """Calculate the true std. dev. of markers"""
        # calculate sigma_r and sigma_En
        sigma_r = dr_std / marker_den - (dr_avg / marker_den) ** 2.0
        sigma_En = En_dr_std / marker_den - (En_dr_avg / marker_den) ** 2.0
        return sigma_r, sigma_En

    def plot_regions(
        self,
        Deff,
        chieff,
        title="log$_{10} \langle (E_n \Delta \sqrt{\psi} - \langle E_n \Delta \sqrt{\psi} \\rangle)^2 \\rangle$",
    ):

        goodinds = Ellipsis  # np.where((En_dr_std>0) & ~np.isinf(En_dr_std) & ~( (self.psintri<1) & (self.RZtri[:,1]<self.eq_x_z) ))[0]

        vmin = np.floor(np.log10(Deff[goodinds]).min())
        vmax = np.ceil(np.log10(Deff[goodinds]).max())
        fig, axs = plt.subplots(2, 1)
        plt.tricontourf(
            self.regpsin[goodinds],
            self.regtheta[goodinds],
            Deff[goodinds],
            100,
            norm=MidpointNormalize(midpoint=0, vmin=vmin, vmax=vmax),
        )
        plt.plot([1.0], [self.thetax], "kx", linewidth=4, markersize=20)
        plt.xlabel("$\psi_N$")
        plt.ylabel("$\\theta$")
        plt.colorbar()
        plt.title(title)
        # En
        plt.tricontourf(
            self.regpsin[goodinds],
            self.regtheta[goodinds],
            chieff[goodinds],
            100,
            norm=MidpointNormalize(midpoint=0, vmin=vmin, vmax=vmax),
        )
        plt.plot([1.0], [self.thetax], "kx", linewidth=4, markersize=20)
        plt.xlabel("$\psi_N$")
        plt.ylabel("$\\theta$")
        plt.colorbar()
        plt.title(title)

    def plot_tri2d(Deff, chieff):
        fig, axs = plt.subplots(2, 1)
        axs[0].tripcolor(self.triObj, Deff)
        axs[1].tripcolor(self.triObj, chieff)


"""
Adios helper functions
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


if __name__ == "__main__":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    diffusion = Diffusion(engine="BP4", channel_name="xgc.tracer_diag.bp")
    diffusion.workflow()
