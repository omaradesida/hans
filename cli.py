
import argparse

def parse_args():

    """Parse command line arguments. Returns an empty Namespace containing the attributes."""

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--iterations", default = 2e5, type=float, help = "Number of iterations to run \n")
    parser.add_argument("-t", "--time", type=float, help = "How much time has been allocated for the program to run")
    parser.add_argument("-w", "--nwalkers", default = 10, type=int, help = "Number of walkers to use \n")
    parser.add_argument("-c", "--nchains", default = 32, type=int, help = "Number of chains in the simulation box \n")
    parser.add_argument("-b", "--nbeads", default = 1, type=int, help = "Number of beads in each chain \n")
    parser.add_argument("-a","--bondangle",default = 109.47,type=float, help ="Angle between consecutive beads on a chain.\n")
    parser.add_argument("-l", "--walklength", default = 80, type=int, help = "Number of sweep for each random walk step \n")
    parser.add_argument("-p","--processes",default = 1,type=int, help = "Number of processes to use in parallel when performing random walks.\n")
    parser.add_argument("-d", "--bondlength", default = 0.4, type=float, help="Distance between hard spheres in a chain")
    parser.add_argument("-R", "--restart", action = "store_true", help = "Whether or not to restart from a previous attempt\n")
    parser.add_argument("-f","--restart_folder",type=str, help =
                        "Folder where the restart file is located. In the folder, a restart.hdf5 file must be present in order to resume the simulation.\n")


    return parser.parse_args()
