import matplotlib.pyplot as plt
import argparse
import scipy as sp
import numpy as np
import os

def calc_pgs(direc, ofile, title):
    np_files = []
    for x in os.listdir(direc):
        if x.endswith(".npy"):
            np_files.append(x)

    sfs_files = [f for _,f in enumerate(np_files) if f.startswith("sfs")]

    sfs_abs = np.load(sfs_files[0])
    sfs_cst = np.load(sfs_files[1])
    sfs_onf = np.load(sfs_files[2])


    #for i in range(4):
     #   plt.hist(sfs_onf[i]/np.sum(sfs_onf[i]), rwidth=0.18)
     #   plt.legend()
      #  plt.savefig(str(i)+ofile)
    onf_nmg = sfs_onf[0][1:100]
    onf_ud = sfs_onf[1][1:100]
    onf_rc = sfs_onf[2][1:100]
    onf_dm = sfs_onf[3][1:100]

    rw = 0.1

    plt.hist(onf_nmg, label='No Migration', rwidth=rw)
    plt.hist(onf_ud, label='Unidirectional Migration', rwidth=rw)
    plt.hist(onf_rc, label='Reciprocal Migration', rwidth=rw)
    plt.hist(onf_dm, label='Directional Migration', rwidth=rw)
    plt.title(title)
    plt.xlabel("Derived Allele Frequency")
    plt.ylabel("Counts")
    plt.legend(loc='upper right')
    plt.savefig(ofile)



def main():
    parser = argparse.ArgumentParser(
        description="Compute number of Segregating sites and SFS plots in Simulations.")

    parser.add_argument("-dir", action="store", dest="dir", type=str,
                        default="/Users/mohamedabdelrahman/Desktop")
    parser.add_argument("-o", action="store", dest="out_file", type=str,
                        default="sfsresults.png")
    parser.add_argument("-t", action="store", dest="title", type=str,
                        default="ONF SFS")
    args = parser.parse_args()

    direc = args.dir
    ofile = args.out_file
    title = args.title
    calc_pgs(direc, ofile, title)

if __name__ == '__main__':
    main()