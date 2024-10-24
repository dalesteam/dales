"""Converts prof.inp and lscale.inp input files to new NetCDF-based input."""
import argparse

import numpy as np
import netCDF4 as nc

PROF_NAMES = ["zh", "thetal", "qt", "ua", "va", "tke"]
LSCALE_NAMES = ["zh", "ug", "vg", "wa", "dqtdxls", "dqtdyls", "tnqt_adv", "tnthetal_rad"]


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--prof", type=str, default="prof.inp.001")
    parser.add_argument("-l", "--lscale", type=str, default="lscale.inp.001")
    parser.add_argument("-o", "--output", type=str, default="init.001.nc")

    return parser.parse_args()


def main():
    args = parse_args()

    prof = np.genfromtxt(args.prof, skip_header=2)
    lscale = np.genfromtxt(args.lscale, skip_header=2)

    with nc.Dataset(args.output, "w") as ds:
        height = ds.createDimension("zh", len(prof))
        heights = ds.createVariable("zh", "f", ("zh",))
        heights[:] = prof[:,0]

        for ivar in range(1, len(prof[0])):
            nc_var = ds.createVariable(PROF_NAMES[ivar], "f", ("zh",))
            nc_var[:] = prof[:,ivar]

        for ivar in range(1, len(lscale[0])):
            nc_var = ds.createVariable(LSCALE_NAMES[ivar], "f", ("zh",))
            nc_var[:] = lscale[:,ivar]


if __name__ == "__main__":
    main()
