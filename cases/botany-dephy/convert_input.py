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

        # Nudging
        dim = ds.createDimension("time")
        var = ds.createVariable("time", "f", ("time",)) 
        var[0] = 0
        var[1] = 1.00000000E+07

        a = 2
        b = 3
        c = 7.4
        nudging_constant = 6 * 3600 + (b * (0.5 * np.pi + np.arctan(a * 0.5 * np.pi * (1 - prof[:,0] / 3000))))**c

        var = ds.createVariable("ua_nud", "f", ("time", "zh"))
        var[0,:] = prof[:,3]
        var[1,:] = prof[:,3]

        var = ds.createVariable("nudging_constant_ua", "f", ("time", "zh"))
        var[0,:] = nudging_constant
        var[1,:] = nudging_constant

        var = ds.createVariable("thetal_nud", "f", ("time", "zh"))
        var[0,:] = prof[:,1]
        var[1,:] = prof[:,1]

        var = ds.createVariable("nudging_constant_thetal", "f", ("time", "zh"))
        var[0,:] = nudging_constant
        var[1,:] = nudging_constant

        var = ds.createVariable("qt_nud", "f", ("time", "zh"))
        var[0,:] = prof[:,2]
        var[1,:] = prof[:,2]

        var = ds.createVariable("nudging_constant_qt", "f", ("time", "zh"))
        var[0,:] = nudging_constant
        var[1,:] = nudging_constant


if __name__ == "__main__":
    main()
