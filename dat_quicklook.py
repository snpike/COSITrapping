#!/usr/bin/env python

### This script takes in a calibrated .dat file and produces plots of the spectra.
### Written by Sean Pike in Augus 2024

import numpy as np
import matplotlib.pyplot as plt
import sys, getopt
import os
import argparse
from depth_helpers import *

sns.set_context('talk', font_scale=1.0)
sns.set_palette('colorblind')


def main(argv):
    infile=''
    emin = 0.0
    emax = 2000.
    opts, args = getopt.getopt(argv, "hi:f", ['infile=','emin=','emax='])
    filelist = []
    file_prefix = ''
    for opt, arg in opts:
        if opt == '-h':
            print('python dat_quicklook.py -i <inputfile> --emin <emin> --emax <emax>')
        elif opt in ('-i', '--infile'):
            filelist.append(arg)
        elif opt=='--emin':
            emin=float(arg)
        elif opt=='--emax':
            emax=float(arg)
        else:
            assert False, "incorrect option"
    if emax < emin:
        assert False, "The maximum energy is less than the minimum"
    data_df = make_df_from_dat(filelist, e_min = 0.0, e_max = 2000.)
    fig_p = plt.figure()
    hist_p, bin_edges_p = plt.hist(data_df['energy_p'], bins = 6000, range=(0.,2000.))
    plt.xrange((emin, emax))
    ymax = np.max(hist_p[argmin(np.abs(bin_edges_p-emin)):argmin(np.abs(bin_edges_p-emax))])
    plt.yrange((0.0, ymax*1.2))
    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts')
    plt.tight_layout()
    plt.show()

    fig_n = plt.figure()
    hist_n, bin_edges_n = plt.hist(data_df['energy_n'], bins = 6000, range=(0.,2000.))
    plt.xrange((emin, emax))
    ymax = np.max(hist_n[argmin(np.abs(bin_edges_n-emin)):argmin(np.abs(bin_edges_n-emax))])
    plt.yrange((0.0, ymax*1.2))
    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts')
    plt.tight_layout()
    plt.show()


if __name__=="__main__":
    main(sys.argv[1:])