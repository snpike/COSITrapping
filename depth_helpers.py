#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numba_stats import norm
from scipy import stats
from scipy.integrate import quad
import pandas as pd


global_gamma = 0.50
global_CoverB = 0.13
global_D = 0.028
global_sigma_ratio = 0.85
source_dict = {'Am-241': 59.5409, 'Cs-137': 661.657,'Co-57': 122.06065, 'Ba-133': 356.0129}


# @njit
def threshold(x, x0, sigma, Eth):
    # return (1+stats.norm.cdf(x+Eth, x0, sigma))
    return stats.norm.cdf(x+Eth, x0, sigma)

# @njit
def shelf(x, x0, sigma):
    return (1.-norm.cdf(x, x0, sigma))

def shelf_scipy(x, x0, sigma):
    return (1.-stats.norm.cdf(x, x0, sigma))

# @njit
def exp_tail(x, x0, gamma):
    return (np.exp(gamma*(x-x0)))

# @njit
def linear_tail(x, x0, sigma, m):
    return (1+m*(x-x0))

# @njit
def gaussian(x,x0,sigma):
    return np.exp(-(x-x0)**2/(2*sigma**2))

# @njit
def gauss_plus_tail(x, BoverA, x0, sigma_gauss, gamma, CoverB, D, sigma_ratio):
    return (gaussian(x,x0,sigma_gauss) + BoverA*exp_tail(x, x0, gamma)*shelf(x,x0, sigma_gauss*sigma_ratio) + \
            BoverA*CoverB*linear_tail(x, x0, sigma_gauss*sigma_ratio, D)*shelf(x,x0, sigma_gauss*sigma_ratio))

# # @njit
# def gauss_plus_tail(x, BoverA, x0, sigma_gauss):
#     return gauss_plus_tail_depth(x, BoverA, x0, sigma_gauss, global_gamma, global_CoverB, global_D, global_sigma_ratio)


def gauss_plus_tail_pdf(x, BoverA, x0, sigma_gauss, gamma, CoverB, D, sigma_ratio, Emin = 640., Emax = 672.):
    return gauss_plus_tail(x, BoverA, x0, sigma_gauss, gamma, CoverB, D, sigma_ratio)/\
    quad(gauss_plus_tail, Emin, Emax, args=(BoverA, x0, sigma_gauss, gamma, CoverB, D, sigma_ratio))[0]

def make_df_from_dat(file, e_min = 640., e_max = 672.):
    
    # SH det  side  strip_num (starting at 1)  did strip trigger (0 or 1)  
           #timing raw AD counter units  corrected AD counter units  ADC  energy  ??  
    # HT x y z energy [cm, cm, cm, keV]   

    rows = []
    col_names = ["ID","det","strip_p","energy_p","strip_n","energy_n","x","y","z"]

    with open(file, "r") as f:
        
        ev_block = []
        
        for line in f:
        #for each line, start a block of lines corresponding to an event
            
            if line.startswith('SE'):
            # If the accumulated block has 6 lines then it's a single-pixel event
                
                if len(ev_block) == 6:
                    if np.product(np.array([("BD" not in l) for l in ev_block])):
                        
                        ID = int(ev_block[0].split(" ")[-1].strip("\n"))

                        if len(ev_block[3].split(" ")) < 9:
                               print(ev_block)
                        if len(ev_block[4].split(" ")) < 9:
                               print(ev_block)

                        energy_p = float(ev_block[3].split(" ")[8])
                        energy_n = float(ev_block[4].split(" ")[8])

                        # select photopeak events
                        if (energy_p < e_max and energy_p > e_min) and (energy_n < e_max and energy_n > e_min):

                            # save info from SH p line
                            det = int(ev_block[3].split(" ")[1])
                            strip_p = int(ev_block[3].split(" ")[3])

                            # save info from SH n line
                            strip_n = int(ev_block[4].split(" ")[3])

                            # save position [cm] info from HT line
                            x = float(ev_block[5].split(" ")[1])
                            y = float(ev_block[5].split(" ")[2])
                            z = float(ev_block[5].split(" ")[3])

                            # save info to df
                            columns = [ID,det,strip_p,energy_p,strip_n,energy_n,x,y,z]
                            rows.append(columns)
                ev_block = []
                
            else:
                ev_block.append(line)
    
    df = pd.DataFrame(rows,columns=col_names)
    return df

def make_df_from_dat_CC(file, e_min = 640., e_max = 672.):
    
    # SH det  side  strip_num (starting at 1)  did strip trigger (0 or 1)  
           #timing raw AD counter units  corrected AD counter units  ADC  energy  ??  
    # HT x y z energy [cm, cm, cm, keV]   

    rows = []
    col_names = ["ID","det","strip_p","energy_p","strip_n","energy_n"]

    with open(file, "r") as f:
        
        ev_block = []
        
        for line in f:
        #for each line, start a block of lines corresponding to an event
            
            if line.startswith('SE'):
            # If the accumulated block has 6 lines then it's a single-pixel event
                
                if len(ev_block) == 8:
                    if np.product(np.array([("StripPairingIncomplete" not in l) for l in ev_block])):
                        
                        ID = int(ev_block[0].split(" ")[-1].strip("\n"))

                        if len(ev_block[3].split(" ")) < 9:
                               print(ev_block)
                        if len(ev_block[4].split(" ")) < 9:
                               print(ev_block)

                        energy_p = float(ev_block[3].split(" ")[8])
                        energy_n = float(ev_block[4].split(" ")[8])

                        # select photopeak events
                        if (energy_p < e_max and energy_p > e_min) and (energy_n < e_max and energy_n > e_min):

                            # save info from SH p line
                            det = int(ev_block[3].split(" ")[1])
                            strip_p = int(ev_block[3].split(" ")[3])

                            # save info from SH n line
                            strip_n = int(ev_block[4].split(" ")[3])

                            # save position [cm] info from HT line
                            # x = float(ev_block[5].split(" ")[1])
                            # y = float(ev_block[5].split(" ")[2])
                            # z = float(ev_block[5].split(" ")[3])

                            # save info to df
                            columns = [ID,det,strip_p,energy_p,strip_n,energy_n]
                            rows.append(columns)
                ev_block = []
                
            else:
                ev_block.append(line)
    
    df = pd.DataFrame(rows,columns=col_names)
    return df