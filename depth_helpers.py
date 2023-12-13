#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numba_stats import norm
from scipy import stats
from scipy.integrate import quad
import pandas as pd
from scipy.interpolate import interp1d, CubicSpline


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
def linear_tail(x, x0, m):
    return (1+m*(x-x0))

# @njit
def gaussian(x,x0,sigma):
    return np.exp(-(x-x0)**2/(2*sigma**2))

# @njit
def gauss_plus_tail(x, BoverA, x0, sigma_gauss, gamma, CoverB, D, sigma_ratio):
    return (gaussian(x,x0,sigma_gauss) + BoverA*exp_tail(x, x0, gamma)*shelf(x,x0, sigma_gauss*sigma_ratio) + \
            BoverA*CoverB*linear_tail(x, x0, D)*shelf(x,x0, sigma_gauss*sigma_ratio))

# # @njit
# def gauss_plus_tail(x, BoverA, x0, sigma_gauss):
#     return gauss_plus_tail_depth(x, BoverA, x0, sigma_gauss, global_gamma, global_CoverB, global_D, global_sigma_ratio)


def gauss_plus_tail_pdf(x, BoverA, x0, sigma_gauss, gamma, CoverB, D, sigma_ratio, Emin = 640., Emax = 672.):
    return gauss_plus_tail(x, BoverA, x0, sigma_gauss, gamma, CoverB, D, sigma_ratio)/\
    quad(gauss_plus_tail, Emin, Emax, args=(BoverA, x0, sigma_gauss, gamma, CoverB, D, sigma_ratio))[0]

class DepthCalibrator_Am241:
    ### TODO: unbinned analysis with noise.
    def __init__(self, AC_param_file, DC_param_file, sim_file, DC_sim_mean = 186.1, AC_sim_mean = -182.6, savefile=None):
        ### take in the Gaussian means from AC- and DC- side illumination with Am241 and compare to simulations to return a CTD->depth mapping function.
        self.AC_params = np.array([[0.0 for p in range(37)] for n in range(37)])
        self.DC_params = np.array([[0.0 for p in range(37)] for n in range(37)])
        self.DC_sim_mean=DC_sim_mean
        self.AC_sim_mean=AC_sim_mean

        ### Read in the Gaussian centroids determined for AC side illumination using dat2CTD.py
        with open(AC_param_file) as file:
            for line in file:
                if '#' not in line:
                    splitline = line.split(', ')
                    p = int(splitline[0])
                    n = int(splitline[1])
                    self.AC_params[p][n] = float(splitline[2])

        ### Read in the Gaussian centroids determined for DC side illumination using dat2CTD.py
        with open(DC_param_file) as file:
            for line in file:
                if '#' not in line:
                    splitline = line.split(', ')
                    p = int(splitline[0])
                    n = int(splitline[1])
                    self.DC_params[p][n] = float(splitline[2])

        ### Calculate the slope or the stretch factor mapping from simulated ctd to measured ctd
        self.slope = ((self.AC_params - self.DC_params)/(AC_sim_mean - DC_sim_mean)) * (self.AC_params!=0.0) * (self.DC_params!=0.0)

        ### Calculate the intercept or offset mapping from simulated ctd to measured ctd
        self.intercept = (((self.AC_params+self.DC_params) - self.slope*(AC_sim_mean + DC_sim_mean))/2.) * (self.AC_params!=0.0) * (self.DC_params!=0.0)
        
        ### Some pixels were unable to be fit due to low counts. Take the mean of adjacent good pixels to guess the slope
        slope_buff = np.array([[0.0 for p in range(39)] for n in range(39)])
        slope_buff[1:-1, 1:-1] = self.slope
        for p in range(37):
            for n in range(37):
                if self.slope[p][n]==0.0:
                    buff_block = slope_buff[p:p+3,n:n+3]
                    mask = (buff_block!=0.0)
                    self.slope[p][n] = np.mean(buff_block[mask])

        ### Some pixels were unable to be fit due to low counts. Take the mean of adjacent good pixels to guess the intercept
        intercept_buff = np.array([[0.0 for p in range(39)] for n in range(39)])
        intercept_buff[1:-1, 1:-1] = self.intercept
        for p in range(37):
            for n in range(37):
                if self.intercept[p][n]==0.0:
                    buff_block = intercept_buff[p:p+3,n:n+3]
                    mask = (buff_block!=0.0)
                    self.intercept[p][n] = np.mean(buff_block[mask])

        ### Save a depth calibration file if a path is provided.
        ### Currently, depth calibration with nuclearizer doesn't work for the spare detector, so don't use this file yet.
        if savefile is not None:
            pix_code = []
            pix_stretch = []
            pix_offset = []
            for p in range(37):
                for n in range(37):
                    ### p=AC=y strips I'm pretty sure...
                    pix_code.append(int(110000 + 100*(n+1) + p+1))
                    pix_stretch.append(self.slope[p][n])
                    pix_offset.append(self.intercept[p][n])
            pix_code = np.array(pix_code)
            pix_stretch = np.array(pix_stretch)[np.argsort(pix_code)]
            pix_offset = np.array(pix_offset)[np.argsort(pix_code)]
            pix_code = np.sort(pix_code)
            with open(savefile, 'w') as file:
                for i in range(len(pix_code)):
                    file.write(str(pix_code[i]) + '   ' + str(pix_stretch[i]) + '   ' + str(pix_offset[i]) + '   0.035   1.0\n')
        
        ### read in the simulated CTD -> depth file.     
        sim_ctd = []
        sim_depth = []
        with open(sim_file) as file:
            for line in file:
                sim_ctd.append(float(line.split(',')[1]))
                sim_depth.append(float(line.split(',')[0]))
        self.sim_ctd = np.array(sim_ctd)
        self.sim_depth = np.array(sim_depth)

        ### Because CTD in not monotonic, there are regions on the edges of the detector where we can't interpolate.
        ### Determine these regions and don't interpolate there.
        sim_ctd_diff = self.sim_ctd[1:] - self.sim_ctd[:-1]
        diff_midpoint = int(len(sim_ctd_diff)/2)
        # print(diff_midpoint)
        
        ### Determine the first index where monotonicity begins
        self.start = 1
        while sim_ctd_diff[:diff_midpoint][-self.start]!=0.0:
            self.start+=1
        self.start = diff_midpoint-self.start + 1
        
        ### Determine the last index to consider before monotonicity stops
        self.end = 0
        while sim_ctd_diff[diff_midpoint:][self.end]!=0.0:
            self.end+=1
        self.end = diff_midpoint+self.end
        
        # print(start)
        # print(end)
        # sim_interp = interp1d(sim_ctd[start:end+1], sim_depth[start:end+1], bounds_error=False)
        self.sim_interp = CubicSpline(np.flip(self.sim_ctd[self.start:self.end+1]), np.flip(self.sim_depth[self.start:self.end+1]), extrapolate=True)

    ### Map the p and n times to a depth. Also returns a stretched and offset ctd
    def depth_from_timing(self, p_strip, n_strip, p_time, n_time):
        p_time = np.array(p_time)
        n_time = np.array(n_time)

        ### map the observed ctd to the simulated curve using the stretch and offset for the specified pixel
        ctd_obs = (p_time + (5.*np.random.rand(*p_time.shape))) - (n_time + (5.*np.random.rand(*n_time.shape)))
        ### Note that strip numbers that come out of nuclearizer are 1-indexed
        ctd_stretch = (ctd_obs-self.intercept[p_strip-1][n_strip-1])/self.slope[p_strip-1][n_strip-1]

        ### Map the simulated ctd to depth
        depth = self.sim_interp(ctd_stretch)

        ### ctds that fall outside the range of monotonicity will be assigned a depth corresponding to 
        ### the center of the region where monotonicity is not maintained.
        edge_high_mask = ctd_stretch < self.sim_ctd[self.end]
        edge_low_mask = ctd_stretch > self.sim_ctd[self.start]
        depth[edge_low_mask] = np.mean(self.sim_depth[:self.start])
        depth[edge_high_mask] = np.mean(self.sim_depth[self.end+1:])
        return ctd_obs, ctd_stretch, depth

    def get_simdata(self):
        return self.sim_depth, self.sim_ctd


def make_df_from_dat(files, e_min = 640., e_max = 672.):
    
    # SH det  side  strip_num (starting at 1)  did strip trigger (0 or 1)  
           #timing raw AD counter units  corrected AD counter units  ADC  energy  ??  
    # HT x y z energy [cm, cm, cm, keV]   

    rows = []
    col_names = ["ID","det","strip_p","energy_p", "time_p", "strip_n","energy_n", "time_n", "x","y","z"]
    for file in files:
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

                            time_p = float(ev_block[3].split(" ")[5])
                            time_n = float(ev_block[4].split(" ")[5])

                            # select photopeak events
                            if (energy_p < e_max and energy_p > e_min) and (energy_n < e_max and energy_n > e_min) and (np.abs(time_p)<400) and (np.abs(time_n)<400):

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
                                columns = [ID,det,strip_p,energy_p,time_p,strip_n,energy_n,time_n,x,y,z]
                                rows.append(columns)
                    ev_block = []
                    
                else:
                    ev_block.append(line)
    
    df = pd.DataFrame(rows,columns=col_names)
    return df

# With updated nuclearizer version, the function below should no longer be needed. 

# def make_df_from_dat_CC(files, e_min = 640., e_max = 672.):
    
#     # SH det  side  strip_num (starting at 1)  did strip trigger (0 or 1)  
#            #timing raw AD counter units  corrected AD counter units  ADC  energy  ??  
#     # HT x y z energy [cm, cm, cm, keV]   

#     rows = []
#     col_names = ["ID","det","strip_p","energy_p","strip_n","energy_n"]

#     for file in files:
#         with open(file, "r") as f:
            
#             ev_block = []
            
#             for line in f:
#             #for each line, start a block of lines corresponding to an event
                
#                 if line.startswith('SE'):
#                 # If the accumulated block has 8 lines then it's a single-pixel event
                    
#                     if (len(ev_block) ==8 and np.product(np.array([("StripPairingIncomplete" not in l) for l in ev_block]))) or (len(ev_block)==9):
                        
#                         ID = int(ev_block[0].split(" ")[-1].strip("\n"))

#                         if len(ev_block[3].split(" ")) < 9:
#                                print(ev_block)
#                         if len(ev_block[4].split(" ")) < 9:
#                                print(ev_block)

#                         energy_p = float(ev_block[3].split(" ")[8])
#                         energy_n = float(ev_block[4].split(" ")[8])

#                         # select photopeak events
#                         if (energy_p < e_max and energy_p > e_min) and (energy_n < e_max and energy_n > e_min):

#                             # save info from SH p line
#                             det = int(ev_block[3].split(" ")[1])
#                             strip_p = int(ev_block[3].split(" ")[3])

#                             # save info from SH n line
#                             strip_n = int(ev_block[4].split(" ")[3])

#                             # save position [cm] info from HT line
#                             # x = float(ev_block[5].split(" ")[1])
#                             # y = float(ev_block[5].split(" ")[2])
#                             # z = float(ev_block[5].split(" ")[3])

#                             # save info to df
#                             columns = [ID,det,strip_p,energy_p,strip_n,energy_n]
#                             rows.append(columns)
#                     ev_block = []
                    
#                 else:
#                     ev_block.append(line)
        
#     df = pd.DataFrame(rows,columns=col_names)
#     return df