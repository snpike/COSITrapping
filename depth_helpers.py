#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numba_stats import norm
from scipy import stats
from scipy.integrate import quad
import pandas as pd
from scipy.interpolate import interp1d, CubicSpline, UnivariateSpline
import dat2CTD
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from iminuit import cost, Minuit
import seaborn as sns

sns.set_context('talk', font_scale=1.0)
sns.set_palette('colorblind')

### Line profile parameters determined via simulations.
global_gamma = 0.50
global_CoverB = 0.13
global_D = 0.028
global_sigma_ratio = 0.85

### Important calibration lines.
source_dict = {'Am-241': 59.5409, 'Cs-137': 661.657,'Co-57': 122.06065, 'Ba-133': 356.0129}


### Functions for producing line profiles.

def threshold(x, x0, sigma, Eth):
    # return (1+stats.norm.cdf(x+Eth, x0, sigma))
    return stats.norm.cdf(x+Eth, x0, sigma)

def shelf(x, x0, sigma):
    return (1.-norm.cdf(x, x0, sigma))

def shelf_scipy(x, x0, sigma):
    return (1.-stats.norm.cdf(x, x0, sigma))

def exp_tail(x, x0, gamma):
    return (np.exp(gamma*(x-x0)))

def linear_tail(x, x0, m):
    return (1+m*(x-x0))

def gaussian(x,x0,sigma):
    return np.exp(-(x-x0)**2/(2*sigma**2))

def gauss_plus_tail(x, BoverA, x0, sigma_gauss, gamma, CoverB, D, sigma_ratio):
    return (gaussian(x,x0,sigma_gauss) + BoverA*exp_tail(x, x0, gamma)*shelf(x,x0, sigma_gauss*sigma_ratio) + \
            BoverA*CoverB*linear_tail(x, x0, D)*shelf(x,x0, sigma_gauss*sigma_ratio))

### We need a normalized version in order to perform unbinned fitting.
def gauss_plus_tail_pdf(x, BoverA, x0, sigma_gauss, gamma, CoverB, D, sigma_ratio, Emin = 640., Emax = 672.):
    return gauss_plus_tail(x, BoverA, x0, sigma_gauss, gamma, CoverB, D, sigma_ratio)/\
    quad(gauss_plus_tail, Emin, Emax, args=(BoverA, x0, sigma_gauss, gamma, CoverB, D, sigma_ratio))[0]


class DepthCalibrator_Am241:
    def __init__(self, AC_param_file, DC_param_file, AC_sim_ev, DC_sim_ev, sim_file, savefile=None):
        ### take in the Gaussian means from AC- and DC- side illumination with Am241 and compare to simulations to return a CTD->depth mapping function.
        self.DC_sim_CTD = []
        self.AC_sim_CTD = []
        self.AC_params = np.array([[[0.0,0.0] for p in range(37)] for n in range(37)])
        self.DC_params = np.array([[[0.0,0.0] for p in range(37)] for n in range(37)])
        self.AC_noise = 0.0
        self.DC_noise = 0.0

        ### Read in the Gaussian centroids determined for AC side illumination using dat2CTD.py
        with open(AC_param_file) as file:
            for line in file:
                if '#' not in line:
                    splitline = line.split(', ')
                    p = int(splitline[0])
                    n = int(splitline[1])
                    self.AC_params[p][n] = [float(splitline[2]), float(splitline[3])]

        ### Read in the Gaussian centroids determined for DC side illumination using dat2CTD.py
        with open(DC_param_file) as file:
            for line in file:
                if '#' not in line:
                    splitline = line.split(', ')
                    p = int(splitline[0])
                    n = int(splitline[1])
                    self.DC_params[p][n] = [float(splitline[2]), float(splitline[3])]


        with open(AC_sim_ev) as file:
            for line in file:
                self.AC_sim_CTD.append(float(line.split(',')[1]))
        self.AC_sim_CTD = np.array(self.AC_sim_CTD) 

        with open(DC_sim_ev) as file:
            for line in file:
                self.DC_sim_CTD.append(float(line.split(',')[1]))
        self.DC_sim_CTD = np.array(self.DC_sim_CTD)

        ### For each pixel, make histograms of the simulated CTD adding successively more noise
        ### Fill in the AC_params and DC_params accordingly
        noise_mu_sigma_AC = []
        for i in range(100):
            noise=np.random.normal(loc=0.0, scale=float(i)+10., size=self.AC_sim_CTD.shape)
            temp_sim = self.AC_sim_CTD + noise
            sim_hist, sim_bin_edges = np.histogram(temp_sim, bins=700, range=(-350.,350.), density=True)
            bin_centers = (sim_bin_edges[1:] + sim_bin_edges[:-1])/2.
            p0 = [500., -200., 30.]
            popt, pcov = curve_fit(dat2CTD.gauss, bin_centers, sim_hist, p0=p0)
            noise_mu_sigma_AC.append([(float(i)+10.), popt[1], popt[2]])
        noise_mu_sigma_AC = np.array(noise_mu_sigma_AC)

        noise_mu_sigma_DC = []
        for i in range(100):
            noise=np.random.normal(loc=0.0, scale=float(i)+10., size=self.DC_sim_CTD.shape)
            temp_sim = self.DC_sim_CTD + noise
            sim_hist, sim_bin_edges = np.histogram(temp_sim, bins=700, range=(-350.,350.), density=True)
            bin_centers = (sim_bin_edges[1:] + sim_bin_edges[:-1])/2.
            p0 = [500., 125., 30.]
            popt, pcov = curve_fit(dat2CTD.gauss, bin_centers, sim_hist, p0=p0)
            noise_mu_sigma_DC.append([(float(i)+10.), popt[1], popt[2]])
        noise_mu_sigma_DC = np.array(noise_mu_sigma_DC)

        self.slope = np.zeros((37, 37))
        self.intercept = np.zeros((37,37))

        AC_noise_match = []
        DC_noise_match = []
        for p in range(37):
            for n in range(37):
                AC_mu_obs = self.AC_params[p][n][0]
                DC_mu_obs = self.DC_params[p][n][0]
                AC_sigma_obs = self.AC_params[p][n][1]
                DC_sigma_obs = self.DC_params[p][n][1]
                AC_mu_sim = noise_mu_sigma_AC.T[1][np.argmin(np.abs(noise_mu_sigma_AC.T[2]-AC_sigma_obs))]
                DC_mu_sim = noise_mu_sigma_DC.T[1][np.argmin(np.abs(noise_mu_sigma_DC.T[2]-DC_sigma_obs))]

                if AC_mu_obs!=0 and DC_mu_obs!=0:
                    AC_noise_match.append(noise_mu_sigma_AC.T[0][np.argmin(np.abs(noise_mu_sigma_AC.T[2]-AC_sigma_obs))])
                    DC_noise_match.append(noise_mu_sigma_DC.T[0][np.argmin(np.abs(noise_mu_sigma_DC.T[2]-DC_sigma_obs))])

                self.slope[p][n] = ((AC_mu_obs - DC_mu_obs)/(AC_mu_sim - DC_mu_sim)) * (AC_mu_obs!=0.0) * (DC_mu_obs!=0.0)
                self.intercept[p][n] = (((AC_mu_obs+DC_mu_obs) - self.slope[p][n]*(AC_mu_sim + DC_mu_sim))/2.) * (AC_mu_obs!=0.0) * (DC_mu_obs!=0.0)

        ### Save a noise parameter for the AC and DC sides which corresponds to the mean over all pixels
        ### of the noise required to make the simulations match the observed Am241 CTDs.
        print(str(round(np.mean(AC_noise_match), 1)) + ' +/- ' + str(round(np.std(AC_noise_match), 1)))
        print(str(round(np.mean(DC_noise_match), 1)) + ' +/- ' + str(round(np.std(DC_noise_match), 1)))
        self.AC_noise = np.mean(AC_noise_match)
        self.DC_noise = np.mean(DC_noise_match)
        
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

    def depth_from_timing(self, p_strip, n_strip, p_time, n_time):
        ### Map the p and n times to a depth. Also returns a stretched and offset ctd
        ### This function is superseded by depth_from_timing_prob.
        p_time = np.array(p_time)
        n_time = np.array(n_time)

        ### map the observed ctd to the simulated curve using the stretch and offset for the specified pixel
        ctd_obs = (p_time + (5.*np.random.rand(*p_time.shape))) - (n_time + (5.*np.random.rand(*n_time.shape)))
        ### Note that strip numbers that come out of nuclearizer are 1-indexed
        ctd_stretch = (ctd_obs-self.intercept[p_strip-1][n_strip-1])/(self.slope[p_strip-1][n_strip-1])

        ### Map the simulated ctd to depth
        depth = self.sim_interp(ctd_stretch)

        ### ctds that fall outside the range of monotonicity will be assigned a depth corresponding to 
        ### the center of the region where monotonicity is not maintained.
        edge_high_mask = ctd_stretch < self.sim_ctd[self.end]
        edge_low_mask = ctd_stretch > self.sim_ctd[self.start]
        depth[edge_low_mask] = np.mean(self.sim_depth[:self.start])
        depth[edge_high_mask] = np.mean(self.sim_depth[self.end+1:])
        return ctd_obs, ctd_stretch, depth

    ### Map the p and n times to a depth. Also returns a stretched and offset ctd
    def depth_from_timing_prob(self, p_strip, n_strip, p_time, n_time, energy_p):
        ### Probabilistic depth estimate
        p_time = np.array(p_time)
        n_time = np.array(n_time)

        ### map the observed ctd to the simulated curve using the stretch and offset for the specified pixel
        ### Note that a float between 0 and 5 is added to each because of the CC readout.
        ctd_obs = (p_time + (5.*np.random.rand(*p_time.shape))) - (n_time + (5.*np.random.rand(*n_time.shape)))
        ### Strip numbers that come out of nuclearizer are 1-indexed
        ctd_stretch = (ctd_obs-self.intercept[p_strip-1][n_strip-1])/(self.slope[p_strip-1][n_strip-1])

        ### The noise on each event is energy dependent.
        ### TODO: Double check this relation. Currently just a rough estimate.
        m = (self.AC_noise - 12.)/((1./source_dict['Am-241']) - (1./source_dict['Cs-137']))
        b = self.AC_noise - (m/source_dict['Am-241'])
        noise = b + m/(energy_p.values)
        
        ### Map the simulated ctd to depth
        depth = []
        depth_err = []
        for i, tau in enumerate(ctd_stretch):
            ### Calculate the probability distribution with depth given measured tau and noise
            prob_dist = stats.norm.pdf(self.sim_ctd, loc=tau, scale=noise[i])
            mean_depth = np.sum(prob_dist*self.sim_depth)/np.sum(prob_dist)
            sigma_depth = np.sqrt(np.sum(prob_dist*np.square(self.sim_depth-mean_depth))/np.sum(prob_dist))
            
            prob_norm = (np.sum(prob_dist) - ((prob_dist[0] + prob_dist[-1])/2.))

            ### Assign the depth randomly given the probability distribution
            prob_dist = prob_dist/prob_norm
            cdf = np.concatenate([[0.0], (np.cumsum(prob_dist)-((prob_dist + prob_dist[0])/2.0))[1:]])
            depth.append(interp1d(cdf, self.sim_depth)(np.random.uniform()))
            
            ### Assign 1-sigma depth error given the probability distribution
            depth_err.append(sigma_depth)

        return ctd_obs, ctd_stretch, depth, depth_err


    def get_simdata(self):
        return self.sim_depth, self.sim_ctd



def make_df_from_dat(files, e_min = 640., e_max = 672.):
    
    # SH det  side  strip_num (starting at 1)  did strip trigger (0 or 1)  
           #timing raw AD counter units  corrected AD counter units  ADC  energy  ??  
    # HT x y z energy [cm, cm, cm, keV]   

    rows = []
    col_names = ["ID","det","strip_p","energy_p", "time_p", "strip_n","energy_n", "time_n", "x","y","z", "z_err"]
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
                                columns = [ID,det,strip_p,energy_p,time_p,strip_n,energy_n,time_n,x,y,z, 0.0]
                                rows.append(columns)
                    ev_block = []
                    
                else:
                    ev_block.append(line)
    
    df = pd.DataFrame(rows,columns=col_names)
    return df

def make_depthplot(df, plot_suffix, z_min=0.0, z_max=1.5, num_z = 30, plot_dir="/home/cosilab/CalibrationData/figures/", source='Cs-137'):
    ### Take in a depth-calibrated dataset, bin events in depth, and fit the spectra.
    ### Produce and return depth plots.
    if source not in source_dict:
        print('Source not recognized. Must be one of the following: Am-241, Cs-137, Ba-133, Co-57')
    line_e = source_dict[source]

    z_bins = np.arange(z_min,z_max + (z_max-z_min)/num_z,(z_max-z_min)/num_z)
    df["z_binned"] = pd.cut(df["z"],bins=z_bins)

    centroid_list = []
    centroid_err_list = []

    tail_cnts = []
    gauss_cnts = []
    tail_frac = []

    emin_list = []
    emax_list = []

    fig, axes = plt.subplots(figsize = (12, num_z*2.), nrows=num_z*2, ncols=2, sharex=True, sharey='row', \
                            gridspec_kw={'hspace':0, 'wspace':0, 'height_ratios':np.array([(2,1) for i in range(num_z)]).flatten()})
    for i in range(num_z):
        
        label = "{}-{} cm".format(round(z_bins[i], 1),round(z_bins[i+1], 1))

        # fit both sides
        temp_x0 = []
        temp_x0_err = []
        for j, side in enumerate(['p', 'n']):
            if side=='p':
                carrier='electron'
            else:
                carrier='hole'

            ax = axes[i*2][j]
            resid_ax = axes[(i*2)+1][j]            
            color = 'C' + str(j)
            
            energies = df.loc[df.z.le(z_bins[i+1])&df.z.gt(z_bins[i]), 'energy_'+side].values
            temp_emin = stats.mode(np.floor(energies))[0] - 20.
            temp_emax = stats.mode(np.floor(energies))[0] + 12.
            emin_list.append(temp_emin)
            emax_list.append(temp_emax)
            emask = (energies < temp_emax) * (energies > temp_emin)
            energies = energies[emask]
            hist,binedges,_ = ax.hist(energies, histtype="step",color=color,bins=100,label=carrier + " signal, " + label)
            bin_centers = np.array((binedges[:-1] + binedges[1:]) / 2)

            c = cost.UnbinnedNLL(energies, gauss_plus_tail_pdf)

            m = Minuit(c, BoverA=0.5, x0=stats.mode(np.floor(energies))[0], sigma_gauss=1.2, gamma=global_gamma, CoverB=global_CoverB, D=global_D, sigma_ratio=global_sigma_ratio, Emin=temp_emin, Emax=temp_emax)
            m.limits["x0"] = (temp_emin, temp_emax)
            m.limits["BoverA", "sigma_gauss"] = (0, None)
            m.fixed["gamma", "CoverB", "D", "sigma_ratio", "Emin", "Emax"] = True
            m.migrad()
            m.hesse()

            BoverA, x0, sigma_gauss = m.values[:3]
            A = np.sum(hist)*(bin_centers[1]-bin_centers[0])/\
                quad(gauss_plus_tail, np.min(energies), np.max(energies), args = (BoverA, x0, sigma_gauss, global_gamma, global_CoverB, global_D, global_sigma_ratio))[0]
            # print(m.values)
            # print(m.errors['x0'])
            temp_x0.append(m.values['x0'])
            temp_x0_err.append(m.errors['x0'])
            # print(c(*m.values))
            # print(np.sum(np.log(gauss_plus_tail_pdf(df_dat_det[df_dat_det["z_binned"]==z]["energy_p"].values, *m.values))))
            B = A*BoverA
            C = B*global_CoverB

            ax.plot(bin_centers,A*gauss_plus_tail(bin_centers, BoverA, x0, sigma_gauss, global_gamma, global_CoverB, global_D, global_sigma_ratio),color= color, lw=0.5)
            resid_ax.errorbar(bin_centers[hist>0], \
                               (hist[hist>0] - A*gauss_plus_tail(bin_centers[hist>0],BoverA, x0, sigma_gauss, global_gamma, global_CoverB, global_D, global_sigma_ratio))/np.sqrt(hist[hist>0]), \
                               xerr=(bin_centers[1]-bin_centers[0])/2, yerr=1., ls='', lw=0.6, color=color)
            ax.plot(bin_centers,A*gaussian(bin_centers, x0, sigma_gauss),color= color, ls='--', lw=0.5)
            ax.plot(bin_centers, B*exp_tail(bin_centers, x0, gamma=global_gamma)*shelf(bin_centers, x0, sigma_gauss*global_sigma_ratio),color= color, ls='--', lw=0.5)
            ax.plot(bin_centers, C*linear_tail(bin_centers, x0, global_D)*shelf(bin_centers, x0, sigma_gauss*global_sigma_ratio),color= color, ls='--', lw=0.5)
            ax.axvline(x0, ls='--', color='C3')
            ax.axvline(line_e, ls='--', color='C2')
            if j==0:
                ax.set_ylabel("Counts")
            ax.set_yscale('log')
            ax.set_ylim((1, 10000))
            ax.legend(loc=2, fontsize=12)
            resid_ax.axhline(0, color='red', lw=0.5)

        axes[(i*2)+1][0].set_ylabel(r'$\chi$')
        centroid_list.append(temp_x0)
        centroid_err_list.append(temp_x0_err)

    axes[-1][0].set_xlim(np.min(emin_list), np.max(emax_list))
    axes[-1][1].set_xlim(np.min(emin_list), np.max(emax_list))
    axes[-1][0].set_xlabel("Energy (keV)")
    axes[-1][1].set_xlabel("Energy (keV)")
    plt.tight_layout()

    plt.savefig(plot_dir + 'zplitspectra_uncorr_' + plot_suffix + '.pdf')
    plt.close()

    centroid_list = np.array(centroid_list)
    centroid_err_list = np.array(centroid_err_list)

    return (z_bins[:-1]+z_bins[1:])/2, [centroid_list.T[0], centroid_err_list.T[0]], [centroid_list.T[1], centroid_err_list.T[1]]



def depth_correction(df, z_list, e_trapping, h_trapping, plot_dir="/home/cosilab/CalibrationData/figures/", plot_suffix='', source='Cs-137'):
    ### Correct measured energy on an event-by-event basis using the depth plots produced by make_depthplot.
    ### Return the dataset with corrected energies in new columns.
    if source not in source_dict:
        print('Source not recognized. Must be one of the following: Am-241, Cs-137, Ba-133, Co-57')
    line_e = source_dict[source]

    splines = {'p': UnivariateSpline(z_list, e_trapping[0]), 'n': UnivariateSpline(z_list, h_trapping[0])}
    plt.figure()
    z_err = [z_list[0], *((z_list[1:] - z_list[:-1])/2.)]
    xs = np.linspace(z_list[0]-z_err[0], z_list[-1] + z_err[-1])
    plt.plot(xs, splines['p'](xs), zorder=0, color='C0', ls = '--', lw = 0.75)
    plt.plot(xs, splines['n'](xs), zorder=0, color='C1', ls = '--', lw = 0.75)
    plt.errorbar(z_list, e_trapping[0], xerr = z_err, yerr=e_trapping[1], fmt=".", label="electron signal", color='C0')
    plt.errorbar(z_list, h_trapping[0], xerr = z_err, yerr=h_trapping[1], fmt=".", label="hole signal", color='C1')
    plt.axhline(line_e, ls='--', color='C2')
    # plt.legend(loc=4)
    plt.xlabel("Detector Depth (cm)"); plt.ylabel("Centroid Energy (keV)")
    plt.tight_layout()
    plt.savefig(plot_dir + 'e_hole_trapping_' + plot_suffix + '.pdf')
    plt.close()

    fig, axes = plt.subplots(figsize = (15, 18), nrows=2, ncols=2, sharex=True, sharey='row', gridspec_kw={'hspace':0, 'wspace':0})

    emin_list = []
    emax_list = []

    for i, side in enumerate(splines):
        if side=='p':
            carrier='electron'
        else:
            carrier='hole'
        
        color = 'C'+str(i)
        
        ax = axes[0][i]
        energies = df['energy_'+side].values
        temp_emin = stats.mode(np.floor(energies))[0] - 20.
        temp_emax = stats.mode(np.floor(energies))[0] + 12.
        
        emin_list.append(temp_emin)
        emax_list.append(temp_emax)

        emask = (energies < temp_emax) * (energies > temp_emin)
        energies = energies[emask]
        hist,binedges,_ = ax.hist(energies, histtype="step", bins=100, label="Uncorrected " + carrier + " signal")
        bin_centers = np.array((binedges[:-1] + binedges[1:]) / 2)

        c = cost.UnbinnedNLL(energies, gauss_plus_tail_pdf)

        m = Minuit(c, BoverA=0.5, x0=stats.mode(np.floor(energies))[0], sigma_gauss=1.2, gamma=global_gamma, CoverB=global_CoverB, D=global_D, sigma_ratio=global_sigma_ratio, Emin=temp_emin, Emax=temp_emax)
        m.limits["x0"] = (temp_emin, temp_emax)
        m.limits["BoverA", "sigma_gauss"] = (0, None)
        m.fixed["gamma", "CoverB", "D", "sigma_ratio", "Emin", "Emax"] = True
        m.migrad()
        m.hesse()
        print(m.params)

        BoverA, x0, sigma_gauss = m.values[:3]
        A = np.sum(hist)*(bin_centers[1]-bin_centers[0])/\
            quad(gauss_plus_tail, np.min(energies), np.max(energies), args = (BoverA, x0, sigma_gauss, global_gamma, global_CoverB, global_D, global_sigma_ratio))[0]
        B = A*BoverA
        C = B*global_CoverB

        ax.plot(bin_centers,A*gauss_plus_tail(bin_centers, BoverA, x0, sigma_gauss, global_gamma, global_CoverB, global_D, global_sigma_ratio),color= color, lw=0.5)
        ax.plot(bin_centers,A*gaussian(bin_centers, x0, sigma_gauss),color= color, ls='--', lw=0.5)
        ax.plot(bin_centers, B*exp_tail(bin_centers, x0, gamma=global_gamma)*shelf(bin_centers, x0, sigma_gauss*global_sigma_ratio),color= color, ls='--', lw=0.5)
        ax.plot(bin_centers, C*linear_tail(bin_centers, x0, global_D)*shelf(bin_centers, x0, sigma_gauss*global_sigma_ratio),color= color, ls='--', lw=0.5)
        ax.axvline(x0, ls='--', color='C3')
        ax.axvline(line_e, ls='--', color='C2')
        if i==0:
            ax.set_ylabel("Counts")
        ax.set_yscale('log')
        ax.legend(loc=2, fontsize=12)

        ### Correct the measured energies according to the CCE spline
        energies = df["energy_"+side].values * line_e/splines[side](df['z'].values)
        df['depth_corrected_energy_'+side] = energies

        ax = axes[1][i]
        temp_emin = stats.mode(np.floor(energies))[0] - 20.
        temp_emax = stats.mode(np.floor(energies))[0] + 12.
        
        emin_list.append(temp_emin)
        emax_list.append(temp_emax)

        emask = (energies < temp_emax) * (energies > temp_emin)
        energies = energies[emask]
        hist,binedges,_ = ax.hist(energies, histtype="step", bins=100, label="Corrected " + carrier + " signal")
        bin_centers = np.array((binedges[:-1] + binedges[1:]) / 2)

        c = cost.UnbinnedNLL(energies, gauss_plus_tail_pdf)

        m = Minuit(c, BoverA=0.5, x0=stats.mode(np.floor(energies))[0], sigma_gauss=1.2, gamma=global_gamma, CoverB=global_CoverB, D=global_D, sigma_ratio=global_sigma_ratio, Emin=temp_emin, Emax=temp_emax)
        m.limits["x0"] = (temp_emin, temp_emax)
        m.limits["BoverA", "sigma_gauss"] = (0, None)
        m.fixed["gamma", "CoverB", "D", "sigma_ratio", "Emin", "Emax"] = True
        m.migrad()
        m.hesse()
        print(m.params)

        BoverA, x0, sigma_gauss = m.values[:3]
        A = np.sum(hist)*(bin_centers[1]-bin_centers[0])/\
            quad(gauss_plus_tail, np.min(energies), np.max(energies), args = (BoverA, x0, sigma_gauss, global_gamma, global_CoverB, global_D, global_sigma_ratio))[0]
        B = A*BoverA
        C = B*global_CoverB

        ax.plot(bin_centers,A*gauss_plus_tail(bin_centers, BoverA, x0, sigma_gauss, global_gamma, global_CoverB, global_D, global_sigma_ratio),color= color, lw=0.5)
        ax.plot(bin_centers,A*gaussian(bin_centers, x0, sigma_gauss),color= color, ls='--', lw=0.5)
        ax.plot(bin_centers, B*exp_tail(bin_centers, x0, gamma=global_gamma)*shelf(bin_centers, x0, sigma_gauss*global_sigma_ratio),color= color, ls='--', lw=0.5)
        ax.plot(bin_centers, C*linear_tail(bin_centers, x0, global_D)*shelf(bin_centers, x0, sigma_gauss*global_sigma_ratio),color= color, ls='--', lw=0.5)
        ax.axvline(x0, ls='--', color='C3')
        ax.axvline(line_e, ls='--', color='C2')
        if i==0:
            ax.set_ylabel("Counts")
        ax.set_yscale('log')
        ax.set_ylim(bottom=10.)
        ax.legend(loc=2, fontsize=12)

    axes[1][0].set_xlim(np.min(emin_list), np.max(emax_list))
    axes[1][1].set_xlim(np.min(emin_list), np.max(emax_list))
    axes[1][0].set_xlabel("Energy (keV)")
    axes[1][1].set_xlabel("Energy (keV)")
    plt.tight_layout()

    plt.savefig(plot_dir + 'depth_corrected_spectra_' + plot_suffix + '.pdf')
    plt.close()

    return df
