import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


import define_common_blks

def diffusion_scale(t, hole=False):
    # t: drift time in nanoseconds
    # hole: boolean flag indicating whether it's for holes or electrons
    if hole:
        sigma = np.sqrt(2.0 * 29.0 * t)  # micrometers
    else:
        sigma = np.sqrt(2.0 * 25.0 * t)  # micrometers
    return sigma

def repulsion_scale_extended(t, energy, R0, hole=False):
    # t: drift time in nanoseconds
    # energy: interaction energy in keV
    # R0: radius of initial charge cloud in micrometers
    # hole: boolean flag indicating whether it's for holes or electrons
    ncharge = energy * 1000. / 2.96  # number of charge carriers
    if hole:
        eta = (R0**3 + ncharge * 1.13 * t)**(1. / 3.)  # micrometers
    else:
        eta = (R0**3 + ncharge * 0.97 * t)**(1. / 3.)  # micrometers
    return eta

def charge_distribution_extended(x, t, energy, R0, hole=False):
    # x: spatial variable
    # t: drift time in nanoseconds
    # energy: interaction energy in keV
    # R0: radius of initial charge cloud in micrometers
    # hole: boolean flag indicating whether it's for holes or electrons
    if hole:
        sigma = diffusion_scale(t, hole=True)  # diffusion size scale, micrometers
        eta = repulsion_scale_extended(t, energy, R0, hole=True)  # repulsion size scale, micrometers
    else:
        sigma = diffusion_scale(t, hole=False)
        eta = repulsion_scale_extended(t, energy, R0, hole=False)

    a = x - eta  # lower limit of integration
    b = x + eta  # upper limit of integration

    density = (0.75 / (eta**3)) * (0.5) * (eta**2 - x**2 - sigma**2) * \
              (erf(b / (np.sqrt(2.0) * sigma)) - erf(a / (np.sqrt(2.0) * sigma)))
    density += (0.75 / (eta**3)) * (sigma / np.sqrt(2.0 * np.pi)) * \
               ((2.0 * x - a) * np.exp(-0.5 * a**2 / sigma**2) - (2.0 * x - b) * np.exp(-0.5 * b**2 / sigma**2))

    return density

def photopeak(x, a):
    """Fit a Gaussian function to the data."""
    f = a[0] * np.exp(-0.5 * ((x - a[1]) / a[2])**2)
    return f

def tailing_unconstrained(x, b):
    """
    Compute the tailing unconstrained function.

    Parameters:
        x (array-like): Input data.
        b (array-like): Parameters for the function.

    Returns:
        array-like: Evaluated function values.
    """
    f = (b[0] * np.exp(b[3] * (x - b[1])) + b[4]) * \
        (1.0 - erf((x - b[1]) / (b[2] * np.sqrt(2.)))) * \
        (1.0 + erf((x - b[1] + 18.0) / (b[2] * np.sqrt(2.))))
    return f

def tailing3(x, b):
    """
    Compute the tailing3 function.

    Parameters:
        x (array-like): Input data.
        b (array-like): Parameters for the function.

    Returns:
        array-like: Evaluated function values.
    """
    f = b[0] * (np.exp(0.47 * (x - b[1])) + 0.079) * \
        (1.0 - erf((x - b[1]) / (b[2] * np.sqrt(2.)))) * \
        (1.0 + erf((x - b[1] + 18.0) / (b[2] * np.sqrt(2.))))
    return f

def tailing_lab(x, b):
    """
    Compute the tailing_lab function.

    Parameters:
        x (array-like): Input data.
        b (array-like): Parameters for the function.

    Returns:
        array-like: Evaluated function values.
    """
    f = b[0] * np.exp(b[3] * (x - b[1])) + b[4] * (1.0 - erf((x - b[1]) / (b[2] * np.sqrt(2.)))) * \
        (1.0 + erf((x - b[1] + 18.0) / (b[2] * np.sqrt(2.))))
    return f

def tailing_fixedpeak(x, b):
    """
    Compute the tailing_fixedpeak function.

    Parameters:
        x (array-like): Input data.
        b (array-like): Parameters for the function.

    Returns:
        array-like: Evaluated function values.
    """
    f = b[0] * np.exp(b[1] * (x - photon_energy) + b[2]) * \
        (1.0 - erf((x - photon_energy) / (sigma_energy * np.sqrt(2.)))) * \
        (1.0 + erf((x - photon_energy + strip_thresh) / (sigma_energy * np.sqrt(2.))))
    return f

def tailing_fixedpeak_v3(x, b):
    """
    Compute the tailing_fixedpeak_v3 function.

    Parameters:
        x (array-like): Input data.
        b (array-like): Parameters for the function.

    Returns:
        array-like: Evaluated function values.
    """
    f = b[0] * np.exp(b[1] * (x - photon_energy)) * \
        (1.0 - erf((x - photon_energy) / (b[2] * np.sqrt(2.)))) * \
        (1.0 + erf((x - photon_energy + strip_thresh) / (b[2] * np.sqrt(2.))))
    return f

def tailing_unconstrained_v3(x, a):
    """
    Compute the tailing_unconstrained_v3 function.

    Parameters:
        x (array-like): Input data.
        a (array-like): Parameters for the function.

    Returns:
        array-like: Evaluated function values.
    """
    f = a[0] * np.exp(a[2] * (x - a[1])) * \
        (1.0 - erf((x - a[1]) / (a[3] * np.sqrt(2.)))) * \
        (1.0 + erf((x - a[1] + strip_thresh) / (a[3] * np.sqrt(2.))))
    return f

def tailing_fixedpeak_v4(x, b):
    """
    Compute the tailing_fixedpeak_v4 function.

    Parameters:
        x (array-like): Input data.
        b (array-like): Parameters for the function.

    Returns:
        array-like: Evaluated function values.
    """
    f = (b[0] * np.exp(b[1] * (x - photon_energy)) + b[2]) * \
        (1.0 - erf((x - photon_energy) / (b[3] * np.sqrt(2.)))) * \
        (1.0 + erf((x - photon_energy + strip_thresh) / (b[3] * np.sqrt(2.))))
    return f

def tailing_unconstrained_v4(x, a):
    """
    Compute the tailing_unconstrained_v4 function.

    Parameters:
        x (array-like): Input data.
        a (array-like): Parameters for the function.

    Returns:
        array-like: Evaluated function values.
    """
    f = (a[0] * np.exp(a[2] * (x - a[1])) + a[3]) * \
        (1.0 - erf((x - a[1]) / (a[4] * np.sqrt(2.)))) * \
        (1.0 + erf((x - a[1] + strip_thresh) / (a[4] * np.sqrt(2.))))
    return f


def tailpeak_v6(x, params):
    # x: spatial variable
    # params: parameters for fitting
    a0, a1, a2, a3, a4 = params
    f = a0 * np.exp(-0.5 * ((x - a1) / a2)**2) + \
        a3 * (np.exp(0.49 * (x - a1)) + 0.13 * (1. + 0.028 * (x - a1))) * \
        (1.0 - erf((x - a1) / (0.84 * a2 * np.sqrt(2.)))) * \
        (1.0 + erf((x - a1 + strip_thresh) / (0.84 * a2 * np.sqrt(2.))))
    return f

def tailing_v6(x, params):
    # x: spatial variable
    # params: parameters for fitting
    b0, b1, b2, b3 = params
    f = b0 * (np.exp(0.49 * (x - b1)) + 0.13 * (1. + 0.028 * (x - b1))) * \
        (1.0 - erf((x - b1) / (0.84 * b2 * np.sqrt(2.)))) * \
        (1.0 + erf((x - b1 + strip_thresh) / (0.84 * b2 * np.sqrt(2.))))
    return f

def tailing_fixedpeak_v6(x, b):
    """
    Fixed E_0 tailing function with variable sigma_tail,
    short exponential + step function.
    
    Parameters:
        x: array-like
            Independent variable (energy axis).
        b: array-like
            Parameter vector containing [amplitude_exp, decay_rate_exp,
            amplitude_step, slope_step, sigma_tail].
    
    Returns:
        f: array-like
            Evaluated function values.
    """
    # Define the exponential decay term
    exp_term = b[0] * np.exp(b[1] * (x - photon_energy))
    
    # Define the step function term
    step_term = b[2] * (1. + b[3] * (x - photon_energy))
    
    # Define the Gaussian tail terms
    left_tail = 1.0 - erf((x - photon_energy) / (b[4] * np.sqrt(2.)))
    right_tail = 1.0 + erf((x - photon_energy + strip_thresh) / (b[4] * np.sqrt(2.)))
    
    # Compute the function as a combination of the terms
    f = (exp_term + step_term) * left_tail * right_tail
    
    return f


def tail_fit_extended(energy = 661.66, drift_time=250., extended=False)
    # Define global variables
    photon_energy = energy  # keV, photopeak energy
    
    practical_range = (0.8) * (1.e1) * (0.55) * photon_energy * (1. - 0.9841 / (1. + 0.0030 * photon_energy)) / (5.323)
    R0 = 0.0  # range
    tau = drift_time  # drift time, nanoseconds
    strip_thresh = 18.0  # keV, energy threshold for triggering neighbor strip
    sigma_energy = (2.17 + 0.65 * (photon_energy / 1000.) ** 0.5) / 2.35  # keV, intrisic 1-sigma spectral resolution at photopeak energy
    half_pitch = 2000. / 2.0  # micrometers, strip half pitch
    emin = int(photon_energy - 22.)
    emax = emin + 30.
    dechn = 0.2  # keV
    nchn = int((emax - emin) / dechn)
    c1 = 8  # start channel for fitting
    c2 = 136  # end channel for fitting
    spce = np.zeros(nchn)  # electron photopeak spectrum
    losse = np.zeros(nchn)  # electron charge sharing spectrum
    spch = np.zeros(nchn)  # hole photopeak spectrum
    lossh = np.zeros(nchn)  # hole charge sharing spectrum
    echn = np.arange(emin, emax, dechn) + dechn / 2.0
    ndx = 100000  # number of spatial points sampled
    dx = half_pitch / ndx  # micrometers
    x = np.arange(ndx) * dx

    ce=charge_distribution_extended(x,tau,photon_energy,R0,hole=False)


    for i in range(ndx):
        # Compute the energy loss for electrons
        eloss = np.sum(ce[i:]) * dx * photon_energy
        
        # Check if energy loss is below the strip threshold
        if eloss < strip_thresh:
            # Generate a random energy loss within the Gaussian distribution
            nrge = photon_energy - eloss + sigma_energy * np.random.randn()
            
            # Determine the channel number based on energy
            chne = 0
            if nrge >= emin:
                chne = int((nrge - emin) / dechn)
                
                # Update the appropriate spectrum
                if eloss > 0.1:  # Charge loss spectrum
                    losse[chne] += 1
                else:  # Photopeak spectrum
                    spce[chne] += 1


    # Define fitting parameters and weights
    weights = np.ones(nchn)
    b_guess = [max(losse) / 4., 0.4, max(losse) / 40., 0.04, 1.0]  # Initial guess for fitting parameters
    sigb = np.ones_like(b_guess)  # Sigma for fitting parameters

    # Perform fitting for electron signal
    popt, pcov = curve_fit(tailing_fixedpeak_v6, echn[c1:c2], losse[c1:c2], p0=b_guess, sigma=sigb)

    # Plot electron signal
    plt.plot(echn, losse, label='Electron Signal')
    plt.plot(echn, spce, label='Electron Photopeak Spectrum')
    plt.plot(echn, spce + losse, label='Combined Electron Spectrum')
    plt.axvline(x=photon_energy, linestyle='--', color='k', label='Photopeak Energy')
    plt.xlabel('Energy [keV]')
    plt.ylabel('Counts')
    plt.title('Electron Signal')
    plt.legend()
    plt.show()

    # Print fitting results
    print('Electron Fit Parameters:')
    print('b1:', popt[1])
    print('b2/b0:', popt[2] / popt[0])
    print('b3:', popt[3])
    print('Total Counts:', np.sum(losse))
    print('Total Fitted Counts:', np.sum(tailing_fixedpeak_v6(echn[c1:c2], *popt)))
    print('Resolution Ratio (b4/sigma_energy):', popt[4] / sigma_energy)
