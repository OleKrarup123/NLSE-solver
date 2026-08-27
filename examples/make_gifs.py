import sys
sys.path.append('../') #Makes sure that this script can find the ssfm code. 
sys.path.append('../examples') #Makes sure that this script can find the ssfm code. 
from NLSE.ssfm_functions import *
from raman import simulate_raman
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['figure.dpi'] = 200
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False
rcParams['lines.linewidth'] = 3


if __name__ == "__main__":

    # Run the raman.py example with a small number of steps to shorten the time it takes to
    # produce the animation. 
    ssfm_result_list_list = simulate_raman(number_of_steps=2**6,show_plots=False)

    #Make animated 2x2 gif of the local chirp.
    #make_chirp_gif_2x2(ssfm_result_list_list,title_list = ["No Raman", "Agrawal", "Silica Exact", "Custom"],nrange=500,chirp_range_Hz= [-75e12, 75e12],framerate= 30)


    ## Now make animated gif of the evolution of the spectrogram of a single signal.
    
    # It's best practice to first make static plots of the first and last spectrogram while iterating on
    # plotting parameters (spectrum range, time resolution etc.)
    nrange_pulse = 500
    nrange_spectrum = 8000
    dB_cutoff = -40
    time_resolution_s =30e-15 
    plot_first_and_last_spectrogram(ssfm_result_list_list[1],nrange_pulse=nrange_pulse,nrange_spectrum=nrange_spectrum,dB_cutoff=dB_cutoff,time_resolution_s=time_resolution_s)

    # When you are satisfied with your plotting parameters, you can run the following line to produce an animated spectrogram gif.
    # Note that this may take several minutes!  
    #make_spectrogram_gif(ssfm_result_list_list[1],nrange_pulse=nrange_pulse,nrange_spectrum=nrange_spectrum,dB_cutoff=dB_cutoff,time_resolution_s=time_resolution_s,framerate= 30)