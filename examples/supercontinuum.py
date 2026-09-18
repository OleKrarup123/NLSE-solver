import sys
sys.path.append('../') #Makes sure that this script can find the ssfm code. 
from NLSE.ssfm_functions import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['figure.dpi'] = 200
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False
rcParams['lines.linewidth'] = 3







def SC_paper_1():
    """
    In this example, we will re-create the results from the following paper:
    https://www.sciencedirect.com/science/article/pii/S2211379720317228 

    The input signal and the fiber used in the paper have been built in to ssfm_functions.py
    """
    np.random.seed(123)




    input_signal=input_signal_from_SC_paper #Load built-in input signal from the paper

    fiber=fiber_from_SC_paper #Load built-in fiber from the paper

    fiber.plot_beta2_versus_freq(freq_min_Hz=-30e12,freq_max_Hz=50e12)

    #Run simulation
    ssfm_result_list = SSFM(
        fiber=fiber,
        input_signal=input_signal,
        show_progress_flag=True,
        experiment_name='supercontinuum'
    )

    dB_cutoff = -40

    plot_everything_about_result(
        ssfm_result_list,
        dB_cutoff_pulse=dB_cutoff,
        nrange_pulse=1400,
        dB_cutoff_spectrum=dB_cutoff,
        nrange_spectrum=1300,#1200,
        show_3D_plot_flag=False)



def SC_paper_2():
    """
    In this example, we will re-create the results from the following paper:
    https://www.researchgate.net/publication/404127465_Octave-spanning_supercontinuum_generation_in_a_wafer-scale_low_loss_deuterated_silicon_nitride_waveguide 

    It simulates the propagation of high power pulses with durations of around 200fs through a custom-made on-chip waveguide. 

    """
    np.random.seed(123)
    ################ Set up time axis of simulation ################ 
    N = 2 ** 19  # Number of points on the time axis
    dt = 1e-15  # Time resolution [s]. The paper involves sech pulses with a duration of 227fs, so the time resolution should be small compared to this value. 
    center_wavelength_m = 1555e-9
    print(wavelength_to_freq(500e-9)/1e12)
    center_freq = wavelength_to_freq(center_wavelength_m)
    time_freq = TimeFreq(number_of_points=N,
                                time_step_s=dt,
                                center_frequency_Hz=center_freq)

    #Set up the initial input signal as a Gaussian pulse.
    input_signal= InputSignal(time_freq=time_freq,
                              duration_s=227e-15,
                              amplitude_sqrt_W=np.sqrt(330), #Paper simulates both 330W and 1606W pulses
                              pulse_type="sech",
                              FFT_tol=1e-4,
                              noise_stdev_sqrt_W=np.sqrt(0e-3),
                              describe_input_signal_flag=False)

    fiber = FiberSpan(length_m=1.21e-2, #5.21cm
                      beta_list=[-5.35e-26,-3.5e-41,2.5e-55], #[s^2/m,s^3/m,...]  s^(entry+2)/m
                        alpha_dB_per_m=54, #0.54dB/cm
                        gamma_per_W_per_m=0.924,
                        number_of_steps=2**8,
                        use_self_steepening_flag=True)  

    ssfm_result_list= SSFM(fiber=fiber,input_signal=input_signal,show_progress_flag=True)
    nrange=N//2
    dB_cutoff=-70
    plot_first_and_last_spectrum(ssfm_result_list,nrange=nrange,dB_cutoff=dB_cutoff)
    #plot_spectrum_matrix_2D(ssfm_result_list,nrange=nrange,dB_cutoff=dB_cutoff)

if __name__ == "__main__":
    SC_paper_1()
    #SC_paper_2() #Work in progress!