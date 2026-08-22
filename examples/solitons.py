import sys
sys.path.append('../') #Makes sure that this script can find the ssfm code. 
from ssfm_functions import ssfm_functions as sf
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['figure.dpi'] = 200
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False
rcParams['lines.linewidth'] = 3


def simulate_solitons():
    """
    Simulates the propagation of solitons through a fiber with anormalous dispersion.

    """
    ################ Set up time axis of simulation ################ 
    N = 2 ** 15  # Number of points on the time axis
    dt = 10e-15  # Time resolution [s]
    center_wavelength_m = 1550e-9

    center_freq = sf.wavelength_to_freq(center_wavelength_m)
    time_freq = sf.TimeFreq(number_of_points=N,
                                time_step_s=dt,
                                center_frequency_Hz=center_freq)


    #Choose soliton duration here since we need it for calculating the characteristic length of the fiber and the 
    #characteristic amplitude of the soliton
    duration_s = 100*time_freq.time_step_s


    ################ Set up fibers ################ 
    alpha_dB_per_km = 0 #dB/m
    beta_list = [sf.BETA2_AT_1550_NM_TYPICAL_SMF_S2_PER_M]   #[s^2/m,s^3/m,...]  s^(entry+2)/m
    
    gamma_W_per_m =  1e-3 # 1/W/m

    length_m = np.pi/2*duration_s**2/np.abs(beta_list[0]) #m
    number_of_steps = 2**10

    fiber = sf.FiberSpan(
        length_m,
        number_of_steps,
        gamma_W_per_m,
        beta_list,
        alpha_dB_per_km)
    fiber_link = sf.FiberLink([fiber]) #The SSFM solver needs a "FiberLink" class, which is generated from a list of 1 or more fibers. 

    #Set up identical fiber with no nonlinearity for comparison
    fiber_no_NL = sf.FiberSpan(
        length_m,
        number_of_steps,
        0,
        beta_list,
        alpha_dB_per_km)
    fiber_link_no_NL = sf.FiberLink([fiber_no_NL]) #The SSFM solver needs a "FiberLink" class, which is generated from a list of 1 or more fibers. 


    ################ Set up input signals ################ 
    #Set up input signal for an N=2 soliton 
    soliton_order_2 = 2
    A_char_sqrt_W_order_2 =soliton_order_2*np.sqrt(np.abs(fiber.beta_list[0])/fiber.gamma_per_W_per_m)/duration_s  #Carefuly chosen peak field to ensure soliton propagation

    input_signal_order_2= sf.InputSignal(time_freq,
                   amplitude_sqrt_W=A_char_sqrt_W_order_2,
                   duration_s=duration_s,
                   pulse_type="sech",
                   describe_input_signal_flag=False)



    #Set up input signal for an N=4 soliton
    soliton_order_4 = 4
    A_char_sqrt_W_order_4 =soliton_order_4*np.sqrt(np.abs(fiber.beta_list[0])/fiber.gamma_per_W_per_m)/duration_s  #Carefuly chosen peak field to ensure soliton propagation

    input_signal_order_4= sf.InputSignal(time_freq,
                   amplitude_sqrt_W=A_char_sqrt_W_order_4,
                   duration_s=duration_s,
                   pulse_type="sech",
                   describe_input_signal_flag=False)


    ################ Run split step simulations ################ 
    ssfm_result_list_order_2_no_NL = sf.SSFM(fiber_link=fiber_link_no_NL,
                               input_signal=input_signal_order_2,
                               show_progress_flag=True)

    ssfm_result_list_order_2 = sf.SSFM(fiber_link=fiber_link,
                               input_signal=input_signal_order_2,
                               show_progress_flag=True)

    ssfm_result_list_order_4 = sf.SSFM(fiber_link=fiber_link,
                               input_signal=input_signal_order_4,
                               show_progress_flag=True)


    ################ Make plots of results ################ 
    nrange_pulse = 1000
    dB_cutoff_pulse = -40
    nrange_spectrum = 1200
    dB_cutoff_spectrum = -40


    #Without nonlinearity, the initial sech pulse will simply broaden in the time domain
    sf.plot_everything_about_result(ssfm_result_list_order_2_no_NL,
                                    nrange_pulse=nrange_pulse,
                                    dB_cutoff_pulse=dB_cutoff_pulse,
                                    nrange_spectrum=nrange_spectrum,
                                    dB_cutoff_spectrum=dB_cutoff_spectrum)

    #The N=2 soliton oscillates once
    sf.plot_everything_about_result(ssfm_result_list_order_2,
                                    nrange_pulse=nrange_pulse,
                                    dB_cutoff_pulse=dB_cutoff_pulse,
                                    nrange_spectrum=nrange_spectrum,
                                    dB_cutoff_spectrum=dB_cutoff_spectrum)

    #The N=4 soliton oscillates 3 times
    sf.plot_everything_about_result(ssfm_result_list_order_4,
                                    nrange_pulse=nrange_pulse,
                                    dB_cutoff_pulse=dB_cutoff_pulse,
                                    nrange_spectrum=nrange_spectrum,
                                    dB_cutoff_spectrum=dB_cutoff_spectrum)

    
if __name__ == "__main__":
    simulate_solitons()


    