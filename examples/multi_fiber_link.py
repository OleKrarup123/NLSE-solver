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


def simulate_multiple_fibers():

    """
    Use Chirped Pulse Amplification (CPA) to illustrate 
    how to set up a link consisting of multiple fibers.
    """
    ################ Set up time axis of simulation ################ 
    N = 2 ** 15  # Number of points on the time axis
    dt = 100e-15  # Time resolution [s]
    center_wavelength_m = 1550e-9 

    center_freq = wavelength_to_freq(center_wavelength_m)
    time_freq = TimeFreq(number_of_points=N,
                                time_step_s=dt,
                                center_frequency_Hz=center_freq)



    ################ Set up fibers ################ 
    
    alpha_dB_per_m = 0 #dB/m
    beta_list = [BETA2_AT_1550_NM_TYPICAL_SMF_S2_PER_M]  #[s^2/m,s^3/m,...]  s^(entry+2)/m
    gamma_W_per_m =  0 # 1/W/m
    length_m = 100e3 #100km
    number_of_steps = 2**6


    #In CPA, we first disperse the pulse anormalously, then amplify it and finally use
    #normal dispersion to re-compress it
    fiber_anorm_disp = FiberSpan(
        length_m,
        number_of_steps,
        gamma_W_per_m,
        beta_list,
        alpha_dB_per_m)

    fiber_gain = FiberSpan(
        length_m,
        number_of_steps,
        gamma_W_per_m,
        [0],
        alpha_dB_per_m=0.1e-3)

    fiber_norm_disp = FiberSpan(
        length_m,
        number_of_steps,
        gamma_W_per_m,
        [-beta_list[0]],
        alpha_dB_per_m)

    fiber_list_CPA = [fiber_anorm_disp,fiber_gain,fiber_norm_disp]
    fiber_link_CPA = FiberLink(fiber_list_CPA)

    #Set up input signal
    A_sqrt_W = np.sqrt(1e-3)
    duration_s=60e-12

    input_signal= InputSignal(time_freq,
                   amplitude_sqrt_W=A_sqrt_W,
                   duration_s=duration_s,
                   pulse_type="sinc",
                   describe_input_signal_flag=False,
                   FFT_tol=1e-5)


    ssfm_result_list_CPA = SSFM(fiber_link=fiber_link_CPA,input_signal=input_signal,show_progress_flag=True)

    nrange_pulse = 1500
    dB_cutoff_pulse=-40
    plot_first_and_last_pulse(ssfm_result_list_CPA,nrange=nrange_pulse,dB_cutoff=dB_cutoff_pulse)
    plot_pulse_matrix_2D(ssfm_result_list_CPA,nrange=nrange_pulse,dB_cutoff=dB_cutoff_pulse)
    make_chirp_gif(ssfm_result_list_CPA,nrange=nrange_pulse,chirp_range_Hz=[-10e9,10e9])


if __name__ == "__main__":
      simulate_multiple_fibers()