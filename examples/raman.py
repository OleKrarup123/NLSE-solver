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


def plot_result_2D(ssfm_result_list,
                                    nrange_pulse=1000,
                                    dB_cutoff_pulse=-40,
                                    nrange_spectrum=8000,
                                    dB_cutoff_spectrum=-40):
    """Quick helper function to plot results for this example

    Args:
        ssfm_result_list (_type_): _description_
        nrange_pulse (int, optional): _description_. Defaults to 1000.
        dB_cutoff_pulse (int, optional): _description_. Defaults to -40.
        nrange_spectrum (int, optional): _description_. Defaults to 8000.
        dB_cutoff_spectrum (int, optional): _description_. Defaults to -40.
    """
    plot_pulse_matrix_2D(ssfm_result_list,nrange=nrange_pulse,dB_cutoff=dB_cutoff_pulse)
    plot_spectrum_matrix_2D(ssfm_result_list,nrange=nrange_spectrum,dB_cutoff=dB_cutoff_spectrum)


def simulate_raman(number_of_steps=2**8,show_plots=False):
    """
    Simulates the Raman effect.

    """
    ################ Set up time axis of simulation ################ 
    N = 2 ** 15  # Number of points on the time axis
    dt = 2.5e-15  # Time resolution [s]
    center_wavelength_m = 1300e-9 

    center_freq = wavelength_to_freq(center_wavelength_m)
    time_freq = TimeFreq(number_of_points=N,
                                time_step_s=dt,
                                center_frequency_Hz=center_freq)



    ################ Set up fibers ################ 
    alpha_dB_per_km = 0 #dB/m
    beta_list = [BETA2_AT_1550_NM_TYPICAL_SMF_S2_PER_M/50]  #[s^2/m,s^3/m,...]  s^(entry+2)/m
    
    gamma_W_per_m =  1e-3 # 1/W/m

    length_m = 10 #m
    number_of_steps = number_of_steps #=2**8 by default in this example

    #Define fibers with different Raman responses, which are otherwise identical. 
    fiber_no_raman = FiberSpan(
        length_m,
        number_of_steps,
        gamma_W_per_m,
        beta_list,
        alpha_dB_per_km,
        raman_model="none")

    fiber_agrawal = FiberSpan(
        length_m,
        number_of_steps,
        gamma_W_per_m,
        beta_list,
        alpha_dB_per_km,
        raman_model="agrawal")

    fiber_silica = FiberSpan(
        length_m,
        number_of_steps,
        gamma_W_per_m,
        beta_list,
        alpha_dB_per_km,
        raman_model="silica_exact")

    #Define a fiber, which is identical to the one above but approximates the Raman effect using the time derivative of the local power.
    #This shortens the run time, which is useful for iterating on your simulation parameters, but is inaccurate for pulse durations
    #comparable to the duration of the Raman response function (i.e. pulse durations on the scale of 10s of fs.)
    fiber_silica_approx = FiberSpan(
        length_m,
        number_of_steps,
        gamma_W_per_m,
        beta_list,
        alpha_dB_per_km,
        raman_model="silica_exact",
        approximate_raman_flag=True)

    #Custom Raman response implemented below. Can also be added directly inside the FiberSpan class. 
    fiber_custom = FiberSpan(
        length_m,
        number_of_steps,
        gamma_W_per_m,
        beta_list,
        alpha_dB_per_km,
        raman_model="custom")

    fiber_custom.raman_in_time_domain_func = lambda t_delay: ( (np.sin(t_delay/10e-15)+0.2*np.sin(t_delay/3e-15))*np.exp(-(t_delay/50e-15)**2))/1.1755126497259335e-14
    fiber_custom.relative_raman_contribution = 0.15

    ## When defining a custom raman response function, it must be normalized to have an area of 1. 
    ## The following 4 lines were used to determine the normalization constant of 1.1755126497259335e-14 that appears in the expression above.  
    #t_delay = np.linspace(0,1000e-15,2000) #Molecular vibrations rarely last longer than 1000 fs.
    #response_custom = fiber_custom.raman_in_time_domain_func(t_delay)
    #area=np.trapezoid(response_custom,t_delay)
    #print(area) #= 1.1755126497259335e-14


    #Built-in functions for plotting the raman response in the time domain
    #fiber_none.plot_raman_response()
    #fiber_agrawal.plot_raman_response()
    #fiber_silica.plot_raman_response()
    #fiber_custom.plot_raman_response()

    #Do the same thing manually to get both curves in the same chart
    t_delay = np.linspace(0,1000e-15,2000) #Molecular vibrations rarely last longer than 1000 fs.
    response_agrawal = fiber_agrawal.raman_in_time_domain_func(t_delay)
    response_silica = fiber_silica.raman_in_time_domain_func(t_delay)
    response_custom = fiber_custom.raman_in_time_domain_func(t_delay)
    if show_plots:
        fig,ax=plt.subplots()
        ax.plot(t_delay/1e-15,response_agrawal/1e15,label=fiber_agrawal.raman_model)
        ax.plot(t_delay/1e-15,response_silica/1e15,label=fiber_silica.raman_model)
        ax.plot(t_delay/1e-15,response_custom/1e15,label=fiber_custom.raman_model)
        ax.set_xlabel("T_delay [fs]")
        ax.set_ylabel("Response [1/fs]")
        ax.legend()
        plt.show()

    #The SSFM solver needs a "FiberLink" class, which is generated from a list of 1 or more fibers. 
    fiber_link_no_raman = FiberLink([fiber_no_raman]) 
    fiber_link_agrawal = FiberLink([fiber_agrawal])
    fiber_link_silica = FiberLink([fiber_silica])  
    fiber_link_silica_approx = FiberLink([fiber_silica_approx])  
    fiber_link_custom = FiberLink([fiber_custom]) 

    ################ Set up input signals ################ 
    A_sqrt_W = np.sqrt(5e3) #5kW peak power
    duration_s = 15e-15 #15fs duration

    input_signal= InputSignal(time_freq,
                   amplitude_sqrt_W=A_sqrt_W,
                   duration_s=duration_s,
                   pulse_type="gauss",
                   describe_input_signal_flag=False)

    ## Run all the simulations and plot the pulse and spectrum evolutions. You may comment out any of these lines to make
    ## it easier to compare two cases individually. 
    ssfm_result_list_no_raman       = SSFM(fiber_link=fiber_link_no_raman,input_signal=input_signal,show_progress_flag=show_plots)
    ssfm_result_list_agrawal        = SSFM(fiber_link=fiber_link_agrawal,input_signal=input_signal,show_progress_flag=show_plots)
    ssfm_result_list_silica         = SSFM(fiber_link=fiber_link_silica,input_signal=input_signal,show_progress_flag=show_plots)
    ssfm_result_list_silica_approx  = SSFM(fiber_link=fiber_link_silica_approx,input_signal=input_signal,show_progress_flag=show_plots)
    ssfm_result_list_custom         = SSFM(fiber_link=fiber_link_custom,input_signal=input_signal,show_progress_flag=show_plots)

    if show_plots:
        plot_result_2D(ssfm_result_list_no_raman)
        plot_result_2D(ssfm_result_list_agrawal)
        plot_result_2D(ssfm_result_list_silica)
        plot_result_2D(ssfm_result_list_silica_approx)
        plot_result_2D(ssfm_result_list_custom)
        

    ##Observations about result:
    #The Raman effect causes a continuous red-shift in the pulse, while anormalous dispersion (beta2<0) ensures that red light
    #is slower than blue light. The resulting Raman soliton drifts towards later times.  

    #Return list of ssfm_result_lists, so we can use them in the make_gifs.py example to demonstrate how to make animated gifs of signal evolution
    return [ssfm_result_list_no_raman,ssfm_result_list_agrawal,ssfm_result_list_silica,ssfm_result_list_custom]

if __name__ == "__main__":
    simulate_raman(show_plots=True)
