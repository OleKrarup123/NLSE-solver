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







def custom_signal():
    """
    Using the pulse_type argument, the InputSignal class can select from a set of built-in pulse shapes,
    such as Gaussian, sech, square, sinc etc. Sometimes however, you may need to define a custom
    input signal manually. 
    """
    os.chdir(os.path.realpath(os.path.dirname(__file__))) #Change directory to current folder to make sure output is saved in the correct place. 

    ################ Set up time axis of simulation ################ 
    N = 2 ** 15  # Number of points on the time axis
    dt = 10e-12  # Time resolution [s]
    center_wavelength_m = 1550e-9

    center_freq = wavelength_to_freq(center_wavelength_m)
    time_freq = TimeFreq(number_of_points=N,
                                time_step_s=dt,
                                center_frequency_Hz=center_freq)



    #Set up the initial input signal as a Gaussian pulse.
    input_signal= InputSignal(time_freq=time_freq,
                              duration_s=10e-9,
                              amplitude_sqrt_W=np.sqrt(0.05),
                              pulse_type="Gauss",
                              FFT_tol=1e-4)

    #Add a weaker continuous wave signal 3GHz above the carrier of the Gaussian pulse.  
    input_signal.pulse_field += get_pulse(time_freq.t_s(), 
                                          1.0,
                                          0.0,
                                          np.sqrt(1e-3),
                                          "CW",
                                          3e9)

    input_signal.update_spectrum() #IMPORTANT: After editing the field in the time domain, it is recommended to run update_spectrum() to ensure
                                   #           that the field in the frequency domain obtain the changes. 
    
    
    input_signal.describe_input_signal() #Plot input signal info after adding the CW component.

    ################ Set up fiber ################ 
    alpha_dB_per_m = 0 #dB/m
    beta_list = [0]   #[s^2/m,s^3/m,...]  s^(entry+2)/m
    
    gamma_W_per_m =  5e-3 # 1/W/m

    length_m = 20e3 #m
    number_of_steps = 2**7

    fiber = FiberSpan(
        length_m,
        number_of_steps,
        gamma_W_per_m,
        beta_list,
        alpha_dB_per_m)




    ################ Run split step simulations ################ 
    ssfm_result_list = SSFM(fiber=fiber,
                               input_signal=input_signal,
                               show_progress_flag=True)




    ################ Make plots of results ################ 
    nrange_pulse = 1000
    dB_cutoff_pulse = -40
    nrange_spectrum = 5200
    dB_cutoff_spectrum = -40


    #Without nonlinearity, the initial sech pulse will simply broaden in the time domain
    plot_everything_about_result(ssfm_result_list,
                                    nrange_pulse=nrange_pulse,
                                    dB_cutoff_pulse=dB_cutoff_pulse,
                                    nrange_spectrum=nrange_spectrum,
                                    dB_cutoff_spectrum=dB_cutoff_spectrum)

    #Conclusion: Because the Gaussian pulse and the CW signal overlap in the time domain
    #and have different frequecnies, we observe Four-Wave-Mixing in the spectral domain
    #as additional, evenly spaced spikes.

def custom_signal_advanced():
    """
    Now, we will see how to build a typical WDM telecom signal by adding pulses at fixed time spacings
    and at different carrier frequencies. 

    """
    os.chdir(os.path.realpath(os.path.dirname(__file__))) #Change directory to current folder to make sure output is saved in the correct place. 

    np.random.seed(123) #Fixed random seed for repeatable 64QAM modulation

    ################ Set up time axis of simulation ################ 
    N = 2 ** 16  # Number of points on the time axis
    dt = 0.1e-12  # Time resolution [s]
    center_wavelength_m = 1550e-9

    center_freq = wavelength_to_freq(center_wavelength_m)
    time_freq = TimeFreq(number_of_points=N,
                                time_step_s=dt,
                                center_frequency_Hz=center_freq)



    #Set up "empty" initial input signal by using a CW signal with zero power.
    input_signal= InputSignal(time_freq=time_freq,
                              duration_s=1,
                              amplitude_sqrt_W=np.sqrt(0),
                              pulse_type="CW",
                              FFT_tol=1e-3,
                              describe_input_signal_flag=False)


    #Set up WDM signal by specifying number of pulses for each frequency and the number of WDM frequencies
    #as well as time and frequency spacings. We use raised cosine pulses as these have a lot of 
    #nice properties (less inter-symbol interference, square-like spectra etc.). 
    #You can check that the implementation below makes sense by setting 
    #     N_pulses_in_time_domain = 1
    #     N_WDM_freqs = 10
    # or
    # 
    #     N_pulses_in_time_domain = 10
    #     N_WDM_freqs = 1
    #    
    N_pulses_in_time_domain = 10
    N_WDM_freqs = 4
    freq_spacing_Hz = 20e9
    time_spacing_s = 1/freq_spacing_Hz
    pulse_duration_s = time_spacing_s



    fig,ax=plt.subplots()

    #For each WDM frequency, add the specified number of evenly spaced, non-overlapping raised-cosine pulses 
    for freq_idx in range(N_WDM_freqs):
        for time_idx in range(N_pulses_in_time_domain):
            pulse = get_pulse(time_freq.t_s(), 
                              pulse_duration_s,
                              time_spacing_s*(time_idx-N_pulses_in_time_domain//2), #Evenly spaced time shifts roughly centered at zero
                              np.sqrt(1e-3),
                              "raised_cosine",
                              roll_off_factor=0.3,
                              freq_offset_Hz=freq_spacing_Hz*(freq_idx-N_WDM_freqs//2)) #Evenly spaced time shifts roughly centered at carrier

            pulse*= np.random.choice(QAM64) #Apply 64QAM modulation my multiplying by randomly sampled entries in the built-in array. 

            ax.plot(time_freq.t_s()/1e-9,get_power(pulse)) #Plot individual pulse.

            input_signal.pulse_field += pulse #Add pulse to total signal


    
    
    ax.plot(time_freq.t_s()/1e-9,get_power(input_signal.pulse_field),"r-",alpha=0.5,label="Power of total signal")
    ax.set_xlim(-time_spacing_s/1e-9*N_pulses_in_time_domain,time_spacing_s/1e-9*N_pulses_in_time_domain)
    ax.set_xlabel("Time [ns]")
    ax.set_ylabel("Power [W]")  
    ax.legend()
    plt.show()

    input_signal.update_spectrum() #IMPORTANT: After editing the field in the time domain, it is recommended to run update_spectrum() to ensure
                                   #           that the field in the frequency domain obtain the changes. 

    #Plot spectrum of input signal.
    fig,ax=plt.subplots()
    ax.plot(time_freq.f_rel_Hz()/1e12,get_power(input_signal.spectrum_field),"r-",alpha=0.5,label="Power of spectrum of total signal")
    ax.set_xlim(-freq_spacing_Hz/1e12*N_WDM_freqs,freq_spacing_Hz/1e12*N_WDM_freqs)
    ax.set_xlabel("Freq. [THz]")
    ax.set_ylabel("PSD [J/Hz]")
    ax.legend()
    plt.show()

    #Note that when many pulses with many modulation symbols and many frequencies are present, the signal in both the time 
    #and frequency domains begins to resemble white noise. This is similar to a room full of people individually having
    #coherent conversations but collectively sounding like noise. It is also the fundamental insight that allows
    #the so-called "Gaussian Noise Model" to estimate the nonlinear field generated in telecom fibers; instead of 
    #doing exact calculations with the actual field, do them with one resembling white noise and the result will be 
    #approximately the same!

    ################ Set up fiber ################ 

    #Set up a fiber with parameters typical for telecom applications. 
    alpha_dB_per_m = -0.22/1e3 #dB/m
    beta_list = [BETA2_AT_1550_NM_TYPICAL_SMF_S2_PER_M]   #[s^2/m,s^3/m,...]  s^(entry+2)/m
    
    gamma_W_per_m =  1e-3 # 1/W/m

    length_m = 100e3 #m
    number_of_steps = 2**10

    fiber = FiberSpan(
        length_m,
        number_of_steps,
        gamma_W_per_m,
        beta_list,
        alpha_dB_per_m,
        output_amp_dB=-length_m*alpha_dB_per_m, #Add ideal amplifier, which exactly compensates for loss at the end of each span
        output_disp_comp_s2=-beta_list[0]*length_m/2) #Add idealized output dispersion compensation (-b2/2(delta_w)^2) at the end of each span


    N_fibers = 3

    ################ Run split step simulations ################ 
    ssfm_result_list = SSFM(fiber=N_fibers*[fiber], #Concatenates the specified number of fibers together. 
                               input_signal=input_signal,
                               show_progress_flag=True)




    ################ Make plots of results ################ 
    nrange_pulse = 5000
    dB_cutoff_pulse = -40
    nrange_spectrum = 2000
    dB_cutoff_spectrum = -40


    plot_first_and_last_pulse(ssfm_result_list=ssfm_result_list,nrange=nrange_pulse,dB_cutoff=dB_cutoff_pulse)
    plot_pulse_matrix_2D(ssfm_result_list=ssfm_result_list,nrange=nrange_pulse,dB_cutoff=dB_cutoff_pulse)
    plot_first_and_last_spectrum(ssfm_result_list=ssfm_result_list,nrange=nrange_spectrum,dB_cutoff=dB_cutoff_spectrum)


    #Conclusion: When propagated through the fibers, the WDM signal retains its overall shape thanks to output amplifiers and 
    #dispersion compensation. The differences are due to the nonlinearity of the fiber.

    ### BONUS: Let's determine the nonlinear contribution to the field above by running the simulation
    ### again with a fiber where gamma=0 and subtracting the results!
    fiber_no_NL = FiberSpan(
        length_m,
        number_of_steps,
        0,
        beta_list,
        alpha_dB_per_m,
        output_amp_dB=-length_m*alpha_dB_per_m, #Add ideal amplifier, which exactly compensates for loss at the end of each span
        output_disp_comp_s2=-beta_list[0]*length_m/2) #Add idealized output dispersion compensation (-b2/2(delta_w)^2) at the end of each span

    ssfm_result_list_no_NL = SSFM(fiber=N_fibers*[fiber_no_NL], #Concatenates the specified number of fibers together. 
                               input_signal=input_signal,
                               show_progress_flag=True)



    plot_first_and_last_pulse(ssfm_result_list=ssfm_result_list_no_NL,nrange=nrange_pulse,dB_cutoff=dB_cutoff_pulse)
    plot_pulse_matrix_2D(ssfm_result_list=ssfm_result_list_no_NL,nrange=nrange_pulse,dB_cutoff=dB_cutoff_pulse)
    plot_first_and_last_spectrum(ssfm_result_list=ssfm_result_list_no_NL,nrange=nrange_spectrum,dB_cutoff=dB_cutoff_spectrum)


    final_field = ssfm_result_list[-1].pulse_matrix[-1,:]
    final_field_no_NL = ssfm_result_list_no_NL[-1].pulse_matrix[-1,:]
        
    nl_field = final_field-final_field_no_NL

    fig,ax=plt.subplots()
    ax.plot(time_freq.t_s()/1e-9,get_power(final_field),alpha=0.5,label="Final field with NL")
    ax.plot(time_freq.t_s()/1e-9,get_power(final_field_no_NL),alpha=0.5,label="Final field without NL")
    ax.plot(time_freq.t_s()/1e-9,get_power(nl_field),alpha=0.5,label="Nonlinear field")
    ax.set_xlim(-time_spacing_s/1e-9*N_pulses_in_time_domain,time_spacing_s/1e-9*N_pulses_in_time_domain)
    ax.set_xlabel("Time [ns]")
    ax.set_ylabel("Power [W]")  
    ax.legend()
    plt.show()

    #Conclusion: The nonlinear field is fairly small for the default pulse parameters.

if __name__ == "__main__":
    custom_signal()
    custom_signal_advanced()

    