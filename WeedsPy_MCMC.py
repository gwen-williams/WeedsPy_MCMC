# =============================================================================
#                               WeedsPy_MCMC.py                               #
#                     Author: Dr. Gwenllian M. Williams                       #       
# =============================================================================
# Uses emcee to explore parameters of synthetic spectra created by WeedsPy.py #
# =============================================================================


import os
import numpy as np
import emcee
#import pickle
#import corner
import matplotlib
matplotlib.use('TkAgg')
from astropy.io import fits
from colorama import Fore, init
init(autoreset=True)


# =============================================================================
# Generate LTE synthetic models using Weeds, sampling the posterior distribution using emcee
# =============================================================================

DIR = '/Users/phygmw/Documents/Projects/ALMA/G19.01/Cycle2_paper2/emcee_multi/'


main_run = True


def read_params_file(paramfile):
    """
    Function to read in the parameter file:
    """

    params={}
    fileObject = open(paramfile)
    
    # Read parameter file into a dictionary, skipping empty lines and comments
    for line in fileObject:
        if(line=="\n"):
            continue
        if(line[0]=="#"):
           continue
        words=line.split()
        params[words[0].strip()] = words[2].strip()
        # Some catalog entries for molecules have whitespace:
        # In this case, need to merge words[2] and words[3], and include quotation marks so that CLASS can read it.
        if words[0]=='molecule':
            if ',' in words[2]:
                print('There is a comma in the molecule string')
                molname = '\"'+words[2]+' '+words[3]+'\"'
                params[words[0].strip()] = molname
    
    return params




def model(theta,ii=None,jj=None):
    """
    Function to evaluate the model of the data
    # For my case, won't need the age=age bit, as that is x, and my model doesn't depend on x
    """
    a1_N, a2_T, a3_v, a4_dv = theta 

    
    # Model properties
    coldens = a1_N*1e18 
    temp = a2_T 
    v_off = a3_v
    width = a4_dv
    
    # =========================================================================
    # Create the model file that CLASS will read:
    # =========================================================================
    
    # Create arrays to write to the .mdl files, which CLASS will read
    hd1 = np.array(["! species", "Ntot", "Tex", "source_size", "v_off", "width"])
    hd2 = np.array(["!", "(cm-2)","(K)", "('')", "(km/s)", "(km/s)"])
    hd3 = np.array([molecule,coldens,temp,source_size,v_off,width])
    d=[hd1,hd2,hd3]

    # Save .mdl file for CLASS to read.
    np.savetxt(DIR+"mdlfiles/temp_mdlfile.mdl",d,delimiter='\t',fmt='%s')
    
    # =========================================================================
    # Perform the model fit:
    # =========================================================================
    
    # Call the CLASS script to perform the modelling
    Sic.comm('@WeedsPy_MCMC.class')
    
    
    # =========================================================================
    # Read in the model data, and return it:
    # =========================================================================
    
    mod = fits.open(DIR+'synthfiles/temp_synthfile.fits')[0]
    modeldata = np.asarray(mod.data)[0,0,0,0:]
    modeldata[np.where(np.isnan(modeldata)==True)] = 0. # Remove nans from the model
    
    if main_run == True:
        """
        Delete the intermediate files whilst trying to converge on a solution i.e. when main_run == True.
        If wanting to plot the best model after completion of the script, i.e. main_run == False, then don't delete these files
        """
        os.remove(DIR+'mdlfiles/temp_mdlfile.mdl')
        os.remove(DIR+'synthfiles/temp_synthfile.fits')#+modname)
        os.remove(DIR+'synthfiles/temp_synthfile.30m')#+modname[:-5]+'.30m')
    
    if main_run == False:
        """
        Keep the best fitting model written to file, just rename the file to that of the pixel ID
        """
        os.rename(DIR+'synthfiles/temp_synthfile.fits',DIR+'synthfiles/best_synthfile_i{0}_j{1}.fits'.format(ii,jj))
        os.rename(DIR+'synthfiles/temp_synthfile.30m',DIR+'synthfiles/best_synthfile_i{0}_j{1}.30m'.format(ii,jj))
        os.rename(DIR+'mdlfiles/temp_mdlfile.mdl',DIR+'mdlfiles/best_mdlfile_i{0}_j{1}.mdl'.format(ii,jj))
    
    
    return modeldata



def lnlike(theta, x, y, yerr):
    """
    Function calculating the "likelihood" of the model
    """
    return -0.5 * np.sum(((y - model(theta))/yerr) ** 2)



def lnprior(theta):
    """
    Function to set the priors, and check that all the values walked to are within the priors. Flat priors.
    """
    a1_N, a2_T, a3_v, a4_dv = theta #, a5_s = theta
    if 0.05 < a1_N < 3.0 and 80 < a2_T < 230 and 56.0 < a3_v < 64.0 and 2.0 < a4_dv < 8.0:
        return 0.0
    return -np.inf




def lnprob(theta, x, y, yerr):
    """
    Function to combine lnlike() and lnprior(), see if theta is within the priors, and calculate the likelihood of the model.
    """
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)



def main_mcmc(p0,nwalkers,nburn,nprod,ndim,lnprob,data):
    """
    Function to run the actual MCMC
    """
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=data)

    print(Fore.BLUE + "Running burn-in..." + Fore.WHITE)
    pos0, prob0, state0 = sampler.run_mcmc(p0, nburn)
    sampler.reset()

    print(Fore.BLUE + "Running production..." + Fore.WHITE) # Start the production run from the last of the burn in runs, and with a clear sampler
    pos, prob, state = sampler.run_mcmc(pos0, nprod)

    return sampler, pos, prob, state




if __name__ == '__main__':
        
 
    # =============================================================================
    # READ PARAMETER FILE
    # =============================================================================
    paramfile=DIR+'paramfile.txt'
    
    params = read_params_file(paramfile)
        
    # Read input parameters from dictionary into arrays, and set their data types
    catalog_name = params['catalog_name']
    molecule = params['molecule']
    source_size = float(params['source_size'])
    freq_lower = int(params['freq_lower'])
    freq_upper = int(params['freq_upper'])
    mean_freq = float(params['mean_freq'])
    peak_spec = int(params['peak_spec'])
    telescope_diameter = int(params['telescope_diameter'])
    bmaj = float(params['bmaj'])
    bmin = float(params['bmin'])
    rms_noise = float(params['rms_noise'])
    rms_noise_unit = str(params['rms_noise_unit'])
    nburn = int(params['nburns'])                
    nprod = int(params['nprod'])             
    nwalkers = int(params['nwalk'])              
    i0_col = float(params['i0_col'])               
    col_base = params['col_base']         
    i0_tex = float(params['i0_tex'])           
    i0_vsys = float(params['i0_vsys'])            
    i0_dV = float(params['i0_dV'])         


    
    
    # =============================================================================
    # Error:
    # =============================================================================
    if rms_noise_unit == 'Jy/beam':
        mean_freq_GHz = mean_freq/1e3
        rms_kelvin = (1.222e3*rms_noise*1e3)/((mean_freq_GHz**2)*bmaj*bmin) # Error, i.e. the rms noise
    elif rms_noise_unit == 'K':
        print('RMS noise already in temperature units')

    
    
    # Define CLASS variables:
    Sic.comm('define character*128 peakname')
    Sic.comm('define character*128 molecule')
    Sic.comm('define real freq_lower')
    Sic.comm('define real freq_upper')
    Sic.comm('define real telescope_diameter')
    Sic.comm('define character*128 t_cont')
    
    # Pass some of these input parameters to CLASS, and set some plot ranges:
    Sic.comm('let molecule '+str(molecule))
    Sic.comm('use in '+str(catalog_name))
    Sic.comm('set unit f')
    Sic.comm('let freq_lower '+str(freq_lower))
    Sic.comm('let freq_upper '+str(freq_upper))
    Sic.comm('set mod y -1 '+str(peak_spec))
    Sic.comm('let telescope_diameter '+str(telescope_diameter))

    # Pixel coordinates, continuum brightness temperature:
    pixi = np.loadtxt(DIR+'i_coords_200sigma.txt',delimiter=',',dtype='str')
    pixj = np.loadtxt(DIR+'j_coords_200sigma.txt',delimiter=',',dtype='str')
    tcont = np.loadtxt(DIR+'tcont_200sigma.txt',delimiter=',',dtype='str')
    npix = len(pixi)
    

    for ind in range(0,npix):
        
        # =============================================================================
        # Read in the spectrum of the ith, jth pixel: 
        # =============================================================================
        peakname = 'spectra/spec_i{0}_j{1}.30m'.format(pixi[ind],pixj[ind])
    
        obs = fits.open(DIR+peakname[:-4]+'.fits')[0]
        obsdata = np.asarray(obs.data).flatten()
        nchans = np.shape(obsdata)[0]
        
        # Set the continuum brightness, and the .30m file name:
        Sic.comm('let peakname '+peakname)
        Sic.comm('let t_cont '+tcont[ind])
        
        
        # =============================================================================
        # Set up of all the emcee specifics, 
        # =============================================================================
        x = np.arange(0,len(obsdata))
        data = (x,obsdata,rms_kelvin) # Vector theta of the data (x,y,yerr)
        initial = np.array([i0_col,i0_tex,i0_vsys,i0_dV]) # Initial guesses for each parameter of the model
        ndim = len(initial) # Number of model parameters
        p0 = [np.array(initial) + 1e-3 * np.random.randn(ndim) for i in range(nwalkers)] # Set-up beginning walk locations as a gaussian around the parameters.

        
        # Call the main procedure:
        sampler, pos, prob, state = main_mcmc(p0,nwalkers,nburn,nprod,ndim,lnprob,data)
        
        
        
        # =============================================================================
        # Save the sampler output directly to a csv, with already flattened output, rather than pickling the sampler class object
        # =============================================================================
        samples = sampler.flatchain
        np.savetxt(DIR+'samples/samples_flatchain_i{0}_j{1}.csv'.format(pixi[ind],pixj[ind]),samples,delimiter=',')
        
                
        # Save the maximum likelihood model parameters
        theta_max = samples[np.argmax(sampler.flatlnprobability)]
        np.savetxt(DIR+'theta_max/theta_max_flatlnprob_i{0}_j{1}.csv'.format(pixi[ind],pixj[ind]),theta_max,delimiter=',')
        
        # Generate model files of the maximum likelihood model:
        main_run = False
        max_lnlike_model = model(theta_max,ii=pixi[ind],jj=pixj[ind])
        #max_lnlike_model = model(theta_max,ii=pixi[i],jj=pixj[i])
        main_run = True
        
    






