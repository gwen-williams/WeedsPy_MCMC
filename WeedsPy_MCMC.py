# =============================================================================
#                               WeedsPy_MCMC.py                               #
#                     Author: Dr. Gwenllian M. Williams                       #       
# =============================================================================
# Uses emcee to explore parameters of synthetic spectra created by WeedsPy.py #
# =============================================================================
import os
import numpy as np
import emcee
from astropy.io import fits
import matplotlib.pyplot as plt
import corner
#matplotlib.use('TkAgg')


""" 
Below is the WeedsPy_MCMC class object.

Both the class, and the executable script, must be present in the same .py file.
This is because of the nesting of the Python/Gildas Sic commands. 
Sic commands implemented in a separate .py script will not execute.

"""



def read_params_file(paramfile):
    """
    Function to read in the parameter file.
    This function is defined outside of the WeedsPy_MCMC class object, 
    because it passes these are the input arguments that will be passed 
    to the WeedsPy_MCMC class object
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




class WeedsPy_MCMC:
    
    
    def __init__(self, 
        catalog_name, 
        molecule_name,
        source_size,
        rms_noise,
        rms_noise_unit,
        n_iter,
        n_burn,
        n_walkers,
        n_dim,
        column_base,
        priors,
        initial,
        mean_freq,
        bmaj,
        bmin,
        freq_unit,
        freq_low,
        freq_up,
        peak_spec,
        telescope_diameter,
        vel_sys,
        vel_width,
        
    ):
        
        
        """
        
        This class takes in the spectrum to the modelled by Weeds, and uses
        emcee to determine the parameters of the synthetic spectra.
        This is not a best fit spectrum, rather a highest likelihood spectrum.
        
        Parameters:
            
        """
    
        self.catalog_name = catalog_name
        self.molecule_name = molecule_name
        self.source_size = source_size
        self.rms_noise = rms_noise
        self.rms_noise_unit = rms_noise_unit
        self.n_iter = n_iter
        self.n_burn = n_burn
        self.n_walkers = n_walkers
        self.n_dim = n_dim
        self.column_base = column_base
        self.priors = priors
        self.initial = initial
        self.mean_freq = mean_freq
        self.bmaj = bmaj
        self.bmin = bmin
        self.freq_unit = freq_unit
        self.freq_low = freq_low
        self.freq_up = freq_up
        self.peak_spec = peak_spec
        self.telescope_diameter = telescope_diameter
        self.vel_sys = vel_sys
        self.vel_width = vel_width


    def set_gildas_variables(self):
        """
        Function to set the variables in Gildas/CLASS
        """
        # Define CLASS variables:
        Sic.comm('define character*128 peakname')
        Sic.comm('define character*128 molecule')
        Sic.comm('define real freq_lower')
        Sic.comm('define real freq_upper')
        Sic.comm('define real telescope_diameter')
        Sic.comm('define character*128 t_cont')
        
        # Pass some of these input parameters to CLASS, and set some plot ranges:
        Sic.comm('let molecule '+str(self.molecule_name))
        Sic.comm('use in '+str(self.catalog_name))
        Sic.comm('set unit f')
        Sic.comm('let freq_lower '+str(self.freq_low))
        Sic.comm('let freq_upper '+str(self.freq_up))
        Sic.comm('set mod y -1 '+str(self.peak_spec))
        Sic.comm('let telescope_diameter '+str(self.telescope_diameter))
        
        return
    
    
    def make_gildas_mdl_file(self,theta):#,ii=None,jj=None):
        """
        This is a function to make the Gildas/CLASS readable text file (.mdl extension) that contains the synthetic spectrum model parameters walked to by the emcee.
        """
        
        # If you want to model column density, temperature, centroid velocity and velocity width:
        if self.n_dim == 4:
        
            # Unpack the theta array of the walked to parameter values:
            a1_N, a2_T, a3_v, a4_dv = theta
            
            # Model properties
            coldens = a1_N*self.column_base
            temp = a2_T 
            v_off = a3_v
            v_width = a4_dv
            
            # =========================================================================
            # Create the model file that CLASS will read:
            # =========================================================================
            
            # Create arrays to write to the .mdl files, which CLASS will read
            hd1 = np.array(["! species", "Ntot", "Tex", "source_size", "v_off", "width"])
            hd2 = np.array(["!", "(cm-2)","(K)", "('')", "(km/s)", "(km/s)"])
            hd3 = np.array([self.molecule_name,coldens,temp,self.source_size,v_off,v_width])
            mdl=[hd1,hd2,hd3]
        
        # If you want to model all 5 parameters:
        elif self.n_dim == 5:
            
            a1_N, a2_T, a3_v, a4_dv, a5_ss = theta
            
            # Model properties
            coldens = a1_N*self.column_base
            temp = a2_T 
            v_off = a3_v
            v_width = a4_dv
            source_size = a5_ss
            
            # Create arrays to write to the .mdl files, which CLASS will read
            hd1 = np.array(["! species", "Ntot", "Tex", "source_size", "v_off", "width"])
            hd2 = np.array(["!", "(cm-2)","(K)", "('')", "(km/s)", "(km/s)"])
            hd3 = np.array([self.molecule_name,coldens,temp,source_size,v_off,v_width])
            mdl=[hd1,hd2,hd3]
        
        # If you want to model only column density and velocity width:
        elif self.n_dim == 2:
            
            a1_N, a2_T = theta
            
            # Model properties
            coldens = a1_N*self.column_base
            temp = a2_T 
            
            # Create arrays to write to the .mdl files, which CLASS will read
            hd1 = np.array(["! species", "Ntot", "Tex", "source_size", "v_off", "width"])
            hd2 = np.array(["!", "(cm-2)","(K)", "('')", "(km/s)", "(km/s)"])
            hd3 = np.array([self.molecule_name,coldens,temp,self.source_size,self.vel_sys,self.vel_width])
            mdl=[hd1,hd2,hd3]
        
        # If n_dim is not 4, 5 or 2, raise a ValueError
        else:
            raise ValueError("n_dim must be either 4, 5 or 2. Code does not yet support modelling of 3 or only 1 parameters.")

        return mdl
        
        

    def main_emcee(self,xdata,ydata):
        """
        This is the main function that will carry out the emcee process
        """
        
    
        def model(theta,ii=None,jj=None):#,save_best=None):#,best_model_flag=None):
            """
            Function to evaluate the model of the data
            """
            
            """
            # If you want to model column density, temperature, centroid velocity and velocity width:
            if self.n_dim == 4:
            
                # Unpack the theta array of the walked to parameter values:
                a1_N, a2_T, a3_v, a4_dv = theta
                
                # Model properties
                coldens = a1_N*self.column_base
                temp = a2_T 
                v_off = a3_v
                v_width = a4_dv
                
                # =========================================================================
                # Create the model file that CLASS will read:
                # =========================================================================
                
                # Create arrays to write to the .mdl files, which CLASS will read
                hd1 = np.array(["! species", "Ntot", "Tex", "source_size", "v_off", "width"])
                hd2 = np.array(["!", "(cm-2)","(K)", "('')", "(km/s)", "(km/s)"])
                hd3 = np.array([self.molecule_name,coldens,temp,self.source_size,v_off,v_width])
                mdl=[hd1,hd2,hd3]
            
            # If you want to model all 5 parameters:
            elif self.n_dim == 5:
                
                a1_N, a2_T, a3_v, a4_dv, a5_ss = theta
                
                # Model properties
                coldens = a1_N*self.column_base
                temp = a2_T 
                v_off = a3_v
                v_width = a4_dv
                source_size = a5_ss
                
                # Create arrays to write to the .mdl files, which CLASS will read
                hd1 = np.array(["! species", "Ntot", "Tex", "source_size", "v_off", "width"])
                hd2 = np.array(["!", "(cm-2)","(K)", "('')", "(km/s)", "(km/s)"])
                hd3 = np.array([self.molecule_name,coldens,temp,source_size,v_off,v_width])
                mdl=[hd1,hd2,hd3]
            
            # If you want to model only column density and velocity width:
            elif self.n_dim == 2:
                
                a1_N, a2_T = theta
                
                # Model properties
                coldens = a1_N*self.column_base
                temp = a2_T 
                
                # Create arrays to write to the .mdl files, which CLASS will read
                hd1 = np.array(["! species", "Ntot", "Tex", "source_size", "v_off", "width"])
                hd2 = np.array(["!", "(cm-2)","(K)", "('')", "(km/s)", "(km/s)"])
                hd3 = np.array([self.molecule_name,coldens,temp,self.source_size,self.vel_sys,self.vel_width])
                mdl=[hd1,hd2,hd3]
            
            else:
                raise ValueError("n_dim must be either 4, 5 or 2. Code does not yet support modelling of 3 or only 1 parameters.")
            """
            
            # Make the .mdl file that Gildas/CLASS will need to make the synthetic spectrum and save it
            mdl = self.make_gildas_mdl_file(theta)
            np.savetxt("mdlfiles/temp_mdlfile.mdl",mdl,delimiter='\t',fmt='%s')
            
            # =========================================================================
            # Call the Gildas/CLASS script to perform the modelling
            # =========================================================================
            Sic.comm('@WeedsPy_MCMC.class')
            
            # =========================================================================
            # Read in the model data, and return it:
            # =========================================================================
            mod = fits.open('synthfiles/temp_synthfile.fits')[0]
            modeldata = np.asarray(mod.data)[0,0,0,0:]
            modeldata[np.where(np.isnan(modeldata)==True)] = 0. # Remove nans from the model
            modeldata[860:910]=0.  # Remove the maser line from the evaluation of the likelihood
            
            #if main_run == True:
            """
            Delete the intermediate files whilst trying to converge on a solution i.e. when main_run == True.
            If wanting to plot the best model, i.e. main_run == False, then don't delete these files
            """
            #os.remove(DIR+"mdlfiles/CH3OH_{0}cm-2_{1}K.mdl".format(coldens,temp))
            os.remove('mdlfiles/temp_mdlfile.mdl')
            os.remove('synthfiles/temp_synthfile.fits')#+modname)
            os.remove('synthfiles/temp_synthfile.30m')#+modname[:-5]+'.30m')
        
            #if main_run == False:
            #if save_best:
            #    """
            #    Keep the best fitting model written to file, just rename the file:
            #    """
            #    print("Saving the synthetic spectrum with the highest likelihood parameter values")
            #    os.rename('synthfiles/temp_synthfile.fits','synthfiles/best_synthfile_i{0}_j{1}.fits'.format(ii,jj))
            #    os.rename('synthfiles/temp_synthfile.30m','synthfiles/best_synthfile_i{0}_j{1}.30m'.format(ii,jj))
            #    os.rename('mdlfiles/temp_mdlfile.mdl','mdlfiles/best_mdlfile_i{0}_j{1}.mdl'.format(ii,jj))
            
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
            # If you are modelling 4 parameters: column density, temperature, centroid velocity and velocity width

            if self.n_dim == 4:
                a1_N, a2_T, a3_v, a4_dv = theta #, a5_s = theta
                if self.priors[0] < a1_N < self.priors[1] and self.priors[2] < a2_T < self.priors[3] and self.priors[4] < a3_v < self.priors[5] and self.priors[6] < a4_dv < self.priors[7]:
                    return 0.0
            
            # If you're modelling all 5 parameters:
            elif self.n_dim == 5:
                a1_N, a2_T, a3_v, a4_dv, a5_ss = theta #, a5_s = theta
                if self.priors[0] < a1_N < self.priors[1] and self.priors[2] < a2_T < self.priors[3] and self.priors[4] < a3_v < self.priors[5] and self.priors[6] < a4_dv < self.priors[7] and self.priors[8] < a5_ss < self.priors[9]:
                    return 0.0
            
            # If you're modelling only 2 parameters: column density and temperature
            elif self.n_dim == 2:
                a1_N, a2_T = theta #, a5_s = theta
                if self.priors[0] < a1_N < self.priors[1] and self.priors[2] < a2_T < self.priors[3]:
                    return 0.0
            else:
                raise ValueError("n_dim must be either 4, 5 or 2. Code does not yet support modelling of 3 or only 1 parameters.")

            
            return -np.inf
        
        
        def lnprob(theta, x, y, yerr):
            """
            Function to combine lnlike() and lnprior(), see if theta is within the priors, and calculate the likelihood of the model.
            """
            lp = lnprior(theta)
            if not np.isfinite(lp):
                return -np.inf
            return lp + lnlike(theta, x, y, yerr)
    
        
        # Set the starting locations of the walkers as a gaussian around the initial guesses
        p0 = [np.array(self.initial) + 1e-3 * np.random.randn(self.n_dim) for i in range(self.n_walkers)] 

        # Setup the vector of theta (i.e. data)
        data = (xdata, ydata, self.rms_noise)

        sampler = emcee.EnsembleSampler(self.n_walkers, self.n_dim, lnprob, args=data)
    
        # Start burn-in run
        print("Running burn-in...")
        pos0, prob0, state0 = sampler.run_mcmc(p0, self.n_burn)
        sampler.reset()
        
        # Start the production run from the last of the burn in runs, and with a clear sampler
        print("Running production...") 
        pos, prob, state = sampler.run_mcmc(pos0, self.n_iter)
                
        
        return sampler

    
    def run_emcee(self,xdata,ydata,ii,jj):
        """
        Function that runs the main_emcee process, and processes the outcome of the emcee
        """
        
        def save_highest_likelihood_spectrum(theta,ii,jj):
            """
            This is a function to save the highest likelihood spectrum made by Gildas/CLASS, because we delete intermediate spectrum files while the emcee is walking to save space.
            """
            
            mdl = self.make_gildas_mdl_file(theta)
            np.savetxt('mdlfiles/temp_mdlfile.mdl'.format(ii,jj),mdl,delimiter='\t',fmt='%s')
            
            # Make the spectrum:
            Sic.comm('@WeedsPy_MCMC.class')
            
            # Rename the temporary spectrum files:
            os.rename('mdlfiles/temp_mdlfile.mdl','mdlfiles/best_mdlfile_i{0}_j{1}.mdl'.format(ii,jj))
            os.rename('synthfiles/temp_synthfile.fits','synthfiles/best_synthfile_i{0}_j{1}.fits'.format(ii,jj))
            os.rename('synthfiles/temp_synthfile.30m','synthfiles/best_synthfile_i{0}_j{1}.30m'.format(ii,jj))
            
            return
        
        # Run the emcee
        sampler = self.main_emcee(xdata,ydata)
        
        # Get the run chains
        samples = sampler.get_chain()
        
        # Get the flat samples
        samples_flat = sampler.get_chain(flat=True) # thin=10, discard=self.n_burn
        np.savetxt('samples/samples_flatchain_i{0}_j{1}.csv'.format(ii,jj),samples_flat,delimiter=',')
        

        # Find the autocorrelation time of the run:
        try:
            tau = sampler.get_autocorr_time()
        except:
            print(
                f'**Warning** Walker chain is shorter than 50 times the integrated autocorrelation time for {self.n_dim} parameter(s). Proceed with caution and consider running a longer chain if you can.'
            )
            tau = np.inf
        print(
            f'The autocorrelation time is {tau}. Consider running the chains for at least 10 times as long as this.'
        )
        print('----------------------------------------------------------')

        
        # Get the results
        theta_max = sampler.flatchain[np.argmax(sampler.flatlnprobability)]
        np.savetxt('theta_max/theta_max_flatlnprob_i{0}_j{1}.csv'.format(ii,jj),theta_max,delimiter=',')


        # Get the highest likelihood parameters and their errors
        results = []
        results_err_up = []
        results_err_low = []
        # Loop over each dimension
        for i in range(self.n_dim):
            mcmc_params = np.percentile(samples_flat[:,i],[16,50,84])
            errors = np.diff(mcmc_params)
            
            # Store the results to arrays
            results.append(mcmc_params[1])
            results_err_low.append(errors[0])
            results_err_up.append(errors[1])
        
        # Save the best spectrum
        save_highest_likelihood_spectrum(theta_max,ii,jj)
        
        """
        # Highest likelihood parameters:
        best_N = results[0]
        best_T = results[1]
        best_V = results[2]
        best_dV = results[3]
        # Errors on highest likelihood parameters:
        best_N_err_up = results_err_up[0]
        best_N_err_low = results_err_low[0]
        best_T_err_up = results_err_up[1]
        best_T_err_low = results_err_low[1]
        best_V_err_up = results_err_up[2]
        best_V_err_low = results_err_low[2]
        best_dV_err_up = results_err_up[3]
        best_dV_err_low = results_err_low[3]        
       
        return samples, samples_flat, best_N, best_N_err_up, best_N_err_low, best_T, best_T_err_up, best_T_err_low, best_V, best_V_err_up, best_V_err_low, best_dV, best_dV_err_up, best_dV_err_low
        """
        return samples, samples_flat, results, results_err_up, results_err_low
        
    
    def plot_corner(self,samples_flat,ii,jj):
        """
        Function to make a corner plot
        """
        # Change the plot labels based on the number of fitted parameters:

        if self.n_dim == 4:
            labels = [r"N(mol) ($\times$ %s cm$^{-2}$)" % (self.column_base),r"T$_{\mathrm{ex}}$ (K)",r"V$_{\mathrm{cen}}$ (km s$^{-1}$)","$\Delta$V (km s$^{-1}$)"]
        elif self.n_dim == 5:
            labels = [r"N(mol) ($\times$ %s cm$^{-2}$)" % (self.column_base),r"T$_{\mathrm{ex}}$ (K)",r"V$_{\mathrm{cen}}$ (km s$^{-1}$)","$\Delta$V (km s$^{-1}$)", r"$\theta_{\mathrm{source}}$ (arcsec)"]
        elif self.n_dim == 2:
            labels = [r"N(mol) ($\times$ %s cm$^{-2}$)" % (self.column_base),r"T$_{\mathrm{ex}}$ (K)"]
        else:
            raise ValueError("n_dim must be either 4, 5 or 2. Code does not yet support modelling of 3 or only 1 parameters.")

        fig = corner.corner(samples_flat, 
                            labels=labels, 
                            show_titles=True,
                            plot_datapoints=True,
                            quantiles=[0.16,0.5,0.84],
                            smooth = 1.0,
                            levels=(1-np.exp(-0.5),1-np.exp(-1.0),1-np.exp(-2.0)),
                            label_kwargs={'fontsize':11},
                            title_kwargs={'fontsize':11},
                            )
        fig.savefig('plots/corner_i{0}_j{1}.pdf'.format(ii,jj),bbox_inches='tight')
        
        return
    
    
    def plot_chains(self,samples,ii,jj):
        """
        Function to make a plot of the walked chains
        """
        fig, axs = plt.subplots(self.n_dim, figsize=(10,7), sharex=True)
        axsf = axs.flatten()
        
        # Change the plot labels based on the number of fitted parameters:
        if self.n_dim == 4:
            labels = [r"N(mol) ($\times$ %s cm$^{-2}$)" % (self.column_base),r"T$_{\mathrm{ex}}$ (K)",r"V$_{\mathrm{cen}}$ (km s$^{-1}$)","$\Delta$V (km s$^{-1}$)"]
        elif self.n_dim == 5:
            labels = [r"N(mol) ($\times$ %s cm$^{-2}$)" % (self.column_base),r"T$_{\mathrm{ex}}$ (K)",r"V$_{\mathrm{cen}}$ (km s$^{-1}$)","$\Delta$V (km s$^{-1}$)", r"$\theta_{\mathrm{source}}$ (arcsec)"]
        elif self.n_dim == 2:
            labels = [r"N(mol) ($\times$ %s cm$^{-2}$)" % (self.column_base),r"T$_{\mathrm{ex}}$ (K)"]
        else:
            raise ValueError("n_dim must be either 4, 5 or 2. Code does not yet support modelling of 3 or only 1 parameters.")


        # Iterate over each dimension
        for i in range(self.n_dim):
            # Plot the chain traces
            axsf[i].plot(samples[:,:,i],"k",alpha=0.4)
            # Plot vertical line denoting the end of the burn-in
            axsf[i].axvline(self.n_burn,np.min(samples[:,:,i]),np.max(samples[:,:,i]),linestyle='--',color='#bbbbbb')
            # General params
            axsf[i].set_xlim(0,len(samples))
            axsf[i].set_ylabel(labels[i])
            axsf[i].tick_params(top=True,right=True)
            axsf[i].yaxis.set_tick_params(which='minor', right=True)
            axsf[i].xaxis.set_tick_params(which='minor', top=True)
        # Only put x label on the final subplot    
        axsf[i].set_xlabel("Iterations")
        # Save
        fig.savefig('plots/chains_i{0}_j{1}.pdf'.format(ii,jj),bbox_inches='tight')
        
        return fig
    
    
    def convert_noise_to_K(self):
        """
        Convert rms noise from mJy/beam, Jy/beam or uJy/beam to Kelvin
        """
        # Check the units of the frequency, and convert to GHz
        if self.freq_unit == 'MHz':
            mean_freq_GHz = self.mean_freq/1e3
        elif self.freq_unit == 'Hz':
            mean_freq_GHz = self.mean_freq/1e9
        elif self.freq_unit == 'kHz':
            mean_freq_GHz = self.mean_freq/1e6
        else:
            # If not in Hz, kHz or MHz, assume already in GHz
            mean_freq_GHz = self.mean_freq
        
        # Check the unit of the rms noise, and do conversion
        if self.rms_noise_unit == 'Jy/beam':
            rms_kelvin = (1.222e3*self.rms_noise*1e3)/((mean_freq_GHz**2)*self.bmaj*self.bmin)
        elif self.rms_noise_unit == 'mJy/beam':
            rms_kelvin = (1.222e3*self.rms_noise)/((mean_freq_GHz**2)*self.bmaj*self.bmin)
        elif self.rms_noise_unit == 'uJy/beam':
            rms_kelvin = (1.222e3*self.rms_noise/1e3)/((mean_freq_GHz**2)*self.bmaj*self.bmin)
        else:
            # If not Jy/beam, mJy/beam or uJy/beam, assume already in Kelvin
            rms_kelvin = rms_noise
        
        self.rms_noise = rms_kelvin
        self.rms_noise_unit = 'K'
        
        return rms_kelvin
    
    
    
            
if __name__ == '__main__':           
            
    # Read the parameter file:
    params = read_params_file('paramfile.txt')

    # Get parameters from the parameter file and set their types:
    catalog_name = params['catalog_name']
    molecule_name = params['molecule']
    #source_size = float(params['source_size'])
    rms_noise = float(params['rms_noise'])
    rms_noise_unit = str(params['rms_noise_unit'])
    n_burn = int(params['nburn'])
    n_iter = int(params['nprod'])
    n_walkers = int(params['nwalkers'])
    #column_base = int(params['col_base'])
    mean_freq = float(params['mean_freq'])
    bmaj = float(params['bmaj'])
    bmin = float(params['bmin'])
    freq_unit = str(params['freq_unit'])
    freq_low = int(params['freq_lower'])
    freq_up = int(params['freq_upper'])
    peak_spec = int(params['peak_spec'])
    telescope_diameter = int(params['telescope_diameter'])
    
    # The power base of the column density. Easier to run emcee without the power and multiply it by this for the modelling
    column_base = 1e18
    
    # List of the priors, each pairs corresponds to upper and lower bounds for
    # column density, temperature, systemic velocity and velocity width
    priors = [0.05,3.0,80,230]#,56.0,64.0,2.0,8.0]
    
    # Initial locations for walkers:
    # column density, temperature, systemic velocity, velocity width
    initial = np.array([1.0,120])#,60.0,5.8])
    
    # Number of parameters to fit:
    n_dim = len(initial)
    
    # Set the value of the parameters that will NOT be modelled with emcee:
    # If you won't be modelling them, set them to None
    source_size = float(params['source_size'])
    vel_sys = 59.7
    vel_width = 6.0
    
    # Change for now:
    n_burn = 5
    n_iter = 5
    
    # Initialise the class
    tt = WeedsPy_MCMC(catalog_name = catalog_name, 
                      molecule_name = molecule_name, 
                      source_size = source_size, 
                      rms_noise = rms_noise, 
                      rms_noise_unit = rms_noise_unit, 
                      n_iter = n_iter, 
                      n_burn = n_burn, 
                      n_walkers = n_walkers,
                      n_dim = n_dim,
                      column_base = column_base, 
                      priors = priors, 
                      initial = initial,
                      mean_freq = mean_freq,
                      bmaj = bmaj,
                      bmin = bmin,
                      freq_unit = freq_unit,
                      freq_low = freq_low,
                      freq_up = freq_up,
                      peak_spec = peak_spec,
                      telescope_diameter = telescope_diameter,
                      vel_sys = vel_sys,
                      vel_width = vel_width)
    
    
    # Set the variables in GILDAS
    tt.set_gildas_variables()
    
    
    # Flag to make plots or not (True = yes, False = no)
    make_plots = True
    
    # Read in the data    
    pixi = np.array([288])
    pixj = np.array([296])
    tcont = np.array([22.062])
    npix = len(pixi)
    
    
    # Loop over each of the pixels
    for ind in range(0,npix):
        
        peakname = 'spec_i{0}_j{1}.30m'.format(pixi[ind],pixj[ind])
        obs = fits.open(peakname[:-4]+'.fits')[0]
        ydata = np.asarray(obs.data).flatten()
        xdata = np.arange(0,len(ydata))
        
        # Set the continuum brightness, and the .30m file name:
        Sic.comm('let peakname '+peakname)
        Sic.comm('let t_cont '+str(tcont[ind]))

        # Convert the noise to Kelvin:
        tt.convert_noise_to_K()
        
        # don't need the x and y data defined in the self, can define them here
        #samples, samples_flat, N, N_err_up, N_err_low, T, T_err_up, T_err_low, V, V_err_up, V_err_low, dV, dV_err_up, dV_err_low = tt.run_emcee(xdata,ydata,pixi[ind],pixj[ind])
        samples, samples_flat, results, results_err_up, results_err_low = tt.run_emcee(xdata,ydata,pixi[ind],pixj[ind])
        
        if make_plots == True:
            # Corner plot
            tt.plot_corner(samples_flat,pixi[ind],pixj[ind])
            # Trace plot of chains
            tt.plot_chains(samples,pixi[ind],pixj[ind])
            
        
                     
       
