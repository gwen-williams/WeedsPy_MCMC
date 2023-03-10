# WeedsPy_MCMC

Repository and README still a work in progress!

## Dependencies

WeedsPy_MCMC requires the installation of the GILDAS package, instructions for which can be found [here](https://www.iram.fr/IRAMFR/GILDAS/gildasli2.html). 

These scripts have been tested in Python 3.9.9, and GILDAS version nov21a, only.

The following Python packages are required:

* os
* numpy
* astropy
* emcee
* corner
* matplotlib

For users of Anaconda on Mac: CLASS' preferred install method on Mac is through MacPorts. These scripts utilise CLASS' ability to talk to Python -- CLASS does not talk to the Anaconda-installed Python, nor the Anacaonda `pip`. It is advised that you have a MacPorts-installed Python as well, and use the system `pip` to install Python packages, otherwise these scripts will not run.


## Overview of script methodology

The [Weeds extension](https://www.iram.fr/IRAMFR/GILDAS/doc/pdf/weeds.pdf) of the CLASS software (part of the Gildas package) is a radiative transfer code that produces synthetic emission spectra for various molecules assuming local thermodynamic equilibrium [(Maret et al. 2011)](https://ui.adsabs.harvard.edu/abs/2011A%26A...526A..47M/abstract). Model parameters that are set by the user are the column density, excitation temperature, centroid velocity, FWHM velocity width, and source size. Weeds however does not allow the user the efficiently explore this parameter space, leaving the user to manually change them. The `WeedsPy_MCMC` scripts we present here use the in-built functionality of Gildas/CLASS to talk to Python to facilitate a Bayesian analysis of the parameter space in the LTE synthetic spectra using the `emcee` Python package.

We include two generalised scripts, called `WeedsPy_MCMC.py` and `WeedsPy_MCMC.class`. Most variables are set by the user in `paramfile.txt`, and are read in by the `read_param_file` function. We use the `emcee` Python package to explore the parameter space of the synthetic spectra. All functions that carry out the emcee are located in a Python class in `WeedsPy_MCMC.py`.

These scripts give the user the option of exploring either 2, 4 or 5 free parameters in the synthetic spectra. The parameters explored in each case are as follows:
* two free parameters : column density and temperature
* four free parameters : column density, temperature, centroid velocity, velocity width
* five free parameters : column density, temperature, centroid velocity, velocity width, source size.

Parameters not listed as free in each of these use-cases are fixed by the user.


## User guide

We include two scripts, called `WeedsPy_MCMC.py` and `WeedsPy_MCMC.class`. Example scripts for different use cases are included in the `example_scripts/` directory.


### Directories and files

Place the `WeedsPy_MCMC.py` and `WeedsPy_MCMC.class` scripts, and the parameter file `paramfile.txt` in the same parent directory. Within this directory, create the following sub-directories:

* `mdlfiles/` : to contain the `.mdl` files containing the model parameters that CLASS reads to conduct the modelling.
* `synthfiles/` : to contain the `.30m` and `.fits` files of the highest likelihood synthetic spectra.
* `spectra/` : to contain your observational data files.
* `plots/` : to contain the `.pdf` plots of the emcee trace plot and corner plot.
* `samples/` : to contain `.csv` files that store the flat chains of the emcee sampler.
* `theta_max/` : to contain `.csv` files that store the highest likelihood values.


### Parameters

Most parameters required by the scripts are set in a parameter file, and read in to a dictionary by the `read_params_file` function within `WeedsPy_MCMC.py`. A template parameter file is given here called `paramfile.txt`. Here, we detail each of the required parameters:

* `catalog_name` : the name of the molecular line catalog for CLASS to use. Either the name of your own offline file, or cdms or jpl.
* `molecule` :  name of molecule to be modelled, as it appears in the catalog.
* `source_size` : size of source, in arcseconds.
* `mean_freq` : mean frequency of the spectral window, in MHz.
* `telescope_diameter` : diameter of telescope, in metres.
* `bmaj` : beam major axis, in arcseconds.
* `bmin` : beam minor axis, in arcseconds.
* `rms_noise` : RMS noise of the spectrum.
* `rms_noise_unit`:  unit of the rms noise, either mJy/beam, Jy/beam, uJy/beam or K.
* `nburn` : number of emcee burn iterations
* `nprod` :	number of emcee production iterations
* `nwalkers` : number of emcee walkers
* `freq_unit` : the unit of the frequency, either Hz, kHz, MHz or GHz.

Other parameters that (in my experience) are regularly tweaked are set within the body of your main Python scripts, rather than the parameter file. Examples of how to do this are shown in the `example_scripts` sub-directory, and each one is detailed here:

* Initial walker positions : this should be a list, N elements long, where N is the number of free parameters.
* The priors of your free parameters : this should be a list, Nx2 elements long. Each pair of numbers corresponds to the lower and upper bounds respectively of one of the free parameters. For example, for the case where there are 2 free parameters, 

	```
	priors = [0.5, 1.5, 100., 200.]
	```
	where `priors[0]` and `priors[1]` correspond to the lower and upper bounds of the column density, and `priors[2]` and `priors[3]` correspond to the lower and upper bounds of the temperature.

* Number of free parameters : This is set in the script as the length of the initial walker position list.
* Base power of the column density : Here you should set the base value of the column density. For example, `column_base = 1e18`. The `emcee` sampler is run on values of the column density / column_base, and the base is applied outside of `emcee`.
* Model parameters that are fixed : `source_size`, `vel_sys` and `vel_width` refer to the source size, centroid velocity, and wvelocity width. They should be set to `None` if you intend them to be free parameters, or set to your fixed (float) values if you intend them to be fixed parameters.

These include the values of any fixed parameters in your modelling, the value range of your priors, and the initial location of your walkers. Examples of how to do this are shown in the `example_scripts` sub-directory.


### Observed spectra

The scripts are currently setup to require you to place your observed spectra in the `spectra/` sub-directory. Currrently, your spectra must already be in the CLASS file format of `.30m` (the GILDAS/CLASS handbook details how to do the conversion), and must already to converted to brightness temperature (units of Kelvin). Your filename must follow the convention `spec_iXX_jYY.30m`, where `XX` and `YY` should be replaced by the integer values of the x and y pixel coordinates respectively.



### Saved results

In all of the following filenames, `XX` and `YY` refer to the integer value of the x and y pixel coordinate of the current spectrum respectively. 

The `make_gildas_mdl_file` function creates a CLASS-readable file that contains the name of the molecules to be modelled, the value of the parameters that the `emcee` sampler has walked to and (if any) the value of the fixed parameters. During the `emcee` iterations, these files are saved as temporary files that continually get overwritten, are placed in the `mdlfiles/` sub-directory, and are called `temp_mdlfile.mdl`. After the `emcee` has walked all of its iterations, the highest likelihood model files is made by the `save_highest_likelihood_spectrum` function and is saved as `best_mdlfile_iXX_jYY.mdl` where XX and YY are replaced by the integer value of the x and y pixel coordinates respectively.

The synthetic spectrum produced by CLASS at the parameters within the CLASS-readable .mdl files are saved in temporary files that continually get overwritten during the `emcee` iterations, and are saved in the `synthfiles/` sub-directory with filenames `temp_synthfile.fits` and `temp_synthfile.30m`. After the `emcee` has walked all of its iterations, the synthetic spectrum with the highest likelihood parameters are saved by the `save_highest_likelihood_spectrum` function to filenames `best_synthfile_iXX_jYY.fits` and `best_synthfile_iXX_jYY.30m`.

The highest likelihood model parameters are saved in the `theta_max/` sub-directory in files called `theta_max_flatlnprob_iXX_jYY.csv`. 

The full flatchain of the `emcee` sampler is saved in the `samples/` sub-directory in files called `samples_flatchain_iXX_jYY.csv`. 

Plots of the walker chains, the corner plot, and the spectrum with overplotted highest-likelihood synthetic spectrum are saved in the `plots/` sub-directory, and are called `chains_iXX_jYY.pdf`, `corner_iXX_jYY.pdf`, and `spec_iXX_jYY.eps` respectively.



### Running the scripts

Start a CLASS terminal from your parent directory. To run the scripts, execute the following:

```
PYTHON WeedsPy_MCMC.py
```

All functions to carry out the `emcee` sampling are placed inside a Python class within the `WeedsPy_MCMC.py` Python script. Unfortunately, due to the nesting of the Gildas/Python `Sic` commands, any calls to the Gildas/CLASS terminal with `Sic` will not be recognised if placed in a secondary Python script that is then called by the primary Python script you run in the CLASS terminal. As such, the Python class that contains all the functions to carry out the analysis must be placed in the preamble of your main Python script, rather than used as a separate module. Example scripts are shown in the `example_scripts` directory of how you should structure your analysis script.

The code is written to be run on a number of pixels consecutively, with a continuum brightness value per pixel. The x and y pixel coordinates should be saved in text files and read in in your script. Examples of this are shown in the `example_scripts/` sub-directory. The x and y pixel coordinates are currently required for the filenames of all files saved by the script (as detailed above). The `WeedsPy_MCMC_2params.py` scripts in the `example_scripts/` sub-directory shows an example script set-up for multiple pixels. The `WeedsPy_MCMC_4params.py` scripts in the `example_scripts/` sub-directory shows an example script set-up for one pixel.


## Future work

* These scripts are not parallelised. We are currently looking into parallelisation, and hope a future version will include an option for this if possible.

* We currently rely on the user to create the `.30m` files of their observed spectra themselves. The Gildas/CLASS manual details how to do this, but a future version of this code may include a function to convert `.fits` files to `.30m` files.

* Querying the online CDMS and JPL catalogs for every iteration of every walker is very slow indeed. The code runs significantly faster if querying an offline version of the molecular line catalog. Instructions on how to save an offline copy of a portion of either of the two online databases are given in the CLASS/Weeds manual. Using your own offline copy is supported here (by placing the filename in the `catalog_name` parameter of the paramfile), however we have not included a function to create the offline database copy in the first place. We aim to include a function to do this in a future version of this code.


## Citing the code 

These scripts were developed for the work presented in Williams et al. (2023, MNRAS, to be submitted). We kindly ask for its citation if you find these scripts helpful in your own work.


## Acknowledgements

These scripts were developed with support from the UK's Science and Technology Facilities Council under grant number ST/W00125X/1.

