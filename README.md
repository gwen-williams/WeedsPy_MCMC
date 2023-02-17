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


## Script methodology

The estimation of the physical parameters of an astromonical source, such as column density and temperature, is possible with the modelling of a synthetic spectrum assuming LTE using the Weeds extension (Maret et al. 2011) of the CLASS software (part of the GILDAS package). The spectrum is produced by assuming a certain column density, excitation temperature, source size, systemic velocity and FWHM velocity width, whilst the Einstein coefficients, partition function, and upper energy level degeneracy and energy are taken from either the JPL or CDMS catalogs.  Weeds however does not allow the exploration of parameter space efficiently, leaving the user to manually change the fitting parameters. The `WeedsPy_MCMC` scripts we present here use the in-built functionality of CLASS to talk to Python to facilitate a Bayesian analysis of the parameter space in the LTE modelling using the `emcee` Python package.

We include two generalised scripts, called `WeedsPy_MCMC.py` and `WeedsPy_MCMC.class`. Variables are set by the use in `paramfile.txt`, and are read in by the `read_param_file` function. We use the `emcee` Python package to explore the parameter space of the CLASS/Weeds LTE synthetic spectra. All function that carry out the emcee are located in a Python class object. Variables from the paramfile, and others set by the use in their main script, are passed to the class object. 

These scripts give the user the option of exploring either 2, 4 or 5 free parameters in the synthetic spectra. The parameters explored are as follows:
* two free parameters : column density and temperature
* four free parameters : column density, temperature, centroid velocity, velocity width
* five free parameters : column density, temperature, centroid velocity, velocity width, source size.

Missing parameters in these use-cases are fixed by the user.


## User guide

We include two scripts, called `WeedsPy_MCMC.py` and `WeedsPy_MCMC.class`. 


### Directories and files

Place the `WeedsPy_MCMC.py` and `WeedsPy_MCMC.class` scripts, and the parameter file `paramfile.txt` in the same parent directory. Within this directory, create the following sub-directories:

* `mdlfiles/` : to contain the `.mdl` files containing the model parameters that CLASS reads to conduct the modelling.
* `synthfiles/` : to contain the `.30m` and `.fits` files of the highest likelihood synthetic spectra.
* `spectra/` : to contain your observational data files.
* `plots/` : to contain the `.pdf` plots of the emcee trace plot and corner plot.
* `samples/` : to contain `.csv` files that store the flat chains of the emcee sampler.
* `theta_max/` : to contain `.csv` files that store the highest likelihood values.


### Parameter file

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
* `freq_unit` : the unit of the frequency, either MHz, GHz, or Hz.

Other parameters that require significant tweaking should be set within the body of your main Python scripts, rather than the parameter file. These include the values of any fixed parameters in your modelling, the value range of your priors, and the initial location of your walkers. Examples of how to do this are shown in the `example_scripts` sub-directory.


### Observed spectra

The scripts are currently setup to require you to place your observed spectra in the `spectra/` sub-directory. Currrently, your spectra must already be in the CLASS file format of `.30m` (the GILDAS/CLASS handbook details how to do the conversion), and your filename must follow the convention `spec_i100_j200.30m`, where `i100` refers to the x pixel coordinate of 100, and `j200` refers to the y pixel coordinate of 200.



### Running the scripts

Start a CLASS terminal from your parent directory. To run the scripts, execute the following:

```
PYTHON WeedsPy_MCMC.py
```

All functions to carry out the `emcee` sampling are placed inside a Python class within the `WeedsPy_MCMC.py` Python script. Unfortunately, due to the nesting of the Gildas/Python `Sic` commands, any calls to the Gildas/CLASS terminal with `Sic` will not be recognised if placed in a secondary Python script that is then called by a primary Python script. As such, the Python class object that contains all the functions to carry out the analysis must be placed in the preamble of your main Python script. Example scripts are shown in the `example_scripts` directory of how you should structure your analysis script.


## Future work

* The biggest bottle-neck with this code is the run time. We are currently looking into parallelisation, and hope a future version will include an option for this if possible.

* We currently rely on the user to create the `.30m` files of their observed spectra themselves. A future version of this script may include a function to convert `.fits` files to `.30m` files.


## Citing the code 

These scripts were developed for the work presented in Williams et al. (2023, MNRAS, to be submitted). We kindly ask for its citation if you find these scripts helpful in your own work.


## Acknowledgements

These scripts were developed with support from the UK's Science and Technology Facilities Council under grant number ST/W00125X/1.

