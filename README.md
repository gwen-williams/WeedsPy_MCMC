# WeedsPy_MCMC

Repository still a work in progress!

## Dependencies

WeedsPy_MCMC requires the installation of the GILDAS package, instructions for which can be found [here](https://www.iram.fr/IRAMFR/GILDAS/gildasli2.html). 

These scripts have been tested in Python 3.9.9 only.

The following Python packages are required:

* numpy
* astropy
* emcee
* corner
* matplotlib
* colorama

For Mac and Anaconda users: CLASS' preferred install method on Mac is through MacPorts. These scripts utilise CLASS' ability to talk to Python -- CLASS does not talk to the Anaconda-installed Python, nor the Anacaonda `pip`. It is advised that you have a MacPorts-installed Python as well, and use the system `pip` to install Python packages, otherwise these scripts will not run. To check whch `pip` is currently in use, type `which pip` into a terminal - if the path printed to the terminal is that to your local bin, then the system's pip is already in use. If the path is instead to the Anaconda bin, what worked for me was to temporarily comment out the Anaconda path in my `~/.bashrc` file to force the system `pip` to be used. Once the packages were installed I then went back and undid the comment.


## Script methodology

The estimation of the physical parameters of an astromonical source, such as column density and temperature, is possible with the modelling of a synthetic spectrum assuming LTE using the Weeds extension (Maret et al. 2011) of the CLASS software (part of the GILDAS package). The spectrum is produced by assuming a certain column density, excitation temperature, source size, systemic velocity and FWHM velocity width, whilst the Einstein coefficients, partition function, and upper energy level degeneracy and energy are taken from either the JPL or CDMS catalogs.  Weeds however does not allow the exploration of parameter space efficiently, leaving the user to manually change the fitting parameters. The `WeedsPy_MCMC` scripts we present here use the in-built functionality of CLASS to talk to Python to facilitate a Bayesian analysis of the parameter space in the LTE modelling using the `emcee` Python package.

We include two generalised scripts, called `WeedsPy_MCMC.py` and `WeedsPy_MCMC.class`. Variables are set by the use in `paramfile.txt`, and are read in by the `read_param_file` function. We use the `emcee` Python package to explore the parameter space of the CLASS/Weeds LTE synthetic spectra. The scripts are currently set up to explore the parameter space of four free parameters: column density, excitation temperature, systemic velocity, and FWHM velocity width.


## User guide

Here we provide some guidance on the use of the scripts.


### Directories and files



### Parameter file


### Running the scripts

Start a CLASS terminal from the parent directory. To run the scripts, execute the following:

```
PYTHON WeedsPy_MCMC.py
```


## Future work

* Future versions of this script will include the option to parallelise the code.

* We currently rely on the user to create the `.30m` files of their observed spectra themselves. A future version of this script will include a function that will do this for the user.


## Citing the code 

These scripts were developed for the work presented in Williams et al. (2023, MNRAS, to be submitted). We kindly ask for its citation if you find these scripts helpful in your own work.


## Acknowledgements

These scripts were developed with support from the UK's Science and Technology Facilities Council under grant number ST/W00125X/1.

