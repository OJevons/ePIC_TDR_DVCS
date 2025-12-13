# ePIC_TDR_DVCS
This repository contains the analysis and plotting code for the analysis of Deeply Virtual Compton Scattering events, simulated for the ePIC experiment. This analysis was performed as part of the Exclusive, Diffractive and Tagging physics working group within ePIC.

The creation of plots for the ePIC TDR will require two steps:
- running the maim analysis script in order to fill ROOT histograms, then
- plotting the created histograms.

The main analysis takes the form of an `ePIC_DVCS_TASK` object, user-defined in `ePIC_DVCS_TASK.h`. THsi specific form which this analysis object will take is defined in `ePIC_DVCS_TDR.cxx`. The file `run_ePIC_DVCS.C` is provided as a wrapper macro to run these analysis objects. This script does the following:
- Initialises a DVCS analysis object.
- Gives, as input to the object, a list of files which the analysis is to run over (standard filelists for the `25.10.2` campaign files for all three 'default' EIC energy settings are provided in this repo).
- Defines the output file nam, which will store the created histograms.
- Sets the behaviour with which to analyse the provided events.
- Tells the object to run its analysis.

The input file lists should point to the location in which DVCS simulation campaign files are stored on the ePIC filesystem hosted at Jefferson Lab. As such, the scripts provided in this repository **must** be run within the [`eic-shell`](https://github.com/eic/eic-shell) environment, using the ROOT analysis framework. In order to run the analysis macro, the syntax is
```
root 'run_ePIC_DVCSC("5x41/10x100/18x275")'
```
where the parameter inside the brackets declare which beam energy to run the analysis on. By default (if the brackets are left empty when executing the run macro) the analysis is set to run using the 10x100 GeV beam energy setting.



~The analysis script runs a light version of the full DVCS analysis over the provided files, and creates three plots as output (which are stored in the `figs` directory). These plots are:
- the track psuedorapidity distributions for all expected final state particles (scattered electron, proton and photon), for both generated and reconstructed particles,
- the distribution of the Mandelstam t variable for both generated and reconstructed DVCS events, on the condition that the full final state is reconstructed, and~
- the difference between the predicted and measured track theta for the detected DVCS photon, for all reconstructed photons (as long as only one photon is identified).
