# ePIC_TDR_DVCS
This repository contains the analysis and plotting code for the ePIC DVCS analysis (done as part of the Exclusive, Diffractive and Tagging PWG).

In order to run the analysis, one must provide a list of input files for the main analysis code (`DVCS_TDR.cxx`) to run over. The list provided in this repository points to the output of the ePIC simulation campaign for October 2024 (listed as the campaign to be used for the TDR physics simulation plots), using the 10x100 beam energy.

The scripts provided in this repository **must** be run within the [`eic-shell`](https://github.com/eic/eic-shell) environment, using the ROOT analysis framework. The syntax for running the analysis code is as follows:
```
root 'DVCS_TDR.cxx("$PATH/$TO/$FILELIST")'
```

The analysis script runs a light version of the full DVCS analysis over the provided files, and creates three plots as output (which are stored in the `figs` directory). These plots are:
- the track psuedorapidity distributions for all expected final state particles (scattered electron, proton and photon), for both generated and reconstructed particles,
- the distribution of the Mandelstam t variable for both generated and reconstructed DVCS events, on the condition that the full final state is reconstructed, and
- the difference between the predicted and measured track theta for the detected DVCS photon, for all reconstructed photons (as long as only one photon is identified).
