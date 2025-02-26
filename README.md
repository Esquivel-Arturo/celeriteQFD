# Stellar Flares Detection From Photometric Data Using Hidden Markov Models

This repository contains our work used for our [paper](https://doi.org/10.48550/arXiv.2404.13145](https://iopscience.iop.org/article/10.3847/1538-4357/ad95f6/meta) on how to use a Gaussian Process and a hidden Markov model (HMM) to conduct flare detection. We call our model `celeriteQFD`. It simultaneously implements detrending and flare detection. We use `Celerite`, a Gaussian Process adopted from Python package [exoplanet-celerite2](https://github.com/exoplanet-dev/celerite2), for the former, and a three-state HMM for the latter. The R script [celeriteQFD](https://github.com/Esquivel-Arturo/celeriteQFD/blob/main/Res/CeleriteQFD/031381302/celeriteQFD.R) can be used to fit the model and produce most of the plots shown in the paper. As part of the analysis, Kepler flares were simulated using [simuFlare.R](https://github.com/Esquivel-Arturo/celeriteQFD/blob/main/R/simuFlares.R) and injected into the light curve of M dwarf TIC 031381302. Multiple injection and recovery experiments were performed using [flare_injection_CHTC.R](https://github.com/Esquivel-Arturo/celeriteQFD/blob/main/Res/Injection_recover/031381302/flare_injection_CHTC.R). An example of one experiment, including a comparison of our method with 1-3$`\sigma`$, is shown below:

![Alt text](https://github.com/Esquivel-Arturo/celeriteQFD/blob/main/Res/Injection_recover/flare_inj.png)

Performance metrics for each experiment were computed using [combine_data.R](https://github.com/Esquivel-Arturo/celeriteQFD/blob/main/Res/Injection_recover/031381302/combine_data.R).

The code for this project was mostly developed by [Yunyi Shen](https://github.com/YunyiShen?tab=repositories). His Github repository contains further content of the development phase of the model. 

