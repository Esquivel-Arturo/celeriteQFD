# Stellar Flares Detection From Photometric Data Using Hidden Markov Models

This repo contains our work, shown in ---, on how to use a Gaussian Process and a hidden Markov model (HMM) to conduct flare detection. We call our model `celeriteQFD`. It simultaneously implements detrending and flare detection. We use `Celerite`, a Gaussian Process adopted from Python package [exoplanet-celerite2](https://github.com/exoplanet-dev/celerite2), for the former, and a three-state HMM for the latter. The R script [celeriteQFD](https://github.com/Esquivel-Arturo/celeriteQFD/blob/main/Res/CeleriteQFD/031381302/celeriteQFD.R) can be used to fit the model and produce most of the plots shown in the paper.

![Alt text](https://github.com/Esquivel-Arturo/celeriteQFD/blob/main/Res/Injection_recover/flare_inj.png)
