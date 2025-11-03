# drift_controls_kelp
This repository contains all (draft) code, data, and figures associated with "_Sea urchin consumption of kelp controlled by the density of drift algae_". 

## Overview 
Kelp forests can rapidly shift between diverse canopy-forming ecosystems and urchin-dominated barrens when grazing pressure intensifies. 
This study experimentally quantified how purple urchin (_Strongylocentrotus purpuratus_) consumption of giant kelp (_Macrocystis pyrifera_) and detached drift varies with their relative availability. 
We found that urchins strongly prefer drift when it is abundant, but switch to grazing live kelp once drift becomes scarce. 
Quantitative thresholds indicate that maintaining at least ~1 g of drift per 1 g of kelp effectively prevents grazing on attached kelp. 
These results provide practical guidance for managers considering the use of drift subsidies to stabilize or restore kelp forests.

## data :open_file_folder:
- `drift_kelp_loss.csv` -- the raw original experimental data as entered from field data sheets.
- `juvenile_kelp_SNI_subtidal_timeseries.csv` -- the San Nicolas Island subtidal time series data from [Randell et. al., 2021](https://www.pnas.org/doi/10.1073/pnas.2103483119) used to visualize juvenile giant kelp individuals (Figure S2). 

## code :open_file_folder:
- `RunMe.R` -- control script used to select Logistic vs van Leeuwen model, read libraries, and `source()` to invoke and run all other scripts
- `analyze_visualize.R` -- visualize posteriors (Figure S7) and perform preference calculations (Figures 4, S9)
- `empirical.R` -- visualize empirical data (Figures 2, S3, and S4)
- `fit_stan_model.R` -- use cmdstan to compile and run stan models in R
- `format_data.R` -- format empirical data; prep Rdata file for stan
- `plot_juvenile_kelp.R` -- visualize juvenile kelp from _in situ_ monitoring to contrast with experimental densities 
- `plot_simulation.R` -- plot simulated temporal dynamics (Figures 3 and S8)
- `simulate.R` -- use stan posteriors to simulate temporal dynamics
- `simulate_process.R` -- extract values from simulations in preparation for visualizing Figures 3 and S8
- `stan_model_Logistic.stan` -- stan model containing the Logistic formulation 
- `stan_model_vanLeeuwen.stan` -- stan model containing the van Leeuwen formulation

## reproducibility
software versions utilized in these analyses: 
- Stan `v.2.32.2`
- _R_ `v.4.4.3`
- cmdstanr `v.0.7.1`

## 
<p align="center">
  <img src="https://github.com/user-attachments/assets/05915d15-01e9-4040-a498-43e06553a643" alt="stockingCage" width="40%" style="object-fit: cover; height: 300px;">
  <img src="https://github.com/user-attachments/assets/0890c64f-3d33-4b6a-ab79-9c2dcacb8aa5" alt="freshCage" width="40%" style="object-fit: cover; height: 300px;">
</p>
<p align="center">
  <img src="https://github.com/user-attachments/assets/4e5ca068-34d1-4a35-8a68-d37d8f017b2a" alt="tuckedAway" width="40%" style="object-fit: cover; height: 300px;">
  <img src="https://github.com/user-attachments/assets/37f7daf0-5d2c-44ed-8124-5042179d1f50" alt="activeGrazing" width="40%" style="object-fit: cover; height: 300px;">
</p>

