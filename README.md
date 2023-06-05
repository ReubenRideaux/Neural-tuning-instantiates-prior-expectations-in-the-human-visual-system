# Neural tuning instantiates prior expectations in the human visual system

The MATLAB scripts/functions stored in this repository are a Key Resource for:

Harrison, Bays & Rideaux (2022) Neural tuning instantiates prior expectations in the human visual system. *bioRxiv*

This code can be used to analyze/generate data, in order to reproduce the results presented in the paper.

**The empirical data can be downloaded from: https://osf.io/5ba9y/**

### Dependencies

To run the analyses in this repository you will need the following dependencies:

Pim Mostert's decoding toolbox from: https://github.com/Pim-Mostert/decoding-toolbox

### Instructions

#### Analyzing empirical data and plotting results

(1) Download the empirical data to a folder ('data') within the same location as the 'empirical_eeg_data_analysis.m' script.

(2) Run the 'empirical_eeg_data_analysis.m' script.

Note: expect the analysis to take a while (~30 mins) depending on your machine.

(3) Run the 'empirical_figure_analyses.m' script to decode orientations and plot results.

Note: expect the analysis to take a while (~5 mins) depending on your machine.

#### Generative modelling

Run the 'generative_modelling.m' script to perform the simulate & analyze neural data, and plot the results.