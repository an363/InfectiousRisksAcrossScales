# InfectiousRisksAcrossScales
This repository contains the scripts for risk assessments of viral spread in macroscopic crowds, anchored in microscopic simulations.

Run:
* 7A to assess the risks in dynamic scenarios:
python3 7A... num_scenario vx_wind vy_wind isotropic_inhalation
* 7D to assess the risks in static scenarios:
python3 7D... num_scenario vx_wind vy_wind isotropic_inhalation
* 7B to plot the histograms of transmission risks across scenarios:
python3 7B...py (num_scenario) (vx_wind) (vy_wind) (0/1 if isotropic inhalation)
