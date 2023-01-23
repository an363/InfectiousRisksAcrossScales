# InfectiousRisksAcrossScales
This repository contains the scripts for risk assessments of viral spread in macroscopic crowds, anchored in microscopic simulations.
******************************************************
******************************************************
The content can be used freely. Please cite the following papers in any work using these scripts:
* MENDEZ, Simon, GARCIA, Willy, et NICOLAS, Alexandre. From microscopic droplets to macroscopic crowds: Crossing the scales in models of short-range respiratory disease transmission, with application to COVID-19. arXiv preprint arXiv:2208.03147, 2022.
* GARCIA, Willy, MENDEZ, Simon, FRAY, Baptiste, et al. Model-based assessment of the risks of viral transmission in non-confined crowds. Safety science, 2021, vol. 144, p. 105453, https://www.sciencedirect.com/science/article/pii/S0925753521002964
******************************************************
******************************************************

0) INSTALLATION / DEPENDENCIES

Clone the GitHub repository ($ git clone ...) or download the Python scripts directly.
These scripts only require Python3, with the following fairly common modules: numpy, pandas, pickle, matplotlib, ...

1) UNZIP CONCENTRATION MAPS

 Uncompress the Diagrams.zip file into a new folder Diagrams.
 You can also uncompress Example_Data.zip into a new folder Example_Data if you want to use our pedestrian trajectory database

2) COLLECT PEDESTRIAN TRAJECTORIES

Create a folder to store the pedestrian trajectories in files titled (group_no)_(ped_no).csv, where
* group_no: an integer label for the group to which the pedestrian belongs 
* ped_no: an integer label for the pedestrian 

	Each file is specific to one pedestrian and should contain the following columns, with no header
	
time;x;y;theta

where
* time: timestamp in seconds associated with this line
* x: x-position (in metres) of the pedestrian at that time
* y: y-position (in metres) of the pedestrian at that time
* theta: angle (in rad) denoting the orientation of the pedestrian's head at that time, in the (x,y) frame

NB: Example data are proposed in the "Example_Data" directory. These correspond to the data of the 2021 paper referenced below; further information can be found in the paper and in the online repository https://zenodo.org/record/4527462

3) RUN RISK ASSESSMENT SCRIPTS

* Amend the Python script "toolbox.py" by linking the adequate input files (FILES section, lines 11-12) and adjusting the model parameters if need be (MODEL PARAMETERS, lines 15-17)
* Run the relevant script (scenario with moving people: 1A_Assess_risks_dynamic.py, mostly static crowd: 1B_Assess_risks_static.py):

$ python3 1A_Assess_risks_dynamic.py vx_wind vy_wind isotropic_inhalation

or

$ python3 1B_Assess_risks_static.py vx_wind vy_wind isotropic_inhalation

where
	* vx_wind: wind speed along the x-direction (in m/s)
	* vy_wind: wind speed along the y-direction (in m/s)
	* isotropic_inhalation: 1 if inhalation is to be considered isotropic, 0 otherwise (anisotropic inhalation)

* If you want to plot the histograms of transmission risks across scenarios, run 2_Plot_risks

$ python3 2_Plot_risks.py

after amending the directory names in the FILES section (around line 20)

4) READ THE RESULTS

Besides the Pickle dump, the results in terms of mean infection rates are saved as CSV files under the "Risks" repository: 
* "Clow_bar": lower bound on the assessed mean rate of new infections (new cases per hour)
* "Chigh_bar": upper bound on the assessed mean rate of new infections (new cases per hour), taking into account possible previous interactions between pedestrians, out of the field of view


