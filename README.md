# InfectiousRisksAcrossScales
This repository contains the scripts for risk assessments of viral spread in macroscopic crowds, anchored in microscopic simulations.
******************************************************
******************************************************
The content can be used freely. Please cite at least one of the first two papers below in any work using these scripts:
* NICOLAS, Alexandre & MENDEZ, Simon. Viral transmission in pedestrian crowds: Coupling an open-source code assessing the risks of airborne contagion with diverse pedestrian dynamics models. Collective Dynamics (2024),  https://doi.org/10.17815/CD.2024.159
​https://doi.org/10.17815/CD.2024.159, https://www.collective-dynamics.eu/index.php/cod/article/view/A159
* MENDEZ, Simon, GARCIA, Willy, et NICOLAS, Alexandre. From Microscopic Droplets to Macroscopic Crowds: Crossing the Scales in Models of Short-Range Respiratory Disease Transmission, with Application to COVID-19, Advanced Science (2023), https://doi.org/10.1002/advs.202205255
* GARCIA, Willy, MENDEZ, Simon, FRAY, Baptiste, et al. Model-based assessment of the risks of viral transmission in non-confined crowds. Safety science (2021), vol. 144, p. 105453, https://www.sciencedirect.com/science/article/pii/S0925753521002964
******************************************************
******************************************************

0) INSTALLATION / DEPENDENCIES

Clone the GitHub repository ($ git clone ...) or download the Python scripts directly.
These scripts only require Python3, with the following fairly common modules: numpy, pandas, pickle, matplotlib, ...

1) UNZIP CONCENTRATION MAPS

 Uncompress the Diagrams.zip file into a new folder Diagrams.
 You can also uncompress Example_Data.zip into a new folder Example_Data if you want to use our pedestrian trajectory database

2) COLLECT PEDESTRIAN TRAJECTORIES

Create a folder to store the pedestrian trajectories; the folder name should end with the suffix "_static" for a static scenario, or "_dynamic" for a dynamic scenario (see the examples in "Example_data").
Within this folder, create one text file for each pedestrian's trajectory and name it [GROUP-ID]_[PED-ID].csv, where
* [GROUP-ID]: an integer label for the group to which the pedestrian belongs 
* [PED-ID]: an integer label for the pedestrian 
E.g., 1_3.csv, for pedestrian 3 belonging to social group 1.

	Each file should contain the following columns, with no header
	
time;x;y;theta

where
* time: timestamp in seconds associated with this line
* x: x-position (in metres) of the pedestrian at that time
* y: y-position (in metres) of the pedestrian at that time
* theta: angle (in rad) denoting the orientation of the pedestrian's head at that time, in the (x,y) frame. theta=0 corresponds to the +x-axis.

NB: Example data are proposed in the "Example_Data" directory. These correspond to the data of the 2021 paper referenced below; further information can be found in the paper and in the online repository https://zenodo.org/record/4527462
Should you want to recover the transmission rate values provided in our paper for these field data, do not forget to correct the size of the (narrow) empirical field of view on line 29 of toolbox.py

3) RUN RISK ASSESSMENT SCRIPTS

* Amend the input file "InputFile.txt" : 
--------------------------------------
DiagramsFolder= [FULL PATH TO FOLDER WITH SPATIO-TEMPORAL DIAGRAMS]

OutputFolder= [FULL PATH TO OUTPUT FOLDER]

TrajectoryFolder= [FULL PATH TO FOLDER WITH TRAJECTORIES]

T0= [characteristic infection time, in seconds]

(vx,vy)= (v_x,v_y) in m/s # external wind speed

ExhalationMode= [choose between: breathing / speaking / large\_droplets]

IsotropicInhalation= [choose between: True / False] # Can agents can inhale aerosols coming hitting the back of their heads? By default, we recommend to set it to False

ContagionAmidGroups= [choose between: True / False] # Can infect other members of their social group?

--------------------------------------

Note that semi-colons can be used to separate multiple input conditions if the user wants to launch multiple runs sequentially with a single input file, for instance (vx,vy)=(0.0,0.0);(-0.2,0.5) or ExhalationMode=speaking;breathing.

* Run the script with Python 3, e.g., by typing in a terminal

$ python3 main.py

4) READ THE RESULTS

The results in terms of infection rates caused by an individual in the OutputFolder repository, with one folder for each set of conditions. Each folder contains a summary of the parameters (parameters.txt), a file detailing how many new cases each distinct agent (ped_ID) would cause per hour, should they be contagious (Risks_by_person_output...dat}), and one file containing the mean number of new cases per hour (Risks_mean_output...dat). Within each file, the columns Clow and Cbar refer to:
* "Clow": lower bound on the assessed mean rate of new infections (new cases per hour)
* "Chigh": upper bound on the assessed mean rate of new infections (new cases per hour), taking into account possible previous interactions between pedestrians, out of the field of view
