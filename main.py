#! coding: utf-8 
from toolbox import *

##########################
# This script takes as input trajectories stored as CSV files and
# assesses the transmission risks with different models.

# REFERENCES
# *  MENDEZ, Simon, GARCIA, Willy, NICOLAS, Alexandre, From Microscopic Droplets to Macroscopic Crowds: Crossing the Scales in Models of Short-Range Respiratory Disease Transmission, with Application to COVID-19. Adv. Sci. 2023, 2205255. https://doi.org/10.1002/advs.202205255 
# * GARCIA, Willy, MENDEZ, Simon, FRAY, Baptiste, NICOLAS, Alexandre. Model-based assessment of the risks of viral transmission in non-confined crowds. Safety science, 2021, vol. 144, p. 105453, https://www.sciencedirect.com/science/article/pii/S0925753521002964
##########################


# RUN THIS SCRIPT USING
# $ python3 main.py


######### FILES #########
InputFile= "InputFile"
#########################




Infection_time_T0= 900 # characteristic time (in seconds) for viral infection if somebody is facing you at 50 cm distance and talking with you 




def MainIt():

	# Reading parameter sets from input file
	(nb_parameter_sets, params)= read_params(InputFile)

	list_param_sets= [ {"TrajectoryFolder": folder, "DiagramsFolder": params["DiagramsFolder"][0], "OutputFolder": params["OutputFolder"][0]} for folder in params["TrajectoryFolder"]]
	# only one possible folder for Diagrams and Ouput
	
	for header in headers:
		if header in ["DiagramsFolder", "TrajectoryFolder"]:
			continue
		
		dimx= len(list_param_sets)
		dimy= len(params[header])
		matrice= np.full((dimx,dimy), {})
		
		for cpt_ps in range(dimx):
			for cpt_h in range(dimy):
				matrice[cpt_ps,cpt_h]= list_param_sets[cpt_ps] | {header: params[header][cpt_h]}
		list_param_sets= np.array(matrice).flatten()

	# CHECKED: the above works okay!
	
	scenarios= {}
	
	# Running script
	for param_set in list_param_sets:
		foldername= param_set["OutputFolder"]+"/Set%i"%time.time()
		try:
			os.makedirs(foldername)
		except:
			pass
			
		with open(foldername+"/parameters.txt",'w') as monfichier:
			mystring="\n Currently analysing case "
			for key, value in param_set.items():
				mystring+= "%s : %s // "%(key, value)
				monfichier.write("%s= %s\n"%(key,value))
			print(mystring)
		
		### load relevant values
		### LOADING DYNAMIC CONCENTRATION MAPS ###
		if param_set["ExhalationMode"] not in spatiotemp_diagrams_all.keys():
			spatiotemp_diagrams_all[ param_set["ExhalationMode"] ]= load_diagrams(param_set) # Loads the dynamics viral concentration maps
			Z_all[  param_set["ExhalationMode"] ]= 960.79508 if (param_set["ExhalationMode"]=="large_droplets") else 1.253489635
		### Renormalisation coefficient Z ###
		
		if param_set["TrajectoryFolder"] not in scenarios.keys():
			scenarios[ param_set["TrajectoryFolder"] ]= getDetections( param_set["TrajectoryFolder"] )

	
		scenario= scenarios[ param_set["TrajectoryFolder"] ]
		mypedestrians= scenario["mypedestrians"]
		
		if isStatic(param_set["TrajectoryFolder"]):
			RISKS= getStaticRisks(scenario, param_set)
			writeStaticResults(foldername,param_set["TrajectoryFolder"].split("/")[-1]+".dat", mypedestrians, RISKS)
		else:
			RISKS= getDynamicRisks(scenario, param_set)
			writeDynamicResults(foldername,param_set["TrajectoryFolder"].split("/")[-1]+".dat", mypedestrians, RISKS)


		
	

if __name__ == "__main__":
	MainIt()
