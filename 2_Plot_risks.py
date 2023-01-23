#! coding: utf-8 
import matplotlib as mpl
from matplotlib import pyplot as plt
import pickle
import pandas as pd

##########################
# This script takes as input the Pickle dump file
# created in 1_ (containing the risks between 2 pedestrians)
# and plots risks of transmission.
# W. Garcia, A. Nicolas (2020), A. Nicolas (2022)
##########################


# RUN THIS SCRIPT USING
# $ python3 2_Plot_risks.py (num_scenario) (vx_wind) (vy_wind) (0/1 if isotropic inhalation)


######### FILES #########

# Link to the Pickle files where risks have been assessed and saved by the 1_ scripts
scenario_files= {0: "path/to/Risks_2.0_1.0_speaking_aniso.pkl",
				1: "path/to/Risks_2.0_0.0_speaking_aniso.pkl" }

# User-friendly labels for the different scenarios
full_names= { 0: "Rhône riverbank",
 			  	  1: "Queue at screening centre"
 			   }

##################
model_parameters= {0 : {"name": "standard", "T0": 900}}
DELTAT= 3600. # total time considered for risks of infection
##################

Risks= {}
for scenario in scenario_files.keys():
	Risks[scenario]= pd.read_pickle(scenario_files[scenario])
	for model in model_parameters.values():
		print(Risks[scenario])
		Risks[scenario]["Ci_%s"%model["name"]]= 0.5 * (Risks[scenario]["Clowi_%s"%model["name"]] + Risks[scenario]["Cbari_%s"%model["name"]])
		

### Uncomment the desired section

"""
### MEAN INFECTION RATE PLOTS
for model in model_parameters.values():
	fig = plt.figure(figsize=(4,4.5))
	model_name=  model["name"]
	T0_model= model["T0"]
	print(" * Analysing scenarios with model %s"%model_name)
	ax = fig.add_subplot(111)
	for scenario in scenario_files.keys():
		ax.bar( float(scenario), Risks[scenario]["Ci_%s"%model_name].mean() * T0_model / DELTAT,
						yerr= 0.5*(Risks[scenario]["Cbari_%s"%model_name].mean() - Risks[scenario]["Clowi_%s"%model_name].mean()) * T0_model / DELTAT,
						width=1,
						ecolor="gray",
						error_kw=dict(lw=2, capsize=4, capthick=2),
						label= "%s"%full_names[scenario])

	#Plotting
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)


	ax.set_ylabel(r"Infection rate $\times T_0$",fontsize=14)
	plt.xticks( [x for x in full_names.keys()], [y for y in full_names.values()] ,fontsize=11, rotation=36)
	plt.tick_params(bottom=False)
	#ax.xticks.set_visible(False)
	#plt.legend( loc= 0, fancybox=True, framealpha=1, shadow=True, borderpad=1)
	plt.tight_layout()
	ax.set_title("Model: %s"%model_name,fontweight= "bold")
	plt.show()
"""


####### BOX PLOTS #######

for model in model_parameters.values():
	model_name= model["name"]
	T0_model= model["T0"]
	print(" * Analysing scenarios with model %s"%model_name)

	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	bplot= ax.boxplot( [ Risks[scenario]["Ci_%s"%model_name].values  * T0_model / DELTAT for scenario in scenario_files.keys() ],
				       labels= [name for name in full_names.values()], 
				       positions=  [float(x) for x in scenario_files.keys()],
				       meanprops= dict(linestyle=':', linewidth=4, color='firebrick'), 
				       medianprops= dict(linestyle='-', linewidth=2, color='black'),
				       boxprops= dict( linewidth=2, facecolor='white', alpha=0.5),
				       meanline= True, 
				       showmeans=True, 
				       patch_artist=True)
				       
	ax.set_ylabel(r"$C_i^{(T_0)}$",fontsize=15)


	plt.xticks( [x for x in full_names.keys()], [y for y in full_names.values()] ,fontsize=11, rotation=36)
	ax.tick_params(axis='y', which='major', labelsize=11)
	plt.tick_params(bottom=False)

	plt.tight_layout()
	plt.show()
