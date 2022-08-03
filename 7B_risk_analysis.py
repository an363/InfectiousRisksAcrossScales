#! coding: utf-8 
from toolbox import *
import matplotlib as mpl

##########################
# This script takes as input the Pickle dump file
# created in 2A (containing the risks between 2 pedestrians)
# and plots risks of transmission.
# W. Garcia, A. Nicolas (2020)
##########################


# run with:
# python3 7B...py (num_scenario) (vx_wind) (vy_wind) (0/1 if isotropic inhalation)

######### FILES #########

group_from_scenario={"Pont_Morand": 0,
 					 "Perrache_Gare": 1, 
 					 "Trottoir_Croix_Rousse": 2, 
 					 "Passerelle_Bouchut": 3, 
 					 "Saint_Jean": 4,
 					 "Centre_Depistage_Gerland": 10,
 					 "Part_Dieu_Gare": 20, 
 					 "Bellecour": 21,
 					 "Croix_Rousse_Marché": 30,
 					 "Terreaux_Terrasses": 40,
 					 "Croix_Rousse_Terrasses": 41}
 
group_names= {1.3: "Streets",
			  4.5: "Queue",
			  6.6: "Stations",
			  8.8: "Market",
			  11.0: "Cafés"
			  }

full_names= { "Pont_Morand": "Rhône riverbank",
 			  "Perrache_Gare": "Forecourt of train station", 
 			  "Trottoir_Croix_Rousse": "Shopping street in Croix-Rousse",
 			  "Passerelle_Bouchut":	 "Street near train station", 
 			  "Saint_Jean": "Street in old town", 
 			  "Centre_Depistage_Gerland": "Queue at a Covid testing site",
 			  "Part_Dieu_Gare": "Lyon central train station", 
 			   "Bellecour": "Platform in metro station",
 			   "Croix_Rousse_Marché": "Market on Croix-Rousse boulevard",
 			   "Terreaux_Terrasses": "Street café in front of Lyon City Hall",
 			   "Croix_Rousse_Terrasses": "Street café on Croix-Rousse boulevard"
 			   }



#Get colormap as matrix of colors (to avoid white as first color)
cmaps = [  mpl.cm.Greys(np.linspace(0,1,7))[2:,:-1], mpl.cm.Blues(np.linspace(0,1,6))[2:,:-1], mpl.cm.Greens(np.linspace(0,1,6))[2:,:-1], mpl.cm.Purples(np.linspace(0,1,6))[2:,:-1], mpl.cm.Oranges(np.linspace(0,1,6))[2:,:-1] ]

couleurs= { scenario: cmaps[group_from_scenario[scenario]//10][group_from_scenario[scenario]%10] for scenario in Scenarios.values()}

abscisses= {scenario_number:  0.5 * ( group_from_scenario[Scenarios[scenario_number]]//10) + scenario_number for scenario_number in Scenarios.keys()}


static_scenarios= ["Centre_Depistage_Gerland", "Terreaux_Terrasses", "Croix_Rousse_Terrasses"]

Risks= {}
for scenario in Scenarios.values():
	#Risks[scenario]= pd.read_pickle(GeneralFolder+scenario+"/RisksEXP_InfAmidGroups.pkl") if (ContagionAmidGroups and scenario not in static_scenarios) else pd.read_pickle(GeneralFolder+scenario+"/RisksEXP.pkl") 
	#thisfilename= "/Risks/Risks_%.1f_%.1f%s.pkl"%(v_wind[0],v_wind[1],Activity)
	thisfilename= "/Risks/Risks_%.1f_%.1f%s%s.pkl"%(v_wind[0],
																	v_wind[1],
																	Activity if Activity!="_breathing" else "", 
																	"_iso" if iso_inh else "_aniso")
	print(thisfilename)
	Risks[scenario]= pd.read_pickle(GeneralFolder+scenario+thisfilename)
	for model in model_parameters.values():
		Risks[scenario]["Ci_%s"%model.name]= 0.5 * (Risks[scenario]["Clowi_%s"%model.name] + Risks[scenario]["Cbari_%s"%model.name])

### ###

print("Pedestrian can%s infect members of their own group"%("" if ContagionAmidGroups else "not"))

### Uncomment the desired section


### MEAN INFECTION RATE PLOTS

for model in model_parameters.values():
	fig = plt.figure(figsize=(4,4.5))
	model_name=  model.name
	T0_model= model.T0
	print(" * Analysing scenarios with model %s"%model_name)
	ax = fig.add_subplot(111)
	#ax.set_ylim(0,0.021)
	for scenario_number in Scenarios.keys():
		scenario= Scenarios[scenario_number]
		ax.bar( abscisses[scenario_number], Risks[scenario]["Ci_%s"%model_name].mean() * T0_model / DELTAT,
						yerr= 0.5*(Risks[scenario]["Cbari_%s"%model_name].mean() - Risks[scenario]["Clowi_%s"%model_name].mean()) * T0_model / DELTAT,
						width=1,
						ecolor="gray",
						error_kw=dict(lw=2, capsize=4, capthick=2),
						color= couleurs[scenario],
						label= "%s"%full_names[scenario])
		#!! Formerly: ax.bar( abscisses[scenario_number], Risks[scenario]["Ci_%s"%model_name].mean(),

	#Plotting
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)



	#!! ax.set_ylabel("New infections (/hour)",fontsize=14)
	ax.set_ylabel(r"Infection rate $\times T_0$",fontsize=14)
	plt.xticks( [x for x in group_names.keys()], [y for y in group_names.values()] ,fontsize=11, rotation=36)
	plt.tick_params(bottom=False)
	#ax.xticks.set_visible(False)
	#plt.legend( loc= 0, fancybox=True, framealpha=1, shadow=True, borderpad=1)
	plt.tight_layout()
	os.chdir(GeneralFolder)
	#ax.set_title("Model: %s"%model_name,fontweight= "bold")
	#!! plt.savefig(GeneralFolder+ ("histo_%s_InfAmidGroups.png"%model_name if ContagionAmidGroups else "histo_%s.png"%model_name), dpi=300)
	plt.savefig(GeneralFolder+ ("histoCFD_%s%s_InfAmidGroups_T0.png"%model_name if ContagionAmidGroups else "histoCFD_%.1f_%.1f%s%s.png"% (v_wind[0],
				  v_wind[1],
				  Activity,
				  "_iso" if iso_inh else "_aniso")), dpi=300) 
	plt.close(fig)

quit()


####### BOX PLOTS #######

for model in model_parameters.values():
	model_name= model.name
	T0_model= model.T0
	print(" * Analysing scenarios with model %s"%model_name)

	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	bplot= ax.boxplot( [ Risks[scenario]["Ci_%s"%model_name].values  * T0_model / DELTAT for scenario in Scenarios.values() ],
				       labels= [scenario for scenario in Scenarios.values()], 
				       positions= abscisses.values(),
				       meanprops= dict(linestyle=':', linewidth=4, color='firebrick'), 
				       medianprops= dict(linestyle='-', linewidth=2, color='black'),
				       boxprops= dict( linewidth=2, facecolor='white', alpha=0.5),
				       meanline= True, 
				       showmeans=True, 
				       patch_artist=True)
				       
	ax.set_ylabel(r"$C_i^{(T_0)}$",fontsize=15)

	for box, whisker, cap, flier, couleur in zip(bplot['boxes'],bplot['whiskers'],bplot['caps'],bplot['fliers'], couleurs.values()):
		box.set_color(couleur)
		flier.set_markeredgecolor(couleur)

	plt.xticks( [x for x in group_names.keys()], [y for y in group_names.values()] ,fontsize=11, rotation=36)
	ax.tick_params(axis='y', which='major', labelsize=11)
	plt.tick_params(bottom=False)

	plt.tight_layout()
	plt.savefig(GeneralFolder+ ("boxplotCFD_%s%s_InfAmidGroups_T0.png"%model_name if ContagionAmidGroups else "boxplotCFD_%.1f_%.1f%s%s.png"% (v_wind[0],
				  v_wind[1],
				  Activity,
				  "_iso" if iso_inh else "_aniso")), dpi=300) 
	plt.close(fig)


quit()



