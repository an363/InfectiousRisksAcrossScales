#! coding: utf-8 
from toolbox import *

##########################
# This script takes as input trajectories stored as CSV files and
# assesses the transmission risks with different models.
# It prints out the results and saves them in a CSV file and in a Pickle dump file.
# W. Garcia, A. Nicolas (2020), A. Nicolas (2022)

# REFERENCES
# * GARCIA, Willy, MENDEZ, Simon, FRAY, Baptiste, et al. Model-based assessment of the risks of viral transmission in non-confined crowds. Safety science, 2021, vol. 144, p. 105453, https://www.sciencedirect.com/science/article/pii/S0925753521002964
# * MENDEZ, Simon, GARCIA, Willy, et NICOLAS, Alexandre. From microscopic droplets to macroscopic crowds: Crossing the scales in models of short-range respiratory disease transmission, with application to COVID-19. arXiv preprint arXiv:2208.03147, 2022.
##########################


# RUN THIS SCRIPT USING
# $ python3 1B_Assess_risks_static.py vx_wind vy_wind isotropic_inhalation

	

print("Pedestrians can%s infect members of their own group"%("" if ContagionAmidGroups else "not"))
DumpFileRisks= FolderTraj + "/Risks" + "/Risks" + suffix + ".pkl" 
SaveRisksFile= FolderTraj + "/Risks" + "/risks" + suffix + ".csv"

if not os.path.isdir(FolderTraj + "/Risks"):
	os.mkdir(FolderTraj + "/Risks")
  


detections_dict= {}

index=0 
for fichier in os.listdir(FolderTraj):
	m= re.match("(\d+)_(\d+).csv",fichier)
	if m==None:
		continue
	group_no= int(m.group(1))
	ped_no= int(m.group(2))
	
	with open(FolderTraj+"/"+fichier,'r') as monfichier:
		for maligne in monfichier.readlines():
			nbs= maligne.strip("\n").split(";")
			detections_dict[ index ]= { "group_no": group_no, "ped_no": ped_no, "time": float(nbs[0]), "x": float(nbs[1]), "y": float(nbs[2]), "theta": float(nbs[3]) }
			index+= 1
			#print(detections_dict[ (group_no,ped_no) ])


detections= pd.DataFrame.from_dict(detections_dict, orient='index')
del detections_dict
#print(detections)



### Add velocities ###
detections["vx"]=0.0
detections["vy"]=0.0
for index, row in detections.iterrows():
	if index+1 < detections.shape[0] and detections.at[index+1,"ped_no"] == row["ped_no"]:
		delta_t= detections.at[index+1,"time"] - row["time"]
		detections.at[index,"vx"] = (detections.at[index+1,"x"] - row["x"]) / delta_t
		detections.at[index,"vy"] = (detections.at[index+1,"y"] - row["y"]) / delta_t
	elif index > 0 and detections.at[index-1,"ped_no"] == row["ped_no"]:
		delta_t= row["time"] - detections.at[index-1,"time"]
		detections.at[index,"vx"] = (row["x"] - detections.at[index-1,"x"]) / delta_t
		detections.at[index,"vy"] = (row["y"] - detections.at[index-1,"y"]) / delta_t
### ###
		
detections["ped_no"]= detections["ped_no"].astype(int)
# Change group labels (because there is a conflict between group numbers from different video sequences)
detections["group_no"]= detections.apply(lambda row: (row["ped_no"]//1000)*1000 + int(row["group_no"]), axis= 1)




times= sorted(detections["time"].unique())
dt= min([times[cpt+1] - times[cpt] for cpt in range(len(times)-1)]) # a greedy way to find the elementary time step
dt= round(dt,3)
print(" The time step is dt = ", dt)



####### Correct for EDGE EFFECTS #######
rescaling_factor= get_rescaling_factor(Lx,Ly)



####### PEDESTRIANS & GROUPS #######


# PEDESTRIANS #
group_from_ped= {} # dictionary such that group_from_ped[ myped ] = mygroup
for row in detections.itertuples():
	group_from_ped[ row.ped_no ] = row.group_no

mypedestrians= pd.DataFrame(data={'ped_no': np.array(sorted(group_from_ped.keys()))} )
print("Total number of pedestrians :",mypedestrians.shape[0])

mypedestrians["ped_no"]= mypedestrians["ped_no"].astype(int)

mypedestrians["tau"]= mypedestrians["ped_no"].apply( lambda ped_no: 0.5 + dt * (detections.loc[ detections["ped_no"]==ped_no ].shape[0] )  )  # time spent in the field of view, with an additional 0.5 s due to mis-reconstruction of beginning and end of trajectory
mypedestrians["group_no"]= mypedestrians["ped_no"].apply(lambda x: group_from_ped[x] )



####### TRANSMISSION RISKS BETWEEN INDIVIDUAL PEDESTRIANS #######
Risks_ij= {}


detections.sort_values(by=["ped_no","time"], inplace=True)



for parameters in model_parameters.values():
	
	# Prepare empty risks: risks_AB is the risk incurred by B (receiver) because of A (emitter)
	risks_AB= {}
	for ped_A in group_from_ped.keys():
		for ped_B in group_from_ped.keys():
			risks_AB[(ped_A,ped_B)]= 0.0
	
	last_hped= -1
	for row_A in detections.itertuples():
		if row_A.ped_no!=last_hped:
			last_hped= row_A.ped_no
			print(" Currently handling ped n.%i"%row_A.ped_no)
			
		# vel_A= row_A.vx # velocity of pedestrian A
		detections_loc= detections.loc[ (detections["ped_no"]!=row_A.ped_no) & (detections["time"]>=row_A.time) & (detections["time"]<row_A.time+tau_max) ]

		for row_B in detections_loc.itertuples():
			r_vec= (row_B.x-row_A.x, row_B.y-row_A.y)
			r= get_norm(r_vec)
			
			tau= row_B.time - row_A.time

			#!! For computational efficiency, discard people more than 4 metres away
			if r>4.0:
				continue
			v_walk= (row_A.vx, row_A.vy)
			
			
			risks_AB[(row_A.ped_no,row_B.ped_no)]+= dt * dt * rescaling_factor.get( int(round(10*r)), 1) * get_nu_dynamic(v_walk, v_wind, r_vec, row_A.theta, row_B.theta, tau, parameters)
			# the rescaling coefficient accounts, at the mean-field level, for the edge effects whereby distant contacts are not seen
			# dt transformed in dt**2 on July 25, 2022
			
				
			
			
	# Drop entries with zero risks
	for key in list(risks_AB.keys()):
		if risks_AB[key]<1e-18:
			del risks_AB[key]
	
	Risks_ij[parameters.name]= pd.DataFrame.from_dict(risks_AB, orient='index', columns=["sum_risk_AtoB"])
	Risks_ij[parameters.name].reset_index(inplace=True,drop=False)
	Risks_ij[parameters.name]["ped_no_A"]= Risks_ij[parameters.name]["index"].apply(lambda t: int(t[0]))
	Risks_ij[parameters.name]["ped_no_B"]= Risks_ij[parameters.name]["index"].apply(lambda t: int(t[1]))
	Risks_ij[parameters.name].drop(columns=['index'], inplace=True)

	
	
	Risks_ij[parameters.name]["group_no_A"]= Risks_ij[parameters.name]["ped_no_A"].apply(lambda x: group_from_ped[x])
	Risks_ij[parameters.name]["group_no_B"]= Risks_ij[parameters.name]["ped_no_B"].apply(lambda x: group_from_ped[x])
	
	# a diseased individual CAN infect members of his/her own group in these static scenarios
	Risks_ij[parameters.name]["exp_risk_AtoB"]= Risks_ij[parameters.name].apply(lambda row: (1. - exp(-DELTAT / float(mypedestrians.loc[ mypedestrians["ped_no"]==row["ped_no_A"] ]["tau"])  * row["sum_risk_AtoB"])), axis=1)
	
	print(Risks_ij[parameters.name].loc[ Risks_ij[parameters.name]["exp_risk_AtoB"]>1e-10 ])
	
#####################################################
	
	

	
with open(SaveRisksFile,'w') as monfichier:

	monfichier.write("Parameters;Clow;Cbar")
	
	for parameters in model_parameters.values():
		## UPPER BOUND on the number of infected individuals by row.ped_no
		myrisks= Risks_ij[parameters.name]
		mypedestrians["Cbari_%s"%parameters.name] = mypedestrians.apply(lambda row: (myrisks.loc[ (myrisks["ped_no_A"]==row["ped_no"]) ])["exp_risk_AtoB"].sum(),	axis=1)
			

		## LOWER BOUND
		mypedestrians["Clowi_%s"%parameters.name] = mypedestrians["Cbari_%s"%parameters.name]
		# the third parameter is the upper bound on nb of individuals who were infected during tau_together
		
		Clow= mypedestrians["Clowi_%s"%parameters.name].mean()
		Cbar= mypedestrians["Cbari_%s"%parameters.name].mean()
		print(" Mean number of individuals infected by a diseased person circulating in the crowd for %i minutes:\n                %.4e - %.4e "
						%(DELTAT/60.,
						Clow,
						Cbar))
		
		
		monfichier.write("\n%s;%.4e;%.4e"%(parameters.name,Clow,Cbar))


####### DUMP DATA #######
mypedestrians.to_pickle(DumpFileRisks,protocol=3)



