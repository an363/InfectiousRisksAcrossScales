#! coding: utf-8 
from toolbox import *

##########################
# This script takes as input the contact tables created with 
# the trajectories extracted from the videos and 
# assesses the transmission risks with different models.
# It prints out the results and saves them in a dump file,
# to be processed by risk_analysis.py
# W. Garcia, A. Nicolas (2020), A. Nicolas (2022)
##########################



print("Pedestrians can%s infect members of their own group"%("" if ContagionAmidGroups else "not"))
DumpFileRisks= PickleFolder + "/Risks" + "/Risks" + suffix + ".pkl" 
DumpFileContacts= PickleFolder +"/Risks" + "/df_Contacts_drisk" + suffix + ".pkl"
SaveRisksFile= PickleFolder + "/Risks" + "/risks" + suffix + ".csv"

if not os.path.isdir(PickleFolder + "/Risks"):
	os.mkdir(PickleFolder + "/Risks")

detections= pd.read_pickle(PickleFolder+"/df.pkl") # dataframe with one-body detections







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




# theta is an angle in an arbitrary absolute frame attached to Earth

"""
print(detections.columns)
detections["v2"]= detections.apply(lambda row: sqrt(row.vx**2+row.vy**2), axis=1)

detections["th"]= detections.apply(lambda row: (row.theta-row.theta_walk+pi)%(2.0*pi)-pi, axis=1)

detections.hist("th", bins=100)
plt.show()
quit()
"""
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
######## Histogram of risks as a function of distance
"""
fig,ax=plt.subplots(5,5, figsize=(15,15))
cpt=0

contacts_less_1m= contacts.loc[ contacts["r"]< 1.0 ]
contacts_1m_2m= contacts.loc[ (contacts["r"]>= 1.0) & (contacts["r"]< 2.0) ]
contacts_more_2m= contacts.loc[ (contacts["r"]> 2.0) ]

for parameters in model_parameters.values():
	ax[cpt//5][cpt%5].set_title(parameters.name)
	ax[cpt//5][cpt%5].hist(contacts["r"].values, bins=np.arange(0,3,0.33), weights=contacts["Risk_%s"%parameters.name].values, density= True)
	
	ratio_less_1m= contacts_less_1m["Risk_%s"%parameters.name].sum() / contacts["Risk_%s"%parameters.name].sum()
	ratio_1m_2m= contacts_1m_2m["Risk_%s"%parameters.name].sum() / contacts["Risk_%s"%parameters.name].sum()
	
	ax[cpt//5][cpt%5].text(1,1,"<1m: %i p.c."%(ratio_less_1m*100.0))
	ax[cpt//5][cpt%5].text(1,0.2,"1m-2m: %i p.c."%(ratio_1m_2m*100.0))
	cpt+=1
plt.tight_layout()
plt.savefig(GeneralFolder+"transmission_risk_vs_distance_%s.png"%(Scenarios[num_scenario]))
quit()
"""
#####################################################
	
	
	
## !! Get rid of pedestrians who stayed for too short a time on the scenario?
"""
print(" Number of pedestrians before: %i"%mypedestrians.shape[0])
mypedestrians= mypedestrians.loc[ mypedestrians["tau"]>=1.0 ]
print(" Number of pedestrians after: %i"%mypedestrians.shape[0])
"""
##
	
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


###### INVESTIGATION OF SUPER-SPREADERS #######
"""
for row in (mypedestrians.loc[ mypedestrians["Cbari_scenario5"]> 20 ]).itertuples():
	ped= row.ped_no
	contacts_loc = contacts.loc[ (contacts["ped_no_A"]==ped) | (contacts["ped_no_B"]==ped) ]
	nb_contacts= contacts_loc["ped_no_A"].nunique() + contacts_loc["ped_no_B"].nunique() - 1
	time_start= (contacts_loc).index.min()
	print("%02d:%02d   -- %.1f s (%i), %i contacts, %i infections over an hour"%(time_start//60, time_start%60, row.tau, ped, nb_contacts, row.Cbari_scenario5) )
"""

####### DUMP DATA #######
mypedestrians.to_pickle(DumpFileRisks,protocol=3)



quit()
	
contacts.drop(columns=["dtheta_A","dtheta_B", "delta_x", "delta_y", "theta_A", "theta_B", "dtheta_A_round", "dtheta_B_round", "r_angle", "TIME"], inplace=True)
for parameters in model_parameters.values():
	contacts["dRisk_%s"%parameters.name]= DELTAT / dt * (contacts["dRisk_AtoB_%s"%parameters.name] + contacts["dRisk_BtoA_%s"%parameters.name])
	contacts.drop(columns=["dRisk_AtoB_%s"%parameters.name, "dRisk_BtoA_%s"%parameters.name], inplace=True) # drop columns that have become useless

aggregation={}
for col in contacts.columns:
	if col[0:5]=="dRisk":
		aggregation[col]= 'sum'
	elif col[0:4]=="flow" or col[0:4]=="filt" or col[0:4]=="stat":
		aggregation[col]= 'mean'

contacts=contacts.groupby('time').agg(aggregation)


## Divide dRisk by number of present pedestrians
# count pedestrians per frame
ped_per_frame= {}
for time in times:
	ped_per_frame[time]= len( detections.loc[ detections["time"]==time ]["ped_no"].unique() )


contacts["nb_peds"]= contacts.index.map(lambda t: float(ped_per_frame[t]))

for parameters in model_parameters.values():
	contacts["dRisk_%s"%parameters.name]= contacts["dRisk_%s"%parameters.name] / contacts["nb_peds"]
	
	
	
# Fill missing times where there is only one pedestrian on the picture
for time in times:
	if time not in contacts.index: # if time is missing in contacts
		new_line= {}
		for parameters in model_parameters.values():
			new_line["dRisk_%s"%parameters.name]= 0.0
		new_line["nb_peds"]= 1.0
		for mykeyword in ['flow_vert_down', 'flow_vert_up', 'filter_dens', 'stationnary','flow_hor_right','flow_hor_left']:
			new_line[mykeyword+"_A"]= float((detections.loc[ (detections["time"]==time) ].iloc[0])[mykeyword])
			new_line[mykeyword+"_B"]= float((detections.loc[ (detections["time"]==time) ].iloc[0])[mykeyword])
		contacts= contacts.append(pd.Series(new_line, name="%.1f"%time))		

contacts.to_pickle(DumpFileContacts,protocol=3)

