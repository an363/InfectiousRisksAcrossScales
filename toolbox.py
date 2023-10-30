#! coding: utf-8 
import numpy as np, time
from math import exp, pi, floor, cos, sin, sqrt, atan
import pandas as pd
from random import random
import os, re, sys, time
import pickle


### BUILT-IN PARAMETERS ###
delta_tau= 0.2 # bins of 0.2 seconds for the concentration maps
delta_r= 0.2 # bins of radius 0.3 m
delta_theta= np.pi/12.0
								  
speed_walks= np.array( list(map(lambda x : round(x,1), np.arange(0,2.0,0.1))))
speeds= np.array([0.0,0.3,1.0,2.0])
speed_max= max(speeds)

phis_degrees= np.array([0, 30., 60., 90., 135., 180.]) # between - pi and pi

DELTAT= 3600. # time window (in seconds) over which infection risks are computed
tau_max= 21.0 # maximum delay between droplet emission and inhalation -- leave as it is, unless you use new dynamic concentration maps
### ###
######## SCENARIO CHARACTERISTICS #######
# The following parameters help compensate an inexhaustive sampling of pedestrians. 
# If the field of view was small and some people thus failed to be detected, we will try to compensate for missed interactions with people off camera by renormalising the interaction frequency
# and to account for past interactions between pedestrians i and j (who might already have infected one another on the premises, before they were observed by us.
# Set these parameters to very large values if all pedestrians could be observed / simulated.
(Lx,Ly)= (200,200) # x- and y- dimensions (in metres) of the field of view in which pedestrian -- if the field of view was small,
	
tau_together= 120. # max time (in seconds) spent together by pedestrians i and j -
	






headers= ["TrajectoryFolder","DiagramsFolder","OutputFolder","T0","(vx,vy)","ExhalationMode","IsotropicInhalation","ContagionAmidGroups"]
spatiotemp_diagrams_all= {}
Z_all= {}
######### MODEL PARAMETERS (read from input file) ######### 

def isStatic(filename):
	if filename[-7:]=='_static':
			return True
	elif filename[-8:]=='_dynamic':
			return False
	else:
		print("Please specify for each scenario if it is static (folder name should end with '_static' suffix) or not ('_dynamic'), notably for %s"%filename)
		quit()
		
def read_params(InputFile):

	myparams= {} # dictionary where the parameters will be stored
	
	with open(InputFile,'r') as monfichier:
		for maligne in monfichier.readlines():
			ligne= maligne.strip("\n").split("=")
			ligne= [ mot.strip(" ") for mot in ligne if mot.strip(" ")!='' ]
			
			header= ligne[0]
			
			if header not in headers:
				print("Unknown header of input file: %s "%header)
				quit()
			if len(ligne)<=1:
				print("Please enter value for header %s of the input file"%header)
				quit()
				
			# read either one value or a list of values separated by semi-colons
			myparams[header]= list([])
			ligne= ligne[1].split("#")[0].strip(" ")
			for param in ligne.split(";"):
				value= param.strip(" ")
				
				if header in ["IsotropicInhalation","ContagionAmidGroups"]:
					if value not in ["True", "False"]:
						print("Parameters 'IsotropicInhalation' and 'ContagionAmidGroups' can only take values True or False.")
						quit()
					else:
						value= True if value=='True' else False
				elif header=="(vx,vy)": 
					try:
						m= re.match("\((-?\d+|-?\d+.\d+),(-?\d+|-?\d+.\d+)\)", value)
						value= [ float(m[1]), float(m[2]) ]
					except:
						print("Wind speed %s is not written adequately; please express it as (vx,vy)"%value)
						quit()
				elif header=="T0":
					try:
						value= float(value)
					except:
						print("Characteristic infection time T0 value '%s' is not a float"%value)
						quit()
				myparams[header].append(value)
	
	for header in headers:
		if header not in myparams.keys():
			print("Missing header in input file: %s"%header)
			quit()
			
	for emode in myparams["ExhalationMode"]:
		if emode not in ["breathing","speaking", "large_droplets"]:
			print("Unknown exhalation mode: %s.\n Please select among: 'breathing','speaking', 'large_droplets'"%emode)
			quit()
	
	# Check whether the scenarios are static or dynamic
	for trajFolder in myparams["TrajectoryFolder"]:
		isStatic(trajFolder)
		
	# count number of parameter sets
	nb_param_sets= 1
	for header in headers:
		#print("%s has %i values"%(header,len(myparams[header])))
		if (header in ["DiagramsFolder","OutputFolder"]) and len(myparams[header])>1:
			print("Please only give one value for parameter %s in input file."%header)
			quit()
		nb_param_sets*= len(myparams[header])
		
	return (nb_param_sets, myparams)
	

	

def getDetections(FolderTraj):
	
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
	detections["group_no"]= detections["group_no"].astype(int)
	detections["ped_no"]= detections["ped_no"].astype(int)
	detections.sort_values(by=["ped_no","time"], inplace= True)
	detections.reset_index(inplace=True)
	del detections_dict


	
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
			

	times= sorted(detections["time"].unique())
	dt= min([times[cpt+1] - times[cpt] for cpt in range(len(times)-1)]) # a greedy way to find the elementary time step
	dt= round(dt,3)
	print(" The time step is dt = ", dt)




	####### PEDESTRIANS & GROUPS #######


	# PEDESTRIANS #
	group_from_ped= {} # dictionary such that group_from_ped[ myped ] = mygroup
	for row in detections.itertuples():
		group_from_ped[ row.ped_no ] = int(row.group_no)

	mypedestrians= pd.DataFrame(data={'ped_no': np.array(sorted(group_from_ped.keys()))} )
	print("Total number of pedestrians :",mypedestrians.shape[0])

	mypedestrians["ped_no"]= mypedestrians["ped_no"].astype(int)

	mypedestrians["tau"]= mypedestrians["ped_no"].apply( lambda ped_no: 0.5 + dt * (detections.loc[ detections["ped_no"]==ped_no ].shape[0] )  )  # time spent in the field of view, with an additional 0.5 s due to mis-reconstruction of beginning and end of trajectory
	mypedestrians["group_no"]= mypedestrians["ped_no"].apply(lambda x: int(group_from_ped[x]) )


	detections.sort_values(by=["ped_no","time"], inplace=True)

	scenario= {"dt": dt, "detections": detections, "group_from_ped": group_from_ped, "mypedestrians": mypedestrians}
	return scenario


######### DEFINITION OF THE TRANSMISSION MODEL #########


def get_norm(v):
	return sqrt( v[0]**2+v[1]**2 )


def make_filename(speed,speed_walk,phi,param_set):
	if speed==0.0:
		return param_set["DiagramsFolder"]+"/diagram_v=0.0_vwalk=%.1f_phi=0_ICS_%s.pkl"%(speed_walk,param_set["ExhalationMode"]) 
	return param_set["DiagramsFolder"]+"/diagram_v=%.1f_vwalk=%.1f_phi=%d_ICS_%s.pkl"%(speed,
														  speed_walk,
														  phi,
                                            param_set["ExhalationMode"])
	
def convert_theta(x, dx):
    return round(dx * round(x/dx,0),2)

def convert_others(x,dx):
    return round(dx * np.floor(x/dx+1e-3),1)


phis_rad= np.array(list(map(lambda phi: convert_theta(pi/180. * phi,delta_theta), phis_degrees)))
phi_max= max(phis_rad)
phi_min= min(phis_rad)

def load_diagrams(param_set):
	# Loads spatiotemporal diagrams of risks from Pickle files and stores them in the following structure
	#  spatiotemp_diagrams[(speed,speed_walk,phi)][delay_cg][(r,theta)]= value
	spatiotemp_diagrams= {} 

	for speed_walk in speed_walks:
		for speed in speeds:
			for cpt_phi in range(len(phis_degrees)):
				st_key=  (speed, speed_walk, phis_rad[cpt_phi])
				spatiotemp_diagrams[ st_key ]= {}
				
				diagram_aux=  pickle.load( open(make_filename(speed,speed_walk,phis_degrees[cpt_phi],param_set),'rb') )
				for key,value in diagram_aux.items():
					delay_cg= key[0] 
					if delay_cg>tau_max or value<1E-20: # Get rid of useless cases
						continue					
					if delay_cg not in spatiotemp_diagrams[st_key].keys():
						 spatiotemp_diagrams[st_key][delay_cg]= {}
					r_cg= key[1]
					if r_cg not in  spatiotemp_diagrams[st_key][delay_cg].keys():
						 spatiotemp_diagrams[st_key][delay_cg][r_cg]= {}
					spatiotemp_diagrams[st_key][delay_cg][r_cg][ key[2] ]= value / delta_tau # division by delta_tau inserted on July, 25th 2022
	
	print(" * I have finished loading the spatiotemporal diagrams of risks")
	return spatiotemp_diagrams





def get_nu_dynamic(v_walk,r_vec,theta_e,theta_r,delay,param_set): # 
	# All coordinates are expressed in a frame centered on the emitter E and of arbitrary (absolute) orientation
	# v_walk: walking velocity in lab frame
	# r_vec= r_Receiver - r_Emitter
	# theta_e: direction of exhalation in LAB FRAME
	# theta_r: receiver head's orienation in LAB FRAME
	# delay is the (non-coarse-grained) delay
	
	# Emitted-related factors
	v_wind= param_set["(vx,vy)"] # v_wind: wind velocity in lab frame
	speed_walk= get_norm(v_walk)
	speed_wind= get_norm(v_wind)
	v_vec= (v_wind[0]-v_walk[0], v_wind[1]-v_walk[1])
	speed= get_norm(v_vec)
	
	
	phi= np.angle(v_walk[0] + 1.0j*v_walk[1]) - np.angle(-v_vec[0] - 1.0j*v_vec[1]) # angle between walking direction and MINUS v
	delta= theta_e - np.angle( v_walk[0] + 1.0j*v_walk[1] ) # angle between head direction and walking direction 
	phiprime= phi + delta  # angle between direction of emission and MINUS v
	phiprime= (phiprime+pi)%(2.0*pi)-pi
	
	r_cg= convert_others( get_norm(r_vec), delta_r)
	theta_ER= np.angle(r_vec[0] + 1.0j*r_vec[1]) - theta_e # angle between receiver and direction of exhalation
	theta_ER= (theta_ER+pi)%(2.0*pi)-pi
	theta_ER_cg= convert_theta(theta_ER, delta_theta)
	
	
	if phiprime<0.: # correct for previous symmetrisation
		theta_ER=-theta_ER
	
	# Case of speed>2m/s is handled similarly to any interpolation with 0<speed<2, where the delay tau is adjusted (via tau*speed= tau2*speed2) and risks are reweighted by speed/speed2
	diagram= interpolate(speed_walk, speed, abs(phiprime), delay, r_cg, param_set["ExhalationMode"]) # symmetrised version for phiprime<0
	
	
	
		
	# Receiver-related factors
	theta_R= theta_r - np.angle(-r_vec[0] - 1.0j*r_vec[1])
	theta_R= (theta_R+pi)%(2.0*pi)-pi

	if param_set["IsotropicInhalation"]:
		receiver_factor= 1.0 # isotropic case
	else:
		receiver_factor= 1.0 if abs(theta_R)<0.5*pi else 0.0 # case in which inhalation vanishes to zero when face is opposed to main direction

		
	return diagram.get(theta_ER_cg, 0.0) * receiver_factor / param_set["T0"] / Z_all[ param_set["ExhalationMode"] ]
#######################################################

############ INTERPOLATION #################

def rotate(diagram, angle):
		
	diagram_rot= {} # diagram_rot[theta]
	for theta_cg, value in diagram.items(): #!! Useless computations
		theta_r= convert_theta(theta_cg-angle, delta_theta)
		diagram_rot[theta_r]= value
	return diagram_rot

def pondere(diagramA, diagramB, coefA, coefB=-1): # returns coefA*diagramA + coefB*diagramB, where coefB=1-coefA, unless explicitly specificied
	if coefB==-1:
		coefB= 1.0-coefA
		
	diagram= {}
	used_keys= set([])
	for key, value in diagramA.items():
		diagram[key]= coefA * value + coefB * diagramB.get(key,0.0)
		used_keys.add(key)
	for key, value in diagramB.items():
		if key not in used_keys:
			diagram[key]= coefA *  diagramA.get(key,0.0) + coefB * value
			used_keys.add(key)
	return diagram
	
		
def interpolate(speed_walk, speed, angle_phi, delay, r_cg, ExhalationMode):
	## Two interpolations are performed: 
	# 1) On the angle, rotate and reweight. 
	# 2) On the speeds, accelerate movie (change tau and add prefactor to compensate for acceleration) and reweight.
	# only consider the diagram shell at r=r_cg
	
	# If delay>tau_crit, rotate reference diagrams to reach the right phi angle imposed by the wind; otherwise, no rotation
	# tau_crit is obtained by considering that emission speed has become negligible after 0.4s at v=0.3m/s
	tau_crit= 0.4 * 0.3 / (speed+0.00001) 
	
	
	# angle_phi should always be between 0 and pi
	speed_walk_rd= round(speed_walk,1)
	if speed_walk_rd> max(speed_walks):
		#print(">>>Requested speed_walk is larger than max available")
		speed_walk_rd= max(speed_walks)
	
	speed1= max( np.where(speeds<=speed,speeds,-1) )
	speed2= min( np.where(speeds>speed,speeds, speed_max) )
	
	
	time_accel1= speed/speed1 if speed1>0.0 else 1.0
	delay1= convert_others(delay * time_accel1, delta_tau)
	time_accel2= speed/speed2 if speed>0.0 else 1.0
	delay2= convert_others(delay * time_accel2, delta_tau)	
	
	diff_phi= np.array([ (phi - angle_phi) for phi in phis_rad])
	
	diff_phiA= max( np.where(diff_phi<=0.0,diff_phi,phi_min-angle_phi) )
	diff_phiB= min( np.where(diff_phi>=0.0,diff_phi,phi_max-angle_phi) ) # angle_phiA and angle_phiB may be equal, if both are equal to angle_phi
	
	angle_phiA= diff_phiA+angle_phi
	angle_phiB= diff_phiB+angle_phi

	# print("Interpolating diagram at ", speed_walk, speed, angle_phi, delay, " with ", speed1,delay1,angle_phiA, " and ", speed2, delay2,angle_phiB)

	spatiotemp_diagrams= spatiotemp_diagrams_all[ ExhalationMode ]	
	try:
		diagram1A= spatiotemp_diagrams.get( (speed1,speed_walk_rd,angle_phiA), -1).get(delay1,{}).get(r_cg,{}) # empty dictionary if key is not found
		diagram1B= spatiotemp_diagrams.get( (speed1,speed_walk_rd,angle_phiB), -1).get(delay1,{}).get(r_cg,{})
	except:
		print("Could not find diagrams1 at: ", speed1,speed_walk_rd,angle_phiA,angle_phiB, r_cg)
		quit()
	#print(speed1,speed_walk_rd,angle_phiA)
	#print(diagram1A, diagram1B)
	
	if delay>tau_crit and speed1>0.0:
		diagram1A= rotate(diagram1A, diff_phiA)
		diagram1B= rotate(diagram1B, diff_phiB)
	coefA= abs(diff_phiB) / ( abs(diff_phiA) + abs(diff_phiB) ) if ( abs(diff_phiA) + abs(diff_phiB) )>0 else 1.0
	
	diagram1= pondere(diagram1A, diagram1B, coefA) 
	
	diagram2A= spatiotemp_diagrams.get( (speed2,speed_walk_rd,angle_phiA), -1).get(delay2,{}).get(r_cg,{})
	diagram2B= spatiotemp_diagrams.get( (speed2,speed_walk_rd,angle_phiB), -1).get(delay2,{}).get(r_cg,{})
	
	
	if delay>tau_crit:
		diagram2A= rotate(diagram2A, diff_phiA)
		diagram2B= rotate(diagram2B, diff_phiB)
	diagram2= pondere(diagram2A, diagram2B, coefA)
	
	#if diagram1A=={} or diagram1B=={} or diagram2A=={} or diagram2B== {}:
	#	print("Could not find associated diagram for delay ", delay1, delay2, " at v=", speed)
	
	#print(angle_phiA, angle_phiB, speed1, speed2 )
	
	coef1= (speed2-speed) / (speed2-speed1) if (speed2>speed1) else 0.0 # condition if (len(diagram1)>0 and speed2>speed1) trimmed on July, 4th
	
	del diagram1A, diagram1B, diagram2A, diagram2B
	return pondere(diagram1, diagram2, coef1/time_accel1, (1.0-coef1)/time_accel2 )

# if speed larger than v2=2m/s, just reason in terms of v1*tau1 = v2*tau2


###########################################


def arrondit(x, precision):
	return round(x/precision)*precision
	
	
	
####### Correct for EDGE EFFECTS #######
# Problem: With a finite-size field of view, some contacts are missed
# Solution: 
# Points are randomly drawn in the field of view, and, for each one, we determine if its potential neighbour at distance r in 
# a random direction is also within the rectangles or out of bounds.
# This defines an R-dependent rescaling factor, with which all contact durations shall be multiplied for rescaling purposes.


def get_rescaling_factor(Lx,Ly):
	Rvec= np.arange(0,30,0.1)
	rescaling_factor= {}
	nb_points= 10**4 # number of random trial points
	
	for cpt in range(len(Rvec)):
		R= Rvec[cpt]
		nb_points_out_of_rectangle= 0
		for k in range(nb_points):
			x_c= random() * Lx
			y_c= random() * Ly
			theta= random() * 2.0 * pi
			x,y= (x_c + R * cos(theta),y_c+ R * sin(theta))
			if x<0 or x>Lx or y<0 or y>Ly: # out of bounds
				nb_points_out_of_rectangle+=1

				
		rescaling_factor[int(10*R)]= max(1.0, 1.0 / (1.0 + 0.05 - nb_points_out_of_rectangle / nb_points)) # 0.05 is to avoid diverging rescaling factors
	return rescaling_factor
	
	
################################ COMPUTE TRANSMISSION RISKS BETWEEN INDIVIDUAL PEDESTRIANS ##################################################
# A) DYNAMIC CASE

def getDynamicRisks(scenario, parameters):

	#return pd.DataFrame([],columns=["ped_no_A","ped_no_B","exp_risk_AtoB"])
	dt= scenario["dt"]
	detections= scenario["detections"]
	group_from_ped= scenario["group_from_ped"]
	mypedestrians= scenario["mypedestrians"]
	
	####### Correct for EDGE EFFECTS #######
	rescaling_factor= get_rescaling_factor(Lx,Ly)


	
	# Prepare empty risks: risks_AB is the risk incurred by B (receiver) because of A (emitter)
	risks_AB= {}
	for ped_A in group_from_ped.keys():
		for ped_B in group_from_ped.keys():
			risks_AB[(ped_A,ped_B)]= 0.0
	
		
	last_hped= -1
	for row_A in detections.itertuples():
	
		if row_A.ped_no!=last_hped:
			last_hped= row_A.ped_no
			# Display a progress bar
			progress_pc= 100*row_A.Index//detections.shape[0]
			progress_string= "Simulation progress   |"
			for cpt in range(20):
				progress_string+=">" if cpt<=progress_pc/5 else " "
			progress_string+= "| (%i%%)"%(progress_pc)
			print(progress_string, end='\r')
			#print(" Currently handling ped n.%i"%row_A.ped_no)

		for row_B in detections.loc[ (detections["ped_no"]!=row_A.ped_no) & (detections["time"]>=row_A.time) & (detections["time"]<row_A.time+tau_max) ].itertuples():
			r_vec= (row_B.x-row_A.x, row_B.y-row_A.y)
			r= get_norm(r_vec)
				
			#!! For computational efficiency, discard people more than 3 metres away
			if r>4.0:
				continue
				
			tau= row_B.time - row_A.time
			v_walk= (row_A.vx, row_A.vy)
				
			risks_AB[(row_A.ped_no,row_B.ped_no)]+= dt * dt * rescaling_factor.get( int(round(10*r)), 1) * get_nu_dynamic(v_walk, r_vec, row_A.theta, row_B.theta, tau, parameters)
			# the rescaling coefficient accounts, at the mean-field level, for the edge effects whereby distant contacts are not seen
			# dt transformed in dt**2 on July 25, 2022
				
				
	# Drop entries with zero risks
	for key in list(risks_AB.keys()):
		if risks_AB[key]<1e-18:
			del risks_AB[key]
		
		
	Risks_ij= pd.DataFrame.from_dict(risks_AB, orient='index', columns=["sum_risk_AtoB"])
	Risks_ij.reset_index(inplace=True,drop=False)
	Risks_ij["ped_no_A"]= Risks_ij["index"].apply(lambda t: t[0])
	Risks_ij["ped_no_B"]= Risks_ij["index"].apply(lambda t: t[1])
	Risks_ij.drop(columns=['index'], inplace=True)
		
		
	Risks_ij["group_no_A"]= Risks_ij["ped_no_A"].apply(lambda x: int(group_from_ped[x]))
	Risks_ij["group_no_B"]= Risks_ij["ped_no_B"].apply(lambda x: int(group_from_ped[x]))
		
	if parameters["ContagionAmidGroups"]:
		Risks_ij["exp_risk_AtoB"]= Risks_ij["sum_risk_AtoB"].apply(lambda x: 1. - exp(-x) ) # this is the mean infection state I_B of B after contacts with diseased A
	else: 
		# a diseased individual cannot infect members of his/her own group
		Risks_ij["exp_risk_AtoB"]= Risks_ij.apply(lambda row: (1. - exp(-row["sum_risk_AtoB"])) if row["group_no_A"]!=row["group_no_B"] else 0.0, axis=1)
			
		
	return (Risks_ij.loc[ Risks_ij["exp_risk_AtoB"]>1e-10 ])

def get_lower_bound_sum(row,myrisks,Cbar):
	# this returns a lower bound on the number of agents infected by agent i (row["ped_no"]=i), under the constraint that agent i has infected at most Cbar individual during "tau_together"
	# get rid / don't get rid of those in the same group
	
	risks_loc= list(myrisks.loc[ (myrisks["ped_no_A"]==row["ped_no"]) ]["exp_risk_AtoB"].values) 
	# first term: if A=row.ped_no is the infected agent; second term: if B=row.ped_no is the infected agent
	
	# sort risks by decreasing order
	risks_loc.sort(reverse=True)
				
	# if agent i has had very few contacts
	if ( len(risks_loc)-1 < int(floor(Cbar)) ):
		return 0.0
		
	# otherwise, remove closest contacts by changing probability associated with next-to... closest agent to (1 - Cbari%1)
	C_floor= int(floor(Cbar))
	somme=  (1.0 - (Cbar%1.0)) * risks_loc[C_floor]
	
	# remove closest-agent terms of the sum because these were potentially already infected
	somme+= np.sum( risks_loc[C_floor+1:] )
	return DELTAT / row["tau"] * somme
	
	
def writeDynamicResults(folder,SaveRisksFile, mypedestrians, myrisks):

	with open(folder+"/Risks_mean_"+SaveRisksFile,'w') as monfichier:

		monfichier.write("Clow;Cbar")
	
		## UPPER BOUND on the number of infected individuals by row.ped_no
		mypedestrians["Cbari"] = mypedestrians.apply(lambda row: DELTAT / row["tau"] * (myrisks.loc[ (myrisks["ped_no_A"]==row["ped_no"]) ])["exp_risk_AtoB"].sum(),	axis=1)
			

		## LOWER BOUND
		mypedestrians["Clowi"] = mypedestrians.apply(lambda row: get_lower_bound_sum(row,
																										myrisks, 
																										tau_together / DELTAT * row["Cbari"]), axis=1)
		# the third parameter is the upper bound on nb of individuals who were infected during tau_together
		
		Clow= mypedestrians["Clowi"].mean()
		Cbar= mypedestrians["Cbari"].mean()
		print(" Mean number of individuals infected by a diseased person circulating in the crowd for %i minutes:\n                %.4e - %.4e "
						%(DELTAT/60.,
						Clow,
						Cbar))
		
		
		monfichier.write("\n%.4e;%.4e"%(Clow,Cbar))

	
	with open(folder+"/Risks_by_person_"+SaveRisksFile,'w') as monfichier:
		monfichier.write("ped_no;Clow;Cbar")
		for row in mypedestrians.itertuples():
			monfichier.write("\n%i;%.4e;%.4e"%(row.ped_no,row.Clowi,row.Cbari))
			
	

# B) STATIC CASE


def getStaticRisks(scenario, parameters):
	
	dt= scenario["dt"]
	detections= scenario["detections"]
	group_from_ped= scenario["group_from_ped"]
	mypedestrians= scenario["mypedestrians"]
	
	####### Correct for EDGE EFFECTS #######
	rescaling_factor= get_rescaling_factor(Lx,Ly)



	# Prepare empty risks: risks_AB is the risk incurred by B (receiver) because of A (emitter)
	risks_AB= {}
	for ped_A in group_from_ped.keys():
		for ped_B in group_from_ped.keys():
			risks_AB[(ped_A,ped_B)]= 0.0
	
	last_hped= -1
	for row_A in detections.itertuples():
		if row_A.ped_no!=last_hped:
			last_hped= row_A.ped_no
			# Display a progress bar
			progress_pc= 100*row_A.Index//detections.shape[0]
			progress_string= "Simulation progress   |"
			for cpt in range(20):
				progress_string+=">" if cpt<=progress_pc/5 else " "
			progress_string+= "| (%i%%)"%(progress_pc)
			print(progress_string, end='\r')
			#print(" Currently handling ped n.%i"%row_A.ped_no)
			
		detections_loc= detections.loc[ (detections["ped_no"]!=row_A.ped_no) & (detections["time"]>=row_A.time) & (detections["time"]<row_A.time+tau_max) ]

		for row_B in detections_loc.itertuples():
			r_vec= (row_B.x-row_A.x, row_B.y-row_A.y)
			r= get_norm(r_vec)
			
			tau= row_B.time - row_A.time

			#!! For computational efficiency, discard people more than 4 metres away
			if r>4.0:
				continue
			v_walk= (row_A.vx, row_A.vy)
			
			
			risks_AB[(row_A.ped_no,row_B.ped_no)]+= dt * dt * rescaling_factor.get( int(round(10*r)), 1) * get_nu_dynamic(v_walk, r_vec, row_A.theta, row_B.theta, tau, parameters)
			# the rescaling coefficient accounts, at the mean-field level, for the edge effects whereby distant contacts are not seen
			# dt transformed in dt**2 on July 25, 2022
			
				
			
			
	# Drop entries with zero risks
	for key in list(risks_AB.keys()):
		if risks_AB[key]<1e-18:
			del risks_AB[key]
	
	Risks_ij= pd.DataFrame.from_dict(risks_AB, orient='index', columns=["sum_risk_AtoB"])
	Risks_ij.reset_index(inplace=True,drop=False)
	Risks_ij["ped_no_A"]= Risks_ij["index"].apply(lambda t: int(t[0]))
	Risks_ij["ped_no_B"]= Risks_ij["index"].apply(lambda t: int(t[1]))
	Risks_ij.drop(columns=['index'], inplace=True)

	
	
	Risks_ij["group_no_A"]= Risks_ij["ped_no_A"].apply(lambda x: int(group_from_ped[x]))
	Risks_ij["group_no_B"]= Risks_ij["ped_no_B"].apply(lambda x: int(group_from_ped[x]))
	
	# a diseased individual MAY or MAYinfect members of his/her own group in these static scenarios
	if parameters["ContagionAmidGroups"]:
		Risks_ij["exp_risk_AtoB"]= Risks_ij.apply(lambda row: (1. - exp(-DELTAT / float(mypedestrians.loc[ mypedestrians["ped_no"]==row["ped_no_A"] ]["tau"])  * row["sum_risk_AtoB"])), axis=1)
	else: 
		# a diseased individual cannot infect members of his/her own group
		Risks_ij["exp_risk_AtoB"]= Risks_ij.apply(lambda row: (1. - exp(-DELTAT / float(mypedestrians.loc[ mypedestrians["ped_no"]==row["ped_no_A"] ]["tau"])  * row["sum_risk_AtoB"])) if row["group_no_A"]!=row["group_no_B"] else 0.0, axis=1)
		
	
	return (Risks_ij.loc[ Risks_ij["exp_risk_AtoB"]>1e-10 ])
	
#####################################################
	
	

def writeStaticResults(folder,SaveRisksFile, mypedestrians, myrisks):

	with open(folder+"/Risks_mean_"+SaveRisksFile,'w') as monfichier:

		monfichier.write("Clow;Cbar")
	
		## UPPER BOUND on the number of infected individuals by row.ped_no
		mypedestrians["Cbari"] = mypedestrians.apply(lambda row: (myrisks.loc[ (myrisks["ped_no_A"]==row["ped_no"]) ])["exp_risk_AtoB"].sum(),	axis=1)
			

		## LOWER BOUND
		mypedestrians["Clowi"] = mypedestrians["Cbari"]
		# the third parameter is the upper bound on nb of individuals who were infected during tau_together
		
		Clow= mypedestrians["Clowi"].mean()
		Cbar= mypedestrians["Cbari"].mean()
		print(" Mean number of individuals infected by a diseased person circulating in the crowd for %i minutes:\n                %.4e - %.4e "
						%(DELTAT/60.,
						Clow,
						Cbar))
		
		monfichier.write("\n%.4e;%.4e"%(Clow,Cbar))

	with open(folder+"/Risks_by_person_"+SaveRisksFile,'w') as monfichier:
		monfichier.write("ped_no;Clow;Cbar")
		for row in mypedestrians.itertuples():
			monfichier.write("\n%i;%.4e;%.4e"%(row.ped_no,row.Clowi,row.Cbari))

