#! coding: utf-8 
import numpy as np, time
from math import exp, pi, floor, cos, sin, sqrt, atan
import pandas as pd
from random import random
import os, re, sys, time
import pickle


######### FILES #########
FolderTraj= "Example_Data/Centre_Depistage_Gerland" # folder containing the trajectories in CSV format
FolderDiagrams= "Diagrams/" # folder containg the dynamic viral concentration maps

######### MODEL PARAMETERS ######### 
Activity= "_speaking" # Mode of emission of respiratory droplets -- choose between (i) "_speaking", (ii) "_breathing", (iii) "_large_droplets" 
ContagionAmidGroups= False # allow (or not) an agent to infect members of his/her group on site
Infection_time_T0= 900 # characteristic time (in seconds) for viral infection if somebody is facing you at 50 cm distance and talking with you 
DELTAT= 3600. # time window (in seconds) over which infection risks are computed
tau_max= 21.0 # maximum delay between droplet emission and inhalation -- leave as it is, unless you use new dynamic concentration maps
	
######## SCENARIO CHARACTERISTICS #######
# The following parameters help compensate an inexhaustive sampling of pedestrians. 
# If the field of view was small and some people thus failed to be detected, we will try to compensate for missed interactions with people off camera by renormalising the interaction frequency
# and to account for past interactions between pedestrians i and j (who might already have infected one another on the premises, before they were observed by us.
# Set these parameters to very large values if all pedestrians could be observed / simulated.
(Lx,Ly)= (200,200) # x- and y- dimensions (in metres) of the field of view in which pedestrian -- if the field of view was small,
tau_together= 120. # max time (in seconds) spent together by pedestrians i and j -
	

######## COMMAND-LINE ARGUMENTS #########
if len(sys.argv)>=3:
	v_wind= (float(sys.argv[1]), float(sys.argv[2]) )
	iso_inh= bool(int(sys.argv[3]))
else:
	print("Aborting: no input scenario nor wind speed. \n Please type in: \n\t  python3 1***.py vx_wind vy_wind isotropic_inhalation\n where\n * vx_wind: wind speed along x (in m/s) \n * vy_wind: wind speed along y (in m/s)\n * isotropic_inhalation: 0 if anisotropic inhalation factor, 1 if isotropic inhalation ")
	quit()
	



suffix= "%s_%.1f_%.1f%s_%s"%("_InfAmidGroups" if ContagionAmidGroups else "",
				  v_wind[0],
				  v_wind[1],
				  Activity,
                  "iso" if iso_inh else "aniso")


### ###

class Parameters(object):	
	
	def __init__(self, Pname, Pd0, PT0, Pthetae0, Pthetar0):
		self.name= Pname # scenario name
		self.d0= Pd0 # characteristic distance for infection
		self.T0= PT0 # characteristic time for infection
		self.thetae0= Pthetae0 # angle of emitter (e), also called psi elsewhere
		self.thetar0= Pthetar0 # angle of receiver (r), also called phi elsewhere
		
######### MODEL PARAMETERS #########
model_parameters= {}
model_parameters_from_name= {}
def set_new_model_parameters(model_parameters, name, d0, T0, thetae0, thetar0):
	nb_models= len(model_parameters)
	model_parameters[nb_models]= Parameters(Pname= "%s"%name, Pd0= d0, PT0= T0, Pthetae0= thetae0, Pthetar0= thetar0,)
	model_parameters_from_name["%s"%name]= model_parameters[nb_models]

set_new_model_parameters(model_parameters, name= 'standard', d0= 0.5, T0= Infection_time_T0, thetae0= 0, thetar0= 0)


######### DEFINITION OF THE TRANSMISSION MODEL #########


def get_norm(v):
	return sqrt( v[0]**2+v[1]**2 )


def make_filename(speed,speed_walk,phi,mode="ICS",activity=Activity):
	if speed==0.0:
		return FolderDiagrams+"diagram_v=0.0_vwalk=%.1f_phi=0_%s%s.pkl"%(speed_walk,mode,activity if activity!="_breathing" else "") 
	return FolderDiagrams+"diagram_v=%.1f_vwalk=%.1f_phi=%d_%s%s.pkl"%(speed,
														  speed_walk,
														  phi,
														  mode,
                                            activity if activity!="_breathing" else "")
	
def convert_theta(x, dx):
    return round(dx * round(x/dx,0),2)

def convert_others(x,dx):
    return round(dx * np.floor(x/dx+1e-3),1)

delta_tau= 0.2 # bins of 0.2 seconds for the concentration maps
delta_r= 0.2 # bins of radius 0.3 m
delta_theta= np.pi/12.0
								  
def load_diagrams():
	# Loads spatiotemporal diagrams of risks from Pickle files and stores them in the following structure
	#  spatiotemp_diagrams[(speed,speed_walk,phi)][delay_cg][(r,theta)]= value
	spatiotemp_diagrams= {} 

	for speed_walk in speed_walks:
		for speed in speeds:
			for cpt_phi in range(len(phis_degrees)):
				st_key=  (speed, speed_walk, phis_rad[cpt_phi])
				spatiotemp_diagrams[ st_key ]= {}
				
				diagram_aux=  pickle.load( open(make_filename(speed,speed_walk,phis_degrees[cpt_phi]),'rb') )
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


### LOADING DYNAMIC CONCENTRATION MAPS ###
speed_walks= np.array( list(map(lambda x : round(x,1), np.arange(0,2.0,0.1))))
speeds= np.array([0.0,0.3,1.0,2.0])
speed_max= max(speeds)

phis_degrees= np.array([0, 30., 60., 90., 135., 180.]) # between - pi and pi
phis_rad= np.array(list(map(lambda phi: convert_theta(pi/180. * phi,delta_theta), phis_degrees)))
phi_max= max(phis_rad)
phi_min= min(phis_rad)

spatiotemp_diagrams= load_diagrams() # Loads the dynamics viral concentration maps

### Renormalisation coefficient Z ###
Z= 1.253489635




def get_nu_dynamic(v_walk,v_wind,r_vec,theta_e,theta_r,delay,parameters): # 
	# All coordinates are expressed in a frame centered on the emitter E and of arbitrary (absolute) orientation
	# v_walk: walking velocity in lab frame
	# v_wind: wind velocity in lab frame
	# r_vec= r_Receiver - r_Emitter
	# theta_e: direction of exhalation in LAB FRAME
	# theta_r: receiver head's orienation in LAB FRAME
	# delay is the (non-coarse-grained) delay
	
	# Emitted-related factors
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
	diagram= interpolate(speed_walk, speed, abs(phiprime), delay, r_cg) # symmetrised version for phiprime<0
	
	
	
		
	# Receiver-related factors
	theta_R= theta_r - np.angle(-r_vec[0] - 1.0j*r_vec[1])
	theta_R= (theta_R+pi)%(2.0*pi)-pi

	if iso_inh:
		receiver_factor= 1.0 # isotropic case
	else:
		receiver_factor= 1.0 if abs(theta_R)<0.5*pi else 0.0 # case in which inhalation vanishes to zero when face is opposed to main direction

		
	return diagram.get(theta_ER_cg, 0.0) * receiver_factor / parameters.T0 / Z
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
	
		
def interpolate(speed_walk, speed, angle_phi, delay, r_cg):
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
	
