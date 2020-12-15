
# Parameters and Properties Definitions

# Simulation Parameters 
import numpy as np
#   Starting and Final Times for Simulation
#   =======================================
ti      = 	0.0      # Initial time for simulation, sec
tf      =	60.0     # Final time for simulation, sec
h = .1               # Define Step Size
t = np.arange(ti,tf+h,h)


tuHis	=np.array([0,5,6,7,60])
deluHis = np.zeros((5,7))

deluHis[0,:] = [0*3*np.pi/180, 0, 0, 0, 0, 0, 0]
deluHis[1,:] = [0*3*np.pi/180, 0, 0, 0, 0, 0, 0]
deluHis[2,:] = [0, 0, 0, 0, 0, 0, 0]
deluHis[3,:] = [0, 0, 0, 0, 0, 0, 0]
deluHis[4,:] = [0, 0, 0, 0, 0, 0, 0]

#	DEFINITION OF THE STATE VECTOR
#   With Euler Angle 
#		x[0]    = 		Body-axis x inertial velocity, ub, m/s
#		x[1]    =		Body-axis y inertial velocity, vb, m/s
#		x[2]    =		Body-axis z inertial velocity, wb, m/s
#		x[3]    =		North position of center of mass WRT Earth, xe, m
#		x[4]    =		East position of center of mass WRT Earth, ye, m
#		x[5]    =		Negative of c.m. altitude WRT Earth, ze = -h, m
#		x[6]    =		Body-axis roll rate, pr, rad/s
#		x[7]    =		Body-axis pitch rate, qr, rad/s
#		x[8]    =		Body-axis yaw rate, rr,rad/s
#		x[9]   =		Roll angle of body WRT Earth, phir, rad
#		x[10]   =		Pitch angle of body WRT Earth, thetar, rad
#		x[11]   =		Yaw angle of body WRT Earth, psir, rad

#	DEFINITION OF THE CONTROL VECTOR
#		u[0]    = 		Elevator, dEr, rad, positive: trailing edge down
#		u[1]    = 		Aileron, dAr, rad, positive: left trailing edge down
#		u[2]    = 		Rudder, dRr, rad, positive: trailing edge left
#		u[3]    = 		Throttle, dT, %
#		u[4]    =		Asymmetric Spoiler, dASr, rad
#		u[5]    =		Flap, dFr, rad
#		u[6]    =		Stabilator, dSr, rad



#	Typical Mass and Inertial Properties
Mass	=	4800 				# Mass, kg
Ixx		=	20950    			# Roll Moment of Inertia, kg-m **2
Iyy		=	49675     			# Pitch Moment of Inertia, kg-m **2
Izz		=	62525     			# Yaw Moment of Inertia, kg-m **2
Ixz		=	-1710    			# Nose-high(low) Product of Inertia, kg-m **2

#	Geometric Properties
cBar	=	3.03 				# Mean Aerodynamic Chord, m
b		=	10 				    # Wing Span, m
S		=	27.77 				# Reference Area, m **2
#ARw		=	5.02 				# Wing Aspect Ratio
#taperw	=	0.507 				# Wing Taper Ratio
#sweepw	=	13 * .01745329 		# Wing 1/4-chord sweep angle, rad
#ARh		=	4 					# Horizontal Tail Aspect Ratio
#sweeph	=	25 * .01745329 		# Horizontal Tail 1/4-chord sweep angle, rad
#ARv		=	0.64 				# Vertical Tail Aspect Ratio
#sweepv	=	40 * .01745329 		# Vert Tail 1/4-chord sweep angle, rad
#lvt		=	4.72 				# Vert Tail Length, m


#	Aerodynamic Properties


#cBar            =   3.03 
#b               =   10 
#S               =   27.77 
lHT             =   5.2 
lVT             =   3.9 
StaticThrust    =   49000 

    
#	Current Longitudinal Characteristics
#	====================================
#	Lift Coefficient
CLo     =   0 
CLa     =   4.92 
CLqhat  =   2.49 
CLdE    =   0.72 
CLdS    =   CLdE 
    
#   Drag Coefficient
CDo     =   0.019 
    
#	Pitching Moment Coefficient
StaticMargin    =   0.2 
Cmo     =   0 
Cma 	=	-CLa*StaticMargin 
Cmqhat  =   -4.3 
CmV     =   0 
CmdE	=	-1.25 
CmdS 	=	CmdE 
    
#	Current Lateral-Directional Characteristics
#	===========================================
#	Side-Force Coefficient
CYo     =   0 
CYb     =   -0.5 
CYphat  =   0 
CYrhat  =   0
CYdA	=	0 
CYdR	=	0.04 
    
#       Total Side-Force Coefficient, w/Mach Correction
CYo     =   0 
    
#	Rolling Moment Coefficient
Clo     =   0 
Clb  	=	-0.066 
Clphat  =   -0.5 
Clrhat  =   -0.5 
CldA 	=	0.12 
CldR 	=	0.03 
    
   
#	Yawing Moment Coefficient
Cno     =   0 
CnBeta  =	0.37 
Cnphat  =   -0.06 
Cnrhat  =   -0.5 
CndA 	=	0 
CndR 	=	-0.2 

# Input Gust Parameters 
# "1 - Cosine" gust input 
gust_amp_1_minus_cos = 0       # max velocity of "1 - cosine" gust   (m/s)
gust_t = .0                   #  fraction of total data time length that is gust    0 - 1 

# TURBULENCE INPUT - random signal between 0 Hz and turb_max_freq Hz
    #   uniform amplitude across random frequency input
turb_amp = 0       # max vertical velocity of turbulence   (m/s)
turb_t = 0.0       #  fraction of total time length that is turbulence    0 - 1 
turb_max_freq = 100  # max frequency of the turbulence (Hz) 

