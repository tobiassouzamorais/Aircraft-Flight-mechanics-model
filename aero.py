# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import math as m

#	Aerodynamic Coefficients of the Aircraft, Thrust Model,
#	and Geometric and Inertial Properties

#	Typical Mass and Inertial Properties
m		=	4536 				# Mass, kg
Ixx		=	35926.5 			# Roll Moment of Inertia, kg-m **2
Iyy		=	33940.7 			# Pitch Moment of Inertia, kg-m **2
Izz		=	67085.5 			# Yaw Moment of Inertia, kg-m **2
Ixz		=	3418.17 			# Nose-high(low) Product of Inertia, kg-m **2

#	Geometric Properties
cBar	=	2.14 				# Mean Aerodynamic Chord, m
b		=	10.4 				# Wing Span, m
S		=	21.5 				# Reference Area, m **2
ARw		=	5.02 				# Wing Aspect Ratio
taperw	=	0.507 				# Wing Taper Ratio
sweepw	=	13 * .01745329 		# Wing 1/4-chord sweep angle, rad
ARh		=	4 					# Horizontal Tail Aspect Ratio
sweeph	=	25 * .01745329 		# Horizontal Tail 1/4-chord sweep angle, rad
ARv		=	0.64 				# Vertical Tail Aspect Ratio
sweepv	=	40 * .01745329 		# Vert Tail 1/4-chord sweep angle, rad
lvt		=	4.72 				# Vert Tail Length, m


def AeroModelAlpha(x,u,Mach,alphar,betar,V):

    global b,cBar,SMI
    
    #	Inertial, Geometric, and Aerodynamic Properties
    #   *** GeoMassAero.m must first be run to save
    #           InerGeo.mat,
    #           DataTable.mat,
    #           RotCont.mat ***
    
    # load InerGeo.mat
    # load DataTable.mat
    # load RotCont.mat
    
    D2R =   m.pi/180 
    R2D =   180/m.pi 
    
    alphadeg	=	R2D * alphar 
    
    #	Thrust Properties
    StaticThrust	=	2*6.49*10 **3 	# Static Thrust @ Sea Level, N
    #	Current Thrust
    airDens,airPres,temp,soundSpeed = Atmos(-x[5]) 
    Thrust			=	u[3] * StaticThrust * (airDens / 1.225) **0.7* (1 - m.exp((-x[5] - 17000) / 2000)) 
    # Thrust at Altitude, N
    
    #	Current Longitudinal Characteristics
    #	====================================
    
    #	Lift Coefficient
    CLStatic    =	interp1(AlphaTable,CLTable,alphadeg) 
    # Static Lift Coefficient
    CLqr        =	CLqHat * cBar/(2*V) 
    # Pitch-Rate Effect, per rad/s
    CLdEr   	=	interp1(AlphaTable,CLdETable,alphadeg) 
    # Elevator Effect, per rad
    CLdSr       =	CLdEr 			# Stabilator Effect, per rad
    CL          =	CLStatic + CLqr*x[7] + CLdSr*u[6] + CLdEr*u[0]
    # Total Lift Coefficient
    #	Drag Coefficient
    CDStatic	=	interp1(AlphaTable,CDTable,alphadeg) 
    # Static Drag Coefficient
    CD          =	CDStatic 		# Total Drag Coefficient
    
    #	Pitching Moment Coefficient
    CmStatic	=	interp1(AlphaTable,CmTable,alphadeg) 
    # Static Pitching Moment Coefficient
    CmdEr		=	interp1(AlphaTable,CmdETable,alphadeg) 
    # Elevator Effect, per rad
    Cmqr        =	-CLqHat*(lHT/cBar) * cBar/(2*V) 
    # Pitch-Rate + Alpha-Rate Effect, per rad/s
    CmdSr       =	CmdEr           # Stabilator Effect, per rad
    Cm          =	CmStatic - CL*SMI + Cmqr*x[7] + CmdSr*u[6] + CmdEr*u[0] 
    # Total Pitching Moment Coefficient
    
    #	Current Lateral-Directional Characteristics
    #	===========================================
    
    #	Rolling Moment Coefficient
    ClBr	=	interp1(AlphaTable,ClBetaTable,alphadeg) 
    # Dihedral Effect, per rad
    ClpHat	=	interp1(AlphaTable,ClpHatTable,alphadeg) 
    Clpr    =   ClpHat * (b / (2 * V)) 
    ClrHat	=	interp1(AlphaTable,ClrHatTable,alphadeg) 
    # Roll-Rate Effect, per rad/s
    Clrr	=	ClrHat * (b / (2 * V)) 
    # Yaw-Rate Effect, per rad/s
    CldAr	=	interp1(AlphaTable,CldATable,alphadeg) 
    # Aileron Effect, per rad
    CldRr	=	interp1(AlphaTable,CldRTable,alphadeg) 
    # Rudder Effect, per rad
    CldASr	=	0                   # Asymmetric Spoiler Effect, per rad
    Cl      =	(ClBr*betar + CldRr*u[2]) + Clrr * x[8] + Clpr * x[6] + (CldAr*u[1] + CldASr*u[4]) 
    # Total Rolling-Moment Coefficient
    
    #	Side-Force Coefficient
    CYBr	=	interp1(AlphaTable,CYBetaTable,alphadeg) 
    # Side-Force Slope, per rad
    CYdAr	=	CYdAo               # Aileron Effect, per rad
    CYdRr	=	0.1574 				# Rudder Effect, per rad
    CYdASr	=	0                   # Asymmetric Spoiler Effect, per rad
    CY      =	(CYBr*betar + CYdRr*u[2]) + (CYdAr*u[1] + CYdASr*u[4]) 
    # Total Side-Force Coefficient
    
    #	Yawing Moment Coefficient
    CnBr	=	interp1(AlphaTable,CnBetaTable,alphadeg) 
    # Directional Stability, per rad
    Cnpr	=	CL * (1 + 3 * taperw)/(12 * (1 + taperw)) * (b / (2 * V)) 
    # Roll-Rate Effect, per rad/s
    Cnrr	=	(-2 * (lVT / b) * CnBr - 0.1 * CL **2) * (b / (2 * V)) 
    CnpHat	=	interp1(AlphaTable,CnpHatTable,alphadeg) 
    Cnpr    =   CnpHat * (b / (2 * V)) 
    CnrHat	=	interp1(AlphaTable,CnrHatTable,alphadeg) 
    # Roll-Rate Effect, per rad/s
    Cnrr	=	CnrHat * (b / (2 * V)) 
    # Yaw-Rate Effect, per rad/s
    # Yaw-Rate Effect, per rad/s
    CndAr	=	interp1(AlphaTable,CndATable,alphadeg) 
    # Aileron Effect, per rad
    CndRr	=	interp1(AlphaTable,CndRTable,alphadeg) 
    # Rudder Effect, per rad
    CndASr	=	0 			# Asymmetric Spoiler Effect, per rad
    Cn      =	(CnBr*betar + CndRr*u[2]) + Cnrr * x[8] + Cnpr * x[6]  + (CndAr*u[1]+ CndASr*u[4]) 
    # Total Yawing-Moment Coefficient




def AeroModelMach(x,u,Mach,alphar,betar,V):
    #	Aerodynamic Coefficients of the Aircraft, Thrust Model,
    #	and Geometric and Inertial Properties
    
    
    global m,Ixx,Iyy,Izz,Ixz,S,b,cBar,GEAR,SPOIL
    #	Thrust Properties
    StaticThrust	=	26243.2 	# Static Thrust @ Sea Level, N
    
    
    #	Current Thrust
    airDens,airPres,temp,soundSpeed = Atmos(-x[5]) 
    Thrust  =	u[3] * StaticThrust * (airDens / 1.225) **0.7* (1 - m.exp((-x[5] - 17000) / 2000)) 
    # Thrust at Altitude, N
    #	Current Mach Effects, normalized to Test Condition B (Mach = 0.1734)
    PrFac			=	1 / (m.sqrt(1 - Mach **2) * 1.015) 
    # Prandtl Factor
    WingMach		=	1 / ((1 + m.sqrt(1 + ((ARw/(2 * m.cos(sweepw))) **2)
        * (1 - Mach **2 * m.cos(sweepw)))) * 0.268249) 
    # Modified Helmbold equation
    HorizTailMach	=	1 / ((1 + m.sqrt(1 + ((ARh/(2 * m.cos(sweeph))) **2)
        * (1 - Mach **2 * m.cos(sweeph)))) * 0.294539) 
    # Modified Helmbold equation
    VertTailMach	=	1 / ((1 + sqrt(1 + ((ARv/(2 * m.cos(sweepv))) **2)
        * (1 - Mach **2 * m.cos(sweepv)))) * 0.480338) 
    # Modified Helmbold equation
    
    #	Current Longitudinal Characteristics
    #	====================================
    
    #	Lift Coefficient
    CLo		=	0.1095 				# Zero-AoA Lift Coefficient (B)
    if GEAR >= 1:
        CLo	=	CLo - 0.0192 		# Gear-down correction
    
    if u[5] >= 0.65:
        CLo	=	CLo + 0.5182 		# 38 deg-flap correction
    
    if SPOIL >= 1:
        CLo	=	CLo - 0.1897 		# 42 deg-Symmetric Spoiler correction

    
    CLar	=	5.6575 				# Lift Slope (B), per rad
    if u[5] >= 0.65:
        CLar	=	CLar - 0.0947 

    
    CLqr    =	4.231 * cBar / (2 * V) 
    # Pitch-Rate Effect, per rad/s
    
    CLdSr	=	1.08 				# Stabilator Effect, per rad
    if u[5] >= 0.65:
        CLdSr	=	CLdSr - 0.4802 	# 38ç-flap correction

    CLdEr	=	0.5774 				# Elevator Effect, per rad
    if u[5] >= 0.65:
        CLdEr	=	CLdEr - 0.2665 	# 38 deg-flap correction

    
    CL	=	CLo + (CLar*alphar + CLqr*x[7] + CLdSr*u[6] + CLdEr*u[0])* WingMach 
    # Total Lift Coefficient, w/Mach Correction
    
    #	Drag Coefficient
    CDo		=	0.0255 				# Parasite Drag Coefficient (B)
    if GEAR >= 1:
        CDo	=	CDo + 0.0191 		# Gear-down correction
    
    if u[5] >= 0.65:
        CDo	=	CDo + 0.0836 		# 38 deg-flap correction
    
    if SPOIL >= 1:
        CDo	=	CDo + 0.0258 		# 42 deg-Symmetric Spoiler correction
    
    
    epsilon	=	0.0718 				# Induced Drag Factor
    if u[5] >= 0.65:
        epsilon	=	0.079 			# 38 deg-flap correction
    
    
    CD	=	CDo * PrFac + epsilon * CL **2 
    # Total Drag Coefficient, w/Mach Correction
    
    #	Pitching Moment Coefficient
    Cmo		=	0 					# Zero-AoA Moment Coefficient (B)
    if GEAR >= 1:
        Cmo	=	Cmo + 0.0255 		# Gear-down correction
    
    if u[5]>= 0.65:
        Cmo	=	Cmo - 0.058 		# 38 deg-flap correction
    
    if SPOIL >= 1:
        Cmo	=	Cmo - 0.0154 		# 42 deg-Symmetric Spoiler correction
    
    
    Cmar	=	-1.231 				# Static Stability (B), per rad
    if u[5] >= 0.65:
        Cmar	=	Cmar + 0.0138 
    
    
    Cmqr    =	 -18.8 * cBar / (2 * V) 
    # Pitch-Rate + Alpha-Rate Effect, per rad/s
    
    CmdSr	=	-2.291 				# Stabilator Effect, per rad
    if u[5] >= 0.65:
        CmdSr	=	CmdSr + 0.121 	# 38 deg-flap correction
    
    
    CmdEr	=	-1.398 				# Elevator Effect, per rad
    if u[5] >= 0.65:
        CmdEr	=	CmdEr + 0.149 	# 38 deg-flap correction
   
    
    Cm	=	Cmo + (Cmar*alphar + Cmqr*x[7] + CmdSr*u[6] + CmdEr*u[0]) * HorizTailMach 
    # Total Pitching Moment Coefficient, w/Mach Correction
    
    #	Current Lateral-Directional Characteristics
    #	===========================================
    
    #	Side-Force Coefficient
    CYBr	=	-0.7162 			# Side-Force Slope (B), per rad
    if u[5]>= 0.65:
        CYBr	=	CYBr + 0.0826 
    
    
    CYdAr	=	-0.00699 			# Aileron Effect, per rad
    
    CYdRr	=	0.1574 				# Rudder Effect, per rad
    if u[5] >= 0.65:
        CYdRr	=	CYdRr - 0.0093 	# 38 deg-flap correction
    
    
    CYdASr	=	0.0264 				# Asymmetric Spoiler Effect, per rad
    if u[5]>= 0.65:
        CYdASr	=	CYdASr + 0.0766 
        # 38 deg-flap correction
    
    
    CY	=	(CYBr*betar + CYdRr*u[2]) * VertTailMach + (CYdAr*u[1] + CYdASr*u[4]) * WingMach 
    # Total Side-Force Coefficient, w/Mach Correction
    
    #	Yawing Moment Coefficient
    CnBr	=	0.1194 				# Directional Stability (B), per rad
    if u[5] >= 0.65:
        CnBr	=	CnBr - 0.0092 
    
    
    Cnpr	=	CL * (1 + 3 * taperw)/(12 * (1 + taperw)) * (b / (2 * V)) 
    # Roll-Rate Effect, per rad/s
    
    Cnrr	=	(-2 * (lvt / b) * CnBr * VertTailMach - 0.1 * CL **2)* (b / (2 * V)) 
    # Yaw-Rate Effect, per rad/s
    
    CndAr	=	0                       # Aileron Effect, per rad
    if u[5] >= 0.65:
        CndAr	=	CndAr + 0.0028 
    
    
    CndRr	=	-0.0713                 # Rudder Effect, per rad
    if u[5] >= 0.65:
        CndRr	=	CndRr - 0.0185      # 38 deg-flap correction
    
    
    CndASr	=	-0.0088                 # Asymmetric Spoiler Effect, per rad
    if u[5] >= 0.65:
        CndASr	=	CndASr - 0.0106 
        # 38 deg-flap correction
    
    
    Cn	=	(CnBr*betar + CndRr*u[2]) * VertTailMach+ Cnrr * x[8] + Cnpr * x[6] + (CndAr*u[1] + CndASr*u[4]) * WingMach 
    # Total Yawing-Moment Coefficient, w/Mach Correction
    
    #	Rolling Moment Coefficient
    ClBr	=	-0.0918                 # Dihedral Effect (B), per rad
    if u[5] >= 0.65:
        ClBr	=	ClBr - 0.0092 
    
    
    Clpr	=	-CLar * (1 + 3 * taperw)/(12 * (1 + taperw))* (b / (2 * V)) 
    # Roll-Rate Effect, per rad/s
    
    Clrr	=	(CL * (1 + 3 * taperw)/(12 * (1 + taperw))* ((Mach * m.cos(sweepw)) **2 - 2) / ((Mach * m.cos(sweepw)) **2 - 1)) * (b / (2 * V)) 
    # Yaw-Rate Effect, per rad/s
    
    CldAr	=	0.1537                  # Aileron Effect, per rad
    if u[5] >= 0.65:
        CldAr	=	CldAr + 0.01178 
    
    
    CldRr	=	0.01208 				# Rudder Effect, per rad
    if u[5] >= 0.65:
        CldRr	=	CldRr + 0.01115 	# 38 deg-flap correction
    
    
    CldASr	=	-0.01496 				# Asymmetric Spoiler Effect, per rad
    if u[5] >= 0.65:
        CldASr	=	CldASr - 0.02376 
        # 38 deg-flap correction
   
    
    Cl	=	(ClBr*betar + CldRr*u[2]) * VertTailMach+ Clrr * x[8] + Clpr * x[6]+ (CldAr*u[1] + CldASr*u[4]) * WingMach 
    # Total Rolling-Moment Coefficient, w/Mach Correction


def AeroModelUser(x,u,Mach,alphar,betar,V):
    
    
    global m,Ixx,Iyy,Izz,Ixz,S,b,cBar
    
    
    cBar            =   3.03 
    b               =   10 
    S               =   27.77 
    lHT             =   5.2 
    lVT             =   3.9 
    StaticThrust    =   49000 
    
    [airDens,airPres,temp,soundSpeed] = Atmos(-x(6)) 
    Thrust   =   u(4)*StaticThrust*(airDens/1.225)*(1 - exp((-x(6) - 18000)/2000)) 
    
    #   Mach Number Effect on All Incompressible Coefficients
    Prandtl			=	1 / (sqrt(1 - Mach **2)) 	# Prandtl Factor
    
    #	Current Longitudinal Characteristics
    #	====================================
    #	Lift Coefficient
    CLo     =   0 
    CLa     =   4.92 
    CLqhat  =   2.49 
    CLq     =	CLqhat*cBar/(2*V) 
    CLdE    =   0.72 
    CLdS    =   CLdE 
    
    #       Total Lift Coefficient, w/Mach Correction
    CL      =	(CLo + CLa*alphar + CLq*x(8) + CLdS*u(7) + CLdE*u(1))*Prandtl 
    
    #   Drag Coefficient
    CDo     =   0.019 
    Epsilon =   0.093*Prandtl 
    
    #       Total Drag Coefficient, w/Mach Correction
    CD      =	CDo*Prandtl + Epsilon*CL **2 
    
    #	Pitching Moment Coefficient
    StaticMargin    =   0.2 
    Cmo     =   0 
    Cma 	=	-CLa*StaticMargin 
    Cmqhat  =   -4.3 
    Cmq     =	Cmqhat*cBar/(2*V) 
    CmV     =   0 
    CmdE	=	-1.25 
    CmdS 	=	CmdE 
    
    #       Total Pitching Moment Coefficient, w/Mach Correction
    Cm      =	(Cmo + Cma*alphar + Cmq*x(8) + CmdS*u(7) + ...
        CmdE*u(1))*Prandtl 
    
    #	Current Lateral-Directional Characteristics
    #	===========================================
    #	Side-Force Coefficient
    CYo     =   0 
    CYb     =   -0.5 
    CYphat  =   0 
    CYp     =   CYphat*(b/(2*V)) 
    CYrhat  =   0 
    CYr     =   CYrhat*(b/(2*V)) 
    CYdA	=	0 
    CYdR	=	0.04 
    
    #       Total Side-Force Coefficient, w/Mach Correction
    CYo     =   0 
    CY      =	(CYo + CYb*betar + CYdR*u(3) + CYdA*u(2) + CYp*x(7) ...
        + CYr*x(9))*Prandtl 
    
    #	Rolling Moment Coefficient
    Clo     =   0 
    Clb  	=	-0.066 
    Clphat  =   -0.5 
    Clp 	=	Clphat*(b/(2*V)) 
    Clrhat  =   -0.5 
    Clr 	=	Clrhat* (b/(2*V)) 
    CldA 	=	0.12 
    CldR 	=	0.03 
    
    #       Total Rolling-Moment Coefficient, w/Mach Correction
    Cl      =	(Clo + Clb*betar + Clp * x(7) + Clr * x(9) ...
        + CldA*u(2) + CldR*u(3))* Prandtl 
    
    #	Yawing Moment Coefficient
    Cno     =   0 
    CnBeta  =	0.37 
    Cnphat  =   -0.06 
    Cnp 	=	Cnphat*(b/(2*V)) 
    Cnrhat  =   -0.5 
    Cnr 	=	Cnrhat*(b/(2*V)) 
    CndA 	=	0 
    CndR 	=	-0.2 
    
    # Total Yawing-Moment Coefficient, w/Mach Correction
    Cn      =	(Cno + CnBeta*betar + Cnp*x(7) + Cnr*x(9) ...
        + CndA*u(2) + CndR*u(3))*Prandtl 




def Atmos(geomAlt):
          
    
    Z   = [-1000,0,2500,5000,10000,11100,15000,20000,47400,51000]
    H   = [-1000,0,2499,4996,9984,11081,14965,19937,47049,50594]
    ppo = [1,1,0.737,0.533,0.262,0.221,0.12,0.055,0.0011,0.0007]
    rro = [1,1,0.781,0.601,0.338,0.293,0.159,0.073,0.0011,0.0007]
    T = [288.15,288.15,271.906,255.676,223.252,216.65,216.65,216.65, 270.65,270.65]
    a = [340.294,340.294,330.563,320.545,299.532,295.069,295.069,295.069,329.799,329.799]
    R = 6367435 # Mean radius of the earth, m
    Dens = 1.225  # Air density at sea level, Kg/m^3
    Pres	= 101300 # Air pressure at sea level, N/m^2
    
    # Geopotential Altitude, m
    #airDens = Dens
    #airPres = Pres
    
    geopAlt	=	R * geomAlt / (R + geomAlt)
    
    #	Linear Interpolation in Geopotential Altitude 
    #	for Temperature and Speed of Sound
    temp		=	np.interp(geopAlt,Z,T)
    soundSpeed	=	np.interp(geopAlt,Z,a)   #ajustado
    
    #	Exponential Interpolation in Geometric Altitude
    #	for Air Density and Pressure for ii in range(0,gt):
    for k in range(1,9):
        if geomAlt <= Z[k]:
            betap	=	math.log(ppo[k] / ppo[k-1]) / (Z[k] - Z[k-1]) #Ajustado
            betar	=	math.log(rro[k] / rro[k-1]) / (Z[k] - Z[k-1]) #Ajustado
            airPres	=	Pres * ppo[k-1] * np.exp(betap * (geomAlt - Z[k-1])) #Ajustado
            airDens	=	Dens * rro[k-1] * np.exp(betar * (geomAlt - Z[k-1])) #Ajustado
        
            return airDens,airPres,temp,soundSpeed
    
    
    [airDens,airPres,temp,soundSpeed] = [-geomAlt,-geomAlt,-geomAlt,-geomAlt]
    
    return airDens,airPres,temp,soundSpeed / (R + geomAlt)
    
    #	Linear Interpolation in Geopotential Altitude 
    #	for Temperature and Speed of Sound
    temp		=	np.interp(geopAlt,Z,T)
    soundSpeed	=	np.interp(geopAlt,Z,a)   #ajustado
    
    #	Exponential Interpolation in Geometric Altitude
    #	for Air Density and Pressure for ii in range(0,gt):
    for k in range(1,9):
        if geomAlt <= Z[k]:
            betap	=	math.log(ppo[k] / ppo[k-1]) / (Z[k] - Z[k-1]) #Ajustado
            betar	=	math.log(rro[k] / rro[k-1]) / (Z[k] - Z[k-1]) #Ajustado
            airPres	=	Pres * ppo[k-1] * np.exp(betap * (geomAlt - Z[k-1])) #Ajustado
            airDens	=	Dens * rro[k-1] * np.exp(betar * (geomAlt - Z[k-1])) #Ajustado
        
    return airDens,airPres,temp,soundSpeed