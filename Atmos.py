def Atmo(geomAlt):
    import numpy as np   
    import math as m    
    
    Z   = np.array([-1000,0,2500,5000,10000,11100,15000,20000,47400,51000])
    H   = np.array([-1000,0,2499,4996,9984,11081,14965,19937,47049,50594])
    ppo = np.array([1,1,0.737,0.533,0.262,0.221,0.12,0.055,0.0011,0.0007])
    rro = np.array([1,1,0.781,0.601,0.338,0.293,0.159,0.073,0.0011,0.0007])
    T = np.array([288.15,288.15,271.906,255.676,223.252,216.65,216.65,216.65, 270.65,270.65])
    a = np.array([340.294,340.294,330.563,320.545,299.532,295.069,295.069,295.069,329.799,329.799])
    R = 6367435 # Mean radius of the earth, m
    Dens = 1.225  # Air density at sea level, Kg/m^3
    Pres	= 101300 # Air pressure at sea level, N/m^2
    
    # Geopotential Altitude, m
    #airDens = Dens
    #airPres = Pres
    
    geopAlt	=	R * (geomAlt) / (R + (geomAlt))
   

    
    #	Linear Interpolation in Geopotential Altitude 
    #	for Temperature and Speed of Sound
    temp		=	(np.interp(geopAlt,Z,T))
    soundSpeed	=	(np.interp(geopAlt,Z,a))   #ajustado
    
    #	Exponential Interpolation in Geometric Altitude
    #	for Air Density and Pressure for ii in range(0,gt):
    ii = 1
    for k in range(1,10):
        if geomAlt <= Z[k] and ii == 1:
            betap	=	m.log(ppo[k] / ppo[k-1]) / (Z[k] - Z[k-1]) #Ajustado
            betar	=	m.log(rro[k] / rro[k-1]) / (Z[k] - Z[k-1]) #Ajustado
            airPres	=	(Pres * ppo[k-1] * np.exp(betap * (geomAlt - Z[k-1]))) #Ajustado
            airDens	=	(Dens * rro[k-1] * np.exp(betar * (geomAlt - Z[k-1]))) #Ajustado
            ii = 0
    return airDens, airPres, temp, soundSpeed
