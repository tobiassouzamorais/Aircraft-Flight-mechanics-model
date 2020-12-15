# Simulation

def simul():
    import numpy as np
    import aircraft_parameters as ap  
    import aero as aero
#    import Trim_airplane as Trim_airplane
#    import EoM_mod as EoM_mod
    
    
    uInc = np.zeros(7)
    u_total = np.array(np.zeros((7)))
    x_t=np.zeros((12,len(ap.t)))
    u_t=np.zeros((7,len(ap.t)))
    Aero_par=np.zeros((7,len(ap.t)))
    gust_t=aero.gust()
    for i in range(0,len(ap.t)):
        if i==0:
            x, u, v = Trim_airplane()  

            x_t=np.zeros((12,len(ap.t)))
            u_t=np.zeros((7,len(ap.t)))
            uInc[0]	=	np.interp(ap.t[i],ap.tuHis,ap.deluHis[:,0])
            uInc[1]	=	np.interp(ap.t[i],ap.tuHis,ap.deluHis[:,1])
            uInc[2]	=	np.interp(ap.t[i],ap.tuHis,ap.deluHis[:,2])
            uInc[3]	=	np.interp(ap.t[i],ap.tuHis,ap.deluHis[:,3])
            uInc[4]	=	np.interp(ap.t[i],ap.tuHis,ap.deluHis[:,4])
            uInc[5]	=	np.interp(ap.t[i],ap.tuHis,ap.deluHis[:,5])
            uInc[6]	=	np.interp(ap.t[i],ap.tuHis,ap.deluHis[:,6])
            x_t[:,i]=x
            u_t[:,i]=u #+uInc
        else:
            uInc[0]	=	np.interp(ap.t[i],ap.tuHis,ap.deluHis[:,0])
            uInc[1]	=	np.interp(ap.t[i],ap.tuHis,ap.deluHis[:,1])
            uInc[2]	=	np.interp(ap.t[i],ap.tuHis,ap.deluHis[:,2])
            uInc[3]	=	np.interp(ap.t[i],ap.tuHis,ap.deluHis[:,3])
            uInc[4]	=	np.interp(ap.t[i],ap.tuHis,ap.deluHis[:,4])
            uInc[5]	=	np.interp(ap.t[i],ap.tuHis,ap.deluHis[:,5])
            uInc[6]	=	np.interp(ap.t[i],ap.tuHis,ap.deluHis[:,6])
            u_total = u + uInc
            xdot1, CD, CL, CY, Cl, Cm, Cn, Thrust = EoM_mod(ap.t[i-1],x_t[:,i-1],u_total,gust_t[i])
            k1 = ap.h*xdot1
            xdot2, CD, CL, CY, Cl, Cm, Cn, Thrust = EoM_mod(ap.t[i-1]+ap.h/2, x_t[:,i-1]+k1/2,u_total,gust_t[i])
            k2 = ap.h*xdot2
            xdot3, CD, CL, CY, Cl, Cm, Cn, Thrust = EoM_mod(ap.t[i-1]+ap.h/2, x_t[:,i-1]+k2/2,u_total,gust_t[i])
            k3 = ap.h*xdot3
            xdot4, CD, CL, CY, Cl, Cm, Cn, Thrust = EoM_mod(ap.t[i-1]+ap.h, x_t[:,i-1]+k3,u_total,gust_t[i])
            k4 = ap.h*xdot4
            x_t[:,i]= x_t[:,i-1] + (k1+2*k2+2*k3+k4)/6
            u_t[:,i]=u_total    
            Aero_par[:,i] =   np.array([CD, CL, CY, Cl, Cm, Cn, Thrust]) 
            # print(Aero_par[:,i])
    time=ap.t
            
    return time, x_t, u_t, Aero_par



def EoM_mod(t,x,u,gust):
    import numpy as np
    import aircraft_parameters as ap  
    import aero as aero
    import Atmos as At
    
    
    D2R =   np.pi/180
    R2D =   180/np.pi
    
        #  Select Aerodynamic Model

        #if MODEL == 0
        #    AeroModel   =   @AeroModelAlpha;
        #end
        #if MODEL == 1
        #   AeroModel   =   @AeroModelMach;
        #end
        #if MODEL == 2
        #   AeroModel   =   @AeroModelUser;
        #end

        #[value,isterminal,direction] = event(x)
    
    #	Atmospheric State

    airDens,airPres,temp,soundSpeed = At.Atmo(-x[5])     
        #if self.airDens == -1:
            #xdot = x*(-1)
            #return xdot,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,self.airDens

    #Earth-to-Body-Axis Transformation Matrix
    HEB = aero.DCM(x[9],x[10],x[11])  #Ajustado
     
    # Body-Axis Wind Field
    windb	=	aero.WindField(-x[5],x[9],x[10],x[11],gust)      #Ajustado
    #print(windb)
        # Body-Axis Gravity Components
    gb		=	np.dot(HEB, np.transpose(np.array([0,0,9.80665]))) #Ajustado
    #print(gb)
        #	Air-Relative Velocity Vector
    x[0]    =   np.max([x[0],0])       #  Limit axial velocity to >= 0 m/s
    
    Va		=	np.array([x[0],x[1],x[2]]) + np.array(windb)
    V		=	np.sqrt(np.dot(np.matrix.transpose(Va),Va)) 
    # print(Va)
    print(V)                    
    alphar	=	np.arctan(Va[2] / np.absolute(Va[0]))
        #alphar  =   min(alphar, (pi/2 - 1e-6));  %   Limit angle of attack to <= 90 deg
    
    alpha 	=	R2D * alphar
    
    betar	= 	np.arcsin(Va[1] / V)
    beta	= 	R2D * betar
    Mach	= 	V / soundSpeed
    qbar	=	0.5 * airDens * V**2
    
    CL, CD, CY, Cm, Cn, Cl, Thrust =	aero.AeroModelUser(x, u, Mach, alphar, betar, V)
    
    qbarS	=	qbar * ap.S

    CX	=	-CD * np.cos(alphar) + CL * np.sin(alphar)  # Body-axis X coefficient
    CZ	= 	-CD * np.sin(alphar) - CL * np.cos(alphar)  # Body-axis Z coefficient

        #	State Accelerations
   
    Xb  =	(CX * qbarS + Thrust) / ap.Mass
    Yb  =	CY * qbarS / ap.Mass
    Zb  =	CZ * qbarS / ap.Mass
    Lb  =	Cl * qbarS * ap.b
    Mb  =	Cm * qbarS * ap.cBar
    Nb  =	Cn * qbarS * ap.b
    #    nz  =	-Zb / 9.80665                                #  Normal load factor

    #	Dynamic Equations
    xd1 = Xb + gb[0] + x[8] * x[1] - x[7] * x[2]        # u_ponto (3.2-47)
    xd2 = Yb + gb[1] - x[8] * x[0] + x[6] * x[2]        # v_ponto (3.2-48)
    xd3 = Zb + gb[2] + x[7] * x[0] - x[6] * x[1]        # w_ponto (3.2-49)

                      
    y	=	np.dot(np.matrix.transpose(HEB),[x[0],x[1],x[2]])
    xd4	=	y[0]                                        # Velocidade inercial x
    xd5	=	y[1]                                        # Velocidade inercial y
    xd6	=	y[2]                                        # Velocidade inercial z         
    xd7 = 	(ap.Izz * Lb + ap.Ixz * Nb - (ap.Ixz * (ap.Iyy - ap.Ixx - ap.Izz) * x[6] + (ap.Ixz**2 + ap.Izz * (ap.Izz - ap.Iyy)) * x[8]) * x[7]) / (ap.Ixx * ap.Izz - ap.Ixz**2)  # p_ponto(3.2-51)
    xd8 = 	(Mb - (ap.Ixx - ap.Izz) * x[6] * x[8] - ap.Ixz * (x[6]**2 - x[8]**2)) / ap.Iyy   # q_ponto(3.2-52)
    xd9 =	(ap.Ixz * Lb + ap.Ixx * Nb + (ap.Ixz * (ap.Iyy - ap.Ixx - ap.Izz) * x[8] + (ap.Ixz**2 + ap.Ixx * (ap.Ixx - ap.Iyy)) * x[6]) * x[7]) / (ap.Ixx * ap.Izz - ap.Ixz**2) # r_ponto(3.2-53)

    cosPitch	=	np.cos(x[10])
                      
    if np.abs(cosPitch)	<=	0.00001:
       cosPitch	=	0.00001 * np.sign(cosPitch)

    tanPitch	=	np.sin(x[10]) / cosPitch

    xd10	=	x[6] + (np.sin(x[9]) * x[7] + np.cos(x[9]) * x[8]) * tanPitch   # phi_ponto(3.2-57)
    xd11	=	np.cos(x[9]) * x[7] - np.sin(x[9]) * x[8]                     # theta_ponto(3.2-58)
    xd12	=	(np.sin(x[9]) * x[7] + np.cos(x[9]) * x[8]) / cosPitch
        
        #States in the RL
    xdot = np.array([xd1, xd2,xd3, xd4, xd5, xd6, xd7, xd8, xd9, xd10, xd11, xd12])

    return xdot, CD, CL, CY, Cl, Cm, Cn, Thrust



def Trim_airplane():
    import numpy as np
    import Atmos as At
        
    u = np.zeros(7)
    x = np.zeros(12)
    OptParam = np.zeros(3)
        

    MIC     =   0.75   #   Initial Mach Number
    hm      =   10500   #   Initial Altitude, m
        
    [airDens,airPres,temp,soundSpeed] = At.Atmo(hm)

    hft     =   hm/0.3048       #   Initial Altitude, ft
    VIC     =   MIC*soundSpeed  #   Initial True Airspeed, m/s
    VTASkt  =   VIC*1.94384     #   Initial True Airspeed, kt
    VKIAS   =   VTASkt*np.sqrt(airDens/1.225)  #    Initial Indicated Airspeed, kt

    hm          =   hft * 0.3048           # Altitude above Sea Level, m
    VmsIAS      =   VKIAS * 0.514444       # Indicated Airspeed, m/s



    #   US Standard Atmosphere, 1976, Table Lookup for I.C.

    [airDens,airPres,temp,soundSpeed] = At.Atmo(hm)

    qBarSL  =   0.5*1.225*VmsIAS**2           # Dynamic Pressure at sea level, N/m^2
    v       =   np.sqrt(2*qBarSL/airDens) 	  # True Airspeed, TAS, m/s
        
        # It was defined in MATLAB 
        
    OptParam[0] = -0.0209
    OptParam[1] = 0.5847
    OptParam[2] = 0.0265
        
    u[0]=  0
    u[1]=    0
    u[2]=    0
    u[3]=    0.58466061293618
    u[4]=    0
    u[5]=    0
    u[6]=    -0.020876191027685

    x[0]   =    223.10093995161
    x[1]   =    0
    x[2]   =    5.9179
    x[3]   =    0
    x[4]   =    0
    x[5]   =    -10500
    x[6]   =    0
    x[7]   =    0
    x[8]   =    0
    x[9]   =    0
    x[10]   =   0.026519
    x[11]   =   0
    
    return x,u,v
    
