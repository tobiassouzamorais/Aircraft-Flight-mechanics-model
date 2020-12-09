# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import math as m

# global GEAR CONHIS SPOIL u x V uInc tuHis deluHis TrimHist SMI MODEL RUNNING

#	DEFINITION OF THE STATE VECTOR
#   With Euler Angle DCM option (QUAT = 0):
#		x(0)    = 		Body-axis x inertial velocity, ub, m/s
#		x(1)    =		Body-axis y inertial velocity, vb, m/s
#		x(2)    =		Body-axis z inertial velocity, wb, m/s
#		x(3)    =		North position of center of mass WRT Earth, xe, m
#		x(4)    =		East position of center of mass WRT Earth, ye, m
#		x(5)    =		Negative of c.m. altitude WRT Earth, ze = -h, m
#		x(6)    =		Body-axis roll rate, pr, rad/s
#		x(7)    =		Body-axis m.pitch rate, qr, rad/s
#		x(8)    =		Body-axis yaw rate, rr,rad/s
#		x(9)   =		Roll angle of body WRT Earth, phir, rad
#		x(10)   =		m.pitch angle of body WRT Earth, thetar, rad
#		x(11)   =		Yaw angle of body WRT Earth, psir, rad
#   With Quaternion DCM option (QUAT = 1):
#		x(0)    = 		Body-axis x inertial velocity, ub, m/s
#		x(1)    =		Body-axis y inertial velocity, vb, m/s
#		x(2)    =		Body-axis z inertial velocity, wb, m/s
#		x(3)    =		North position of center of mass WRT Earth, xe, m
#		x(4)    =		East position of center of mass WRT Earth, ye, m
#		x(5)    =		Negative of c.m. altitude WRT Earth, ze = -h, m
#		x(6)    =		Body-axis roll rate, pr, rad/s
#		x(7)    =		Body-axis m.pitch rate, qr, rad/s
#		x(8)    =		Body-axis yaw rate, rr,rad/s
#		x(9)   =		q1, x Component of quaternion
#		x(10)   =		q2, y Component of quaternion
#		x(11)   =		q3, z Component of quaternion
#		x(12)   =		q4, cos(Euler) Component of quaternion

#	DEFINITION OF THE CONTROL VECTOR
#		u(0)    = 		Elevator, dEr, rad, positive: trailing edge down
#		u(1)    = 		Aileron, dAr, rad, positive: left trailing edge down
#		u(2)    = 		Rudder, dRr, rad, positive: trailing edge left
#		u(3)    = 		Throttle, dT, #
#		u(4)    =		Asymmetric Spoiler, dASr, rad
#		u(5)    =		Flap, dFr, rad
#		u(6)    =		Stabilator, dSr, rad

#   ======================================================================
#	USER INPUTS
#   ===========
#   --  Flags define which analyses or IC/Inputs will be engaged
#	======================================================
#	FLIGHT Flags (1 = ON, 0 = OFF)

MODEL   =   2       # Aerodynamic model selection
# 0: Incompressible flow, high angle of attack
# 1: Compressible flow, low angle of attack
# 2: User-Defined model

TRIM    = 	1 		# Trim flag (= 1 to calculate trim @ I.C.)
LINEAR  = 	0 		# Linear model flag (= 1 to calculate and store F and G)
SIMUL   =	2 		# Flight path flag (= 1 for nonlinear simulation = 2 for external solver)
CONHIS  =	1 		# Control history ON (= 1) or OFF (= 0)
RUNNING =   0       # internal flag, -

GEAR    = 	0 		# Landing gear DOWN (= 1) or UP (= 0)
SPOIL   =	0 		# Symmetric Spoiler DEPLOYED (= 1) or CLOSED (= 0)
dF      = 	0 		# Flap setting, deg

##   Starting and Final Times for Simulation
#   =======================================
ti      = 	0       # Initial time for simulation, sec
tf      =	60     # Final time for simulation, sec

MIC     =   0.75     #   Initial Mach Number
hm      =   10500    #   Initial Altitude, m

airDens,airPres,temp,soundSpeed = Atmos(hm) 

hft     =   hm/0.3048        #   Initial Altitude, ft
VIC     =   MIC*soundSpeed   #   Initial True Airspeed, m/s
VTASkt  =   VIC*1.94384      #   Initial True Airspeed, kt
VKIAS   =   VTASkt*m.sqrt(airDens/1.225)   #   Initial Indicated Airspeed, kt

hm          =   hft * 0.3048     # Altitude above Sea Level, m
VmsIAS      =   VKIAS * 0.514444   # Indicated Airspeed, m/s


#   US Standard Atmosphere, 1976, Table Lookup for I.C.
airDens,airPres,temp,soundSpeed = Atmos(hm) 
qBarSL  =   0.5*1.225*VmsIAS **2   # Dynamic Pressure at sea level, N/m **2
V       =   m.sqrt(2*qBarSL/airDens) 	# True Airspeed, TAS, m/s
TASms   =   V 


D2R =   m.pi/180 
R2D =   180/m.pi 

#	Alphabetical List of Additional Initial Conditions
#   ==================================================
#   These valuees subsume the hft and VKIAS entered above

alpha   =	0       # Angle of attack, deg (relative to air mass)
beta    =	0       # Sideslip angle, deg (relative to air mass)
dA      =	0       # Aileron angle, deg
dAS     =	0       # Asymmetric spoiler angle, deg
dE      =	0       # Elevator angle, deg
dR      =	0       # Rudder angle, deg
dS      = 	0       # Stabilator setting, deg
dT      = 	0       # Throttle setting, # / 100
hdot    =	0       # Altitude rate, m/s
p       =	0       # Body-axis roll rate, deg/s
phi     =	0       # Body roll angle wrt earth, deg
psi     =	0       # Body yaw angle wrt earth, deg
q       =	0       # Body-axis m.pitch rate, deg/sec
r       =	0       # Body-axis yaw rate, deg/s
SMI     =	0       # Static margin increment due to center-of-mass variation

# from reference, #/100
theta   =	alpha   # Body m.pitch angle wrt earth, deg [theta = alpha if hdot = 0]
xe      =	0       # Initial longitudinal position, m
ye      = 	0       # Initial lateral position, m
ze      = 	-hm     # Initial vertical position, m [h: + up, z: + down]

#	Initial Conditions Derived from Prior Initial Conditions
gamma	=	R2D * m.atan(hdot / m.sqrt(V **2 - hdot **2)) 
# Inertial Vertical Flight Path Angle, deg
qbar	= 	0.5 * airDens * V **2 
# Dynamic Pressure, N/m **2
IAS		=	m.sqrt(2 * qbar / 1.225) 
# Indicated Air Speed, m/s
Mach	= 	V / soundSpeed 
uInc    =   [] 

#   Test Inputs at Initial Condition
#   ================================
#	Initial Control Perturbation (Test Inputs: rad or percent)
delu	=	np.array([0,0, 0, 0, 0, 0, 0]).T
#	Initial State Perturbation (Test Inputs: m, m/s, rad, or rad/s)
delx	=	np.array([0, 0, 0,
                      0, 0, 0,
                      0, 0, 0,
                      0, 0, 0]).T

#	Control Perturbation History (Test Inputs: rad or 100#)
#   =======================================================
#   Each control effector represented by a column
#	Each row contains control increment delta-u(t) at time t:

tuHis	=	np.array([0,5,6,7,60]) 

# print('Columns:  Elements of the Control Vector')
# print('Rows:     Value at time, t0, t1, ...')
deluHis	=	np.array([0, 0, 0, 0, 0, 0, 0,
                     3*m.pi/180, 0,  0, 0, 0, 0, 0,
                     3*m.pi/180, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0]) 


#	DEFINITION OF THE CONTROL VECTOR
#		u(0)    = 		Elevator, dEr, rad, positive: trailing edge down
#		u(1)    = 		Aileron, dAr, rad, positive: left trailing edge down
#		u(2)    = 		Rudder, dRr, rad, positive: trailing edge left
#		u(3)    = 		Throttle, dT, #
#		u(4)    =		Asymmetric Spoiler, dASr, rad
#		u(5)    =		Flap, dFr, rad
#		u(6)    =		Stabilator, dSr, rad


#	State Vector and Control Initialization, rad
phir	=	phi * D2R 
thetar	=	theta * D2R 
psir	=	psi * D2R 

windb	=	WindField(-ze,phir,thetar,psir) 
alphar	=	alpha * D2R 
betar	=	beta * D2R 

x	=	np.array([V * m.cos(alphar) * m.cos(betar) - windb[0],
    V * m.sin(betar) - windb[1],
    V * m.sin(alphar) * m.cos(betar) - windb[2],
    xe,
    ye,
    ze,
    p * D2R,
    q * D2R,
    r * D2R,
    phir,
    thetar,
    psir]) 

u	=	np.array([dE * D2R,
    dA * D2R,
    dR * D2R,
    dT,
    dAS * D2R,
    dF * D2R,
    dS * D2R]) 

#	Trim Calculation (for Steady Level Flight at Initial V and h)
#   =============================================================
#   Euler Angles used in trim calculation
#   Trim Parameter Vector (OptParam):
#		1 = Stabilator, rad
#		2 = Throttle, #
#		3 = m.pitch Angle, rad

if TRIM >= 1:
     print('  ')
     print('TRIM Stabilator, Thrust, and m.pitch Angle')
     print('========================================')
     OptParam        =   [] 
     TrimHist        =   [] 
     #       Arbitrary starting values (user-selected, e.g., best guess at solution)
     InitParam		=	np.array([0.0369,0.1892,0.0986]) 
     
     options =   optimset('TolFun',1e-10) 
     OptParam,J,ExitFlag,Output  =	fminsearch(@TrimCost,InitParam,options) 
     
     u	=	np.array([u[0],
     u[1],
     u[2],
     OptParam[1],
     u[4],
     u[5],
     OptParam[0]])
               
     x	=	np.array([V * m.cos(OptParam[2]),
     x[1],
     V * m.sin(OptParam[2]),
     x[3],
     x[4],
     x[5],
     x[6],
     x[7],
     x[8],
     x[9],
     OptParam[2],
     x[11]]) 
     


#	Stability-and-Control Derivative Calculation
#   ============================================
if LINEAR >= 1: 
    print('  ')
    print('Generate and Save LINEAR MODEL')
    print('==============================')
    thresh	=	np.array([.1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1]) 
    xj		=	np.array([x ,u]).T 
    uTemp   =   u   # 'numjac' modifies 'u'  reset 'u' after the call
    xdotj		=	LinModel(ti,xj) 
    dFdX,fac	=	numjac(@LinModel,ti,xj,xdotj,thresh,[],0) 
    u           =   uTemp 
    Fmodel		=	dFdX[1:12,1:12]
    Gmodel		=	dFdX[1:12,13:19]


#	Flight Path Simulation, with Quaternion Option
#   ==============================================
if SIMUL >= 1:
    RUNNING =   1 
    tspan	=	np.array([ti,tf]) 
    xo		=	x + delx 
    u		=	u + delu 
    
    if SIMUL== 1:
        options =   odeset('Events',@event,'RelTol',1e-7,'AbsTol',1e-7) 
        tic
        [t,x]	=	ode45(@EoM,tspan,xo,options) 
        toc
        kHis	=	np.length(t) 
        #       Plot Control History
        
    if SIMUL== 2:
            
        f = @EoM_mod 
        h = 0.1   # Define Step Size
        t = ti:h:tf 
        y = np.zeros(12,numel(t)) 
        y[:,0] = x   # condição trimada
        
        for i in range (1,numel(t)):
            k1 = h*f(t[i-1],y[:,i-1],u) 
            k2 = h*f(t[i-1]+h/2, y[:,i-1]+k1/2,u) 
            k3 = h*f(t[i-1]+h/2, y[:,i-1]+k2/2,u) 
            k4 = h*f(t[i-1]+h, y[:,i-1]+k3,u) 
            y[:,i] = y[:,i-1] + (k1+2*k2+2*k3+k4)/6 
       
        x = y.T 
        kHis	=	np.length(t) 


# figure
# subplot(2,2,1)
# plot(tuHis,R2D*deluHis(:,1) + R2D*delu(1), tuHis, R2D*deluHis(:,7) + R2D*delu(7))
# xlabel('Time, s'), ylabel('Elevator (blue), Stabilator (green), deg'), grid
# title('m.pitch Test Inputs'), legend('Elevator, dE', 'Stabilator, dS')
# subplot(2,2,2)
# plot(tuHis,R2D*deluHis(:,2) + R2D*delu(2), tuHis, R2D*deluHis(:,3) + R2D*delu(3), tuHis, R2D*deluHis(:,5) + R2D*delu(5))
# xlabel('Time, s'), ylabel('Aileron (blue), Rudder (green), Asymmetric Spoiler (red), deg'), grid
# title('Lateral-Directional Test Inputs'), legend('Aileron, dA', 'Rudder, dR', 'Asymmetric Spoiler, dAS')
# subplot(2,2,3)
# plot(tuHis, deluHis(:,4) + delu(4))
# xlabel('Time, s'), ylabel('Throttle Setting'), grid
# title('Throttle Test Inputs')
# subplot(2,2,4)
# plot(tuHis,R2D*deluHis(6) + R2D*delu(6))
# xlabel('Time, s'), ylabel('Flap, deg'), grid
# title('Flap Test Inputs')

# #       Plot State History
# figure
# subplot(2,2,1)
# plot(t,x(:,1))
# xlabel('Time, s'), ylabel('Axial Velocity (u), m/s'), grid
# title('Forward Body-Axis Component of Inertial Velocity, u')
# subplot(2,2,2)
# plot(t,x(:,2))
# xlabel('Time, s'), ylabel('Side Velocity (v), m/s'), grid
# title('Side Body-Axis Component of Inertial Velocity, v')
# subplot(2,2,3)
# plot(t,x(:,3))
# xlabel('Time, s'), ylabel('Normal Velocity (w), m/s'), grid
# title('Normal Body-Axis Component of Inertial Velocity, z')
# subplot(2,2,4)
# plot(t,x(:,1),t,x(:,2),t,x(:,3))
# xlabel('Time, s'), ylabel('u (blue), v (green), w (red), m/s'), grid
# title('Body-Axis Component of Inertial Velocity')
# legend('Axial velocity, u', 'Side velocity, v', 'Normal velocity, w')

# figure
# subplot(3,2,1)
# plot(t,x(:,4))
# xlabel('Time, s'), ylabel('North (x), m'), grid
# title('North Location, x')

# subplot(3,2,2)
# plot(t,x(:,5))
# xlabel('Time, s'), ylabel('East (y), m'), grid
# title('East Location, y')

# subplot(3,2,3)
# plot(t,-x(:,6))
# xlabel('Time, s'), ylabel('Altitude (-z), m'), grid
# title('Altitude, -z')

# subplot(3,2,4)
# plot((sqrt(x(:,4).*x(:,4) + x(:,5).*x(:,5))),-x(:,6))
# xlabel('Ground Range, m'), ylabel('Altitude, m'), grid
# title('Altitude vs. Ground Range')

# subplot(3,2,5)
# plot(x(:,4),x(:,5))
# xlabel('North, m'), ylabel('East, m'), grid
# title('Ground Track, North vs. East')

# subplot(3,2,6)
# plot3(x(:,4),x(:,5),-x(:,6))
# xlabel('North, m'), ylabel('East, m'), zlabel('Altitude, m'), grid
# title('3D Flight Path')

# figure
# subplot(2,2,1)
# plot(t,x(:,7) * R2D)
# xlabel('Time, s'), ylabel('Roll Rate (p), deg/s'), grid
# title('Body-Axis Roll Component of Inertial Rate, p')
# subplot(2,2,2)
# plot(t,x(:,8) * R2D)
# xlabel('Time, s'), ylabel('m.pitch Rate (q), deg/s'), grid
# title('Body-Axis m.pitch Component of Inertial Rate, q')
# subplot(2,2,3)
# plot(t,x(:,9) * R2D)
# xlabel('Time, s'), ylabel('Yaw Rate (r), deg/s'), grid
# title('Body-Axis Yaw Component of Inertial Rate, r')
# subplot(2,2,4)
# plot(t,x(:,7) * R2D,t,x(:,8) * R2D,t,x(:,9) * R2D)
# xlabel('Time, s'), ylabel('p (blue), q (green), r (red), deg/s'), grid
# title('Body-Axis Inertial Rate Vector Components')
# legend('Roll rate, p', 'm.pitch rate, q', 'Yaw rate, r')

# figure
# subplot(2,2,1)
# plot(t,x(:,10) * R2D)
# xlabel('Time, s'), ylabel('Roll Angle (phi), deg'), grid
# title('Earth-Relative Roll Attitude')
# subplot(2,2,2)
# plot(t,x(:,11) * R2D)
# xlabel('Time, s'), ylabel('m.pitch Angle (theta), deg'), grid
# title('Earth-Relative m.pitch Attitude')
# subplot(2,2,3)
# plot(t,x(:,12) * R2D)
# xlabel('Time, s'), ylabel('Yaw Angle (psi, deg'), grid
# title('Earth-Relative Yaw Attitude')
# subplot(2,2,4)
# plot(t,x(:,10) * R2D,t,x(:,11) * R2D,t,x(:,12) * R2D)
# xlabel('Time, s'), ylabel('phi (blue), theta (green), psi (red), deg'), grid
# title('Euler Angles')
# legend('Roll angle, phi', 'm.pitch angle, theta', 'Yaw angle, psi')
    
VAirRel         =   [] 
vEarth          =   [] 
AlphaAR         =   [] 
BetaAR          =   [] 
windBody        =   [] 
airDensHis      =   [] 
soundSpeedHis   =   [] 
qbarHis         =   [] 
GammaHis        =   [] 
XiHis           =   [] 
    
for i in range(kHis):
    windb           =	WindField(-x[i,5],x[i,9],x[i,10],x[i,11]) 
    windBody        =   np.array([windBody,windb]) 
    airDens,airPres,temp,soundSpeed = Atmos(-x[i,5]) 
    airDensHis      =   np.array([airDensHis,airDens]) 
    soundSpeedHis   =   np.array([soundSpeedHis soundSpeed]) 

    
vBody           =   np.array([x[:,0], x[:,1], x[:,2]).T
vBodyAir        =   vBody + windBody 

for i in range(kHis):
    vE          =   DCM(x[i,9],x[i,10],x[i,11]).T* np.array([vBody[0,i],vBody[1,i],vBody[2,i]]) 
    VER         =   np.sqrt(vE[0]**2 + vEvE[1] **2 + vEvE[2] **2) 
    VAR         =   np.sqrt(vBodyAir[0,i] **2 + vBodyAir[1,i] **2 + vBodyAir[2,i] **2) 
    VARB        =   np.sqrt(vBodyAir[0,i] **2 + vBodyAir[2,i] **2) 
    
    if vBodyAir(1,i) >= 0:
        Alphar      =	m.asin(vBodyAir[2,i] / VARB) 
    else:
        Alphar      =	m.pi - m.asin(vBodyAir[2,i] / VARB) 

    AlphaAR     =   np.array([AlphaAR,Alphar]) 
    Betar       = 	m.asin(vBodyAir[1,i] / VAR) 
    BetaAR      =   np.array([BetaAR,Betar]) 
    vEarth      =   np.array( [vEarth,vE] )
    Xir         =   m.asin(vEarth[1,i] / m.sqrt((vEarth[0,i]) **2 + (vEarth[1,i]) **2)) 
    
    if vEarth[0,i] <= 0    and vEarth[1,i] <= 0:
        Xir     = -m.pi - Xir 

    if vEarth[0,i] <= 0 and vEarth[1,i] >= 0:
        Xir     = m.pi - Xir 

    Gammar      =   m.asin(-vEarth[2,i] / VER) 
    GammaHis    =   [GammaHis Gammar] 
    XiHis       =   [XiHis Xir] 
    VAirRel     =   [VAirRel VAR] 

    
    MachHis         =   VAirRel/soundSpeedHis 
    AlphaDegHis 	=	R2D * AlphaAR 
    BetaDegHis      =	R2D * BetaAR 
    qbarHis         =	0.5 * airDensHis* VAirRel*VAirRel 
    GammaDegHis     =   R2D * GammaHis 
    XiDegHis        =   R2D * XiHis 
    
    # figure
    # subplot(3,1,1)
    # plot(t, VAirRel')
    # xlabel('Time, s'), ylabel('Air-relative Speed, m/s'), grid
    # title('True AirSpeed, Vair')
    # subplot(3,1,2)
    # plot(t, MachHis')
    # xlabel('Time, s'), ylabel('M'), grid
    # title('Mach Number, M')
    # subplot(3,1,3)
    # plot(t, qbarHis')
    # xlabel('Time, s'), ylabel('qbar, N/m **2'), grid
    # title('Dynamic Pressure, qbar')
    
    # figure
    # subplot(2,1,1)
    # plot(t, AlphaDegHis', t, BetaDegHis')
    # xlabel('Time, s'), ylabel('Angle of Attack, deg (blue), Sideslip Angle, deg (green)'), grid
    # title('Aerodynamic Angles'), legend('Angle of Attack, alpha', 'Sideslip Angle, beta')
    # subplot(2,1,2)
    # plot(t, GammaDegHis', t, XiDegHis')
    # xlabel('Time, s'), ylabel('Vertical, deg (blue), Horizontal, deg (green)'), grid
    # title('Flight Path Angles'), legend('Flight Path Angle, gamma', 'Heading Angle, psi')
    
    # 'End of FLIGHT Simulation'