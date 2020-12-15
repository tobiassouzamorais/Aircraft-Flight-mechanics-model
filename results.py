# Results
import matplotlib.pyplot as plt
import simulation
import numpy as np
time, x_t, u_t, Aero_par = simulation.simul()
v_0=str('Body-axis x inertial velocity, ub, m/s')
v_1=str('Body-axis y inertial velocity, vb, m/s')
v_2=str('Body-axis z inertial velocity, wb, m/s')
v_3=str('North position of center of mass WRT Earth, xe, m')
v_4=str('East position of center of mass WRT Earth, ye, m')
v_5=str('Negative of c.m. altitude WRT Earth, ze = -h, m')
v_6=str('Body-axis roll rate, pr, rad/s')
v_7=str('Body-axis pitch rate, qr, rad/s')
v_8=str('Body-axis yaw rate, rr,rad/s')
v_9=str('Roll angle of body WRT Earth, phir, rad')
v_10=str('Pitch angle of body WRT Earth, thetar, rad')
v_11=str('Yaw angle of body WRT Earth, psir, rad')
v=np.array([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8, v_9, v_10, v_11])
for i in range(0,12):
    plt.plot(time[:], x_t[i,:])
    plt.ylabel(v[i])
    plt.xlabel('time [s]')
    plt.show()
