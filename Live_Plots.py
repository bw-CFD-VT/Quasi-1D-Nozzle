# AOE 6145
# Homework 2: Quasi-1D Nozzle FVM Code
# Brendan Walsh (PID: bwalsh4)

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec
import numpy as np


# fig, [ax1,ax2,ax3,ax4,ax5] = plt.subplots(5)
fig = plt.figure(constrained_layout=True)
spec2 = gridspec.GridSpec(ncols=4, nrows=2, figure=fig)
ax1 = fig.add_subplot(spec2[1, 0])
ax2 = fig.add_subplot(spec2[1, 1])
ax3 = fig.add_subplot(spec2[1, 2])
ax4 = fig.add_subplot(spec2[1, 3])
ax5 = fig.add_subplot(spec2[0, :])

#NORMS
Iterations = []
L2_rho = []
L2_u = []
L2_p= []  

#Flow Variables
x_cell_center = []
area = []
mach_exact = []
rho_exact = []
press_exact = []
u_exact = []
mach_CFD = []
rho_CFD = []
press_CFD = []
u_CFD = []

#Exact Solution Data
exact_data = np.loadtxt('Exact_Solution_Isentropic.txt', delimiter=",", dtype='float',usecols = [0,1,2,3,4,5])
x_cell_center = exact_data[:,0]
area = exact_data [:,1]
mach_exact = exact_data[:,2]
rho_exact = exact_data[:,3]
u_exact = exact_data[:,4]
press_exact = exact_data[:,5]


def update(i):
    
    #-------------------- CFD Density ------------------------------------#
    with open('rho.txt', 'r') as rho_CFD:
        rho_CFD = rho_CFD.readlines()[-1]
    
    rho_CFD = np.fromstring(rho_CFD, dtype=float, sep=',')
    rho_CFD = rho_CFD[1:]
    rho_CFD = rho_CFD[:-1]
  
    ax1.cla()
    ax1.set(xlabel = 'x (m)',ylabel = '\u03C1 (kg/$m^3$)')
    ax1.plot(x_cell_center,rho_exact,'b', label = '$rho_{Exact}$')
    ax1.plot(x_cell_center,rho_CFD,'k--', label = '$rho_{CFD}$')
    ax1.legend(loc = 'upper right')
    #---------------------------------------------------------------------#

    #-------------------- CFD Velocity (x-comp.) -------------------------#
    with open('u.txt', 'r') as u_CFD:
        u_CFD = u_CFD.readlines()[-1]
    
    u_CFD = np.fromstring(u_CFD, dtype=float, sep=',')
    u_CFD = u_CFD[1:]
    u_CFD = u_CFD[:-1]
    
    ax2.cla()
    ax2.set(xlabel = 'x (m)',ylabel = 'u (m/s)')
    ax2.plot(x_cell_center,u_exact,'b', label = '$u_{Exact}$')
    ax2.plot(x_cell_center,u_CFD, 'k--', label = '$u_{CFD}$')
    ax2.legend(loc = 'upper left')
    #---------------------------------------------------------------------#

    #-------------------- CFD Pressure -----------------------------------#
    with open('press.txt', 'r') as press_CFD:
       press_CFD = press_CFD.readlines()[-1]
    
    press_CFD = np.fromstring(press_CFD, dtype=float, sep=',')
    press_CFD = press_CFD[1:]
    press_CFD = press_CFD[:-1]
  
   
    ax3.cla()
    ax3.set(xlabel = 'x (m)',ylabel = 'p (Pa)')
    ax3.plot(x_cell_center,press_exact,'b', label = '$Pressure_{Exact}$')
    ax3.plot(x_cell_center,press_CFD, 'k--', label = '$Pressure_{CFD}$')
    ax3.legend(loc = 'upper right')
    #---------------------------------------------------------------------#

    #-------------------- CFD Mach ---------------------------------------#
    with open('mach.txt', 'r') as mach_CFD:
        mach_CFD = mach_CFD.readlines()[-1]
    
    mach_CFD = np.fromstring(mach_CFD, dtype=float, sep=',')
    mach_CFD = mach_CFD[1:]
    mach_CFD = mach_CFD[:-1]

    ax4.cla()
    ax4.set(xlabel = 'x (m)',ylabel = 'Mach')
    ax4.plot(x_cell_center,mach_exact,'b', label = '$mach_{Exact}$')
    ax4.plot(x_cell_center,mach_CFD, 'k--', label = '$mach_{CFD}$')
    ax4.legend(loc = 'upper left')
    #---------------------------------------------------------------------#

    #-------------------- Residual Norms ---------------------------------#
    norm_data = np.loadtxt('residual_norm.txt', delimiter=",", dtype='float', usecols = [0,1,2,3])
    Iterations = norm_data[:,0]
    L2_rho = norm_data[:,1]
    L2_u = norm_data[:,2]
    L2_p = norm_data[:,3]

    ax5.cla()
    ax5.set(xlabel = 'Iterations',ylabel = '$L_2$ Norm')
    ax5.semilogy(Iterations,L2_rho, label = '$L_2$: Mass')
    ax5.semilogy(Iterations,L2_u, label = '$L_2$: Momentum')
    ax5.semilogy(Iterations,L2_p, label = '$L_2$: Energy')
    ax5.legend(loc='upper right')
    #---------------------------------------------------------------------#

animation = FuncAnimation(fig, update, interval=100)

plt.show()

