import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time
import csv
import numpy as np


# fig, [ax1,ax2,ax3,ax4,ax5] = plt.subplots(5)
fig,ax = plt.subplots(nrows = 4,ncols = 2)

#NORMS
Iterations = []
L2_rho = []
L2_u = [] 
L2_p= []  

#FLOW Variables
x_cell_center = []
area = []
mach_exact = []
mach_CFD = []
rho_exact = []
rho_CFD = []
press_exact = []
press_CFD = []
u_exact = []
u_CFD = []
mass_flow = []

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
  
    ax[0][0].cla()
    ax[0][0].plot(x_cell_center,rho_exact,'b', label = '$rho_{Exact}$')
    ax[0][0].plot(x_cell_center,rho_CFD,'k--', label = '$rho_{CFD}$')
    ax[0][0].legend(loc = 'upper right')
    #---------------------------------------------------------------------#

    #-------------------- CFD Velocity (x-comp.) -------------------------#
    with open('u.txt', 'r') as u_CFD:
        u_CFD = u_CFD.readlines()[-1]
    
    u_CFD = np.fromstring(u_CFD, dtype=float, sep=',')
    u_CFD = u_CFD[1:]
    u_CFD = u_CFD[:-1]
    
    ax[1][0].cla()
    ax[1][0].plot(x_cell_center,u_exact,'b', label = '$u_{Exact}$')
    ax[1][0].plot(x_cell_center,u_CFD, 'k--', label = '$u_{CFD}$')
    ax[1][0].legend(loc = 'upper left')
    #---------------------------------------------------------------------#

    #-------------------- CFD Pressure -----------------------------------#
    with open('press.txt', 'r') as press_CFD:
       press_CFD = press_CFD.readlines()[-1]
    
    press_CFD = np.fromstring(press_CFD, dtype=float, sep=',')
    press_CFD = press_CFD[1:]
    press_CFD = press_CFD[:-1]
  
   
    ax[2][0].cla()
    ax[2][0].plot(x_cell_center,press_exact,'b', label = '$Pressure_{Exact}$')
    ax[2][0].plot(x_cell_center,press_CFD, 'k--', label = '$Pressure_{CFD}$')
    ax[2][0].legend(loc = 'upper right')
    #---------------------------------------------------------------------#

    #-------------------- CFD Mach ---------------------------------------#
    with open('mach.txt', 'r') as mach_CFD:
        mach_CFD = mach_CFD.readlines()[-1]
    
    mach_CFD = np.fromstring(mach_CFD, dtype=float, sep=',')
    mach_CFD = mach_CFD[1:]
    mach_CFD = mach_CFD[:-1]

    ax[3][0].cla()
    ax[3][0].plot(x_cell_center,mach_exact,'b', label = '$mach_{Exact}$')
    ax[3][0].plot(x_cell_center,mach_CFD, 'k--', label = '$mach_{CFD}$')
    ax[3][0].legend(loc = 'upper left')
    ax[3][0].get_shared_x_axes().join(ax[0], ax[1],ax[2],ax[3])
    #---------------------------------------------------------------------#

    #-------------------- Residual Norms ---------------------------------#
    norm_data = np.loadtxt('residual_norm.txt', delimiter=",", dtype='float', usecols = [0,1,2,3])
    Iterations = norm_data[:,0]
    L2_rho = norm_data[:,1]
    L2_u = norm_data[:,2]
    L2_p = norm_data[:,3]
    ax[4][1].cla()
    ax[4][1].set(xlabel = 'Iterations',ylabel = '$L_2$ Norm')
    ax[4][1].semilogy(Iterations,L2_rho, label = '$L_2$: Mass')
    ax[4][1].semilogy(Iterations,L2_u, label = '$L_2$: Momentum')
    ax[4][1].semilogy(Iterations,L2_p, label = '$L_2$: Energy')
    ax[4][1].legend(loc='upper right')
    #---------------------------------------------------------------------#

animation = FuncAnimation(fig, update, interval=100)

plt.show()

