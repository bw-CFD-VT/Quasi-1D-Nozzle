import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time
import csv
import numpy as np


fig, [ax1,ax2,ax3,ax4] = plt.subplots(4)

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
    
    with open('mach.txt', 'r') as mach_CFD:
        mach_CFD = mach_CFD.readlines()[-1]
    
    mach_CFD = np.fromstring(mach_CFD, dtype=float, sep=',')
    mach_CFD = mach_CFD[1:]
    mach_CFD = mach_CFD[:-1]

    with open('u.txt', 'r') as u_CFD:
        u_CFD = u_CFD.readlines()[-1]
    
    u_CFD = np.fromstring(u_CFD, dtype=float, sep=',')
    u_CFD = u_CFD[1:]
    u_CFD = u_CFD[:-1]
  
    
    ax1.cla()
    ax1.scatter(x_cell_center,u_CFD, s = 20, facecolors ='none', edgecolors ='b', label = '$u_{CFD}$')
    ax1.plot(x_cell_center,u_exact,'b', label = '$u_{Exact}$')
    ax1.legend(loc = 'upper left')

    with open('press.txt', 'r') as press_CFD:
       press_CFD = press_CFD.readlines()[-1]
    
    press_CFD = np.fromstring(press_CFD, dtype=float, sep=',')
    press_CFD = press_CFD[1:]
    press_CFD = press_CFD[:-1]
  
   
    ax3.cla()
    ax3.scatter(x_cell_center,press_CFD, s = 20, facecolors ='none', edgecolors ='b', label = '$Pressure_{CFD}$')
    ax3.plot(x_cell_center,press_exact,'b', label = '$Pressure_{Exact}$')
    ax3.legend(loc = 'upper right')


    with open('rho.txt', 'r') as rho_CFD:
        rho_CFD = rho_CFD.readlines()[-1]
    
    rho_CFD = np.fromstring(rho_CFD, dtype=float, sep=',')
    rho_CFD = rho_CFD[1:]
    rho_CFD = rho_CFD[:-1]
  
   
    ax2.cla()
    ax2.scatter(x_cell_center,rho_CFD, s = 20, facecolors ='none', edgecolors ='b', label = '$rho_{CFD}$')
    ax2.plot(x_cell_center,rho_exact,'b', label = '$rho_{Exact}$')
    ax2.legend(loc = 'upper right')






    norm_data = np.loadtxt('norm.txt', delimiter=",", dtype='float', usecols = [0,1,2,3])
    Iterations = norm_data[:,0]
    L2_rho = norm_data[:,1]
    L2_u = norm_data[:,2]
    L2_p = norm_data[:,3]
    ax4.cla()
    ax4.semilogy(Iterations,L2_rho, label = '$L_2$: Mass')
    ax4.semilogy(Iterations,L2_u, label = '$L_2$: Momentum')
    ax4.semilogy(Iterations,L2_p, label = '$L_2$: Energy')
    ax4.legend(loc='upper right')
  

    # ax2.xlabel('Iterations')
    # ax2.ylabel('$L_2$')
    

animation = FuncAnimation(fig, update, interval=100)

plt.show()

