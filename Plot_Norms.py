import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time
import csv
import numpy as np


fig, [ax1,ax2,ax3] = plt.subplots(3)

#NORMS
Iterations = []
L2_rho = []
L2_u = [] 
L2_p= []  

#FLOW Variables
x_cell_center = []
mach_exact = []
mach_CFD = []
rho_exact = []
rho_CFD = []
press_exact = []
press_CFD = []

exact_data = np.loadtxt('Exact_Solution_Isentropic.txt', delimiter=",", dtype='float',usecols = [0,1,2,3,4,5])
x_cell_center = exact_data[:,0]
mach_exact = exact_data[:,2]
rho_exact = exact_data[:,3]
press_exact = exact_data[:,5]


def update(i):
    
    with open('mach.txt', 'r') as mach_CFD:
        mach_CFD = mach_CFD.readlines()[-1]
    
    mach_CFD = np.fromstring(mach_CFD, dtype=float, sep=',')
    mach_CFD = mach_CFD[1:]
    mach_CFD = mach_CFD[:-1]
  
   
    ax1.cla()
    ax1.plot(x_cell_center,mach_CFD,'bo', label = '$Mach_{CFD}$')
    ax1.plot(x_cell_center,mach_exact,'b', label = '$Mach_{Exact}$')
    ax1.legend(loc = 'upper left')


    with open('rho.txt', 'r') as rho_CFD:
        rho_CFD = rho_CFD.readlines()[-1]
    
    rho_CFD = np.fromstring(rho_CFD, dtype=float, sep=',')
    rho_CFD = rho_CFD[1:]
    rho_CFD = rho_CFD[:-1]
  
   
    ax2.cla()
    ax2.plot(x_cell_center,rho_CFD,'bo', label = '$rho_{CFD}$')
    ax2.plot(x_cell_center,rho_exact,'b', label = '$rho_{Exact}$')
    ax2.legend(loc = 'upper right')






    norm_data = np.loadtxt('norm.txt', delimiter=",", dtype='float', usecols = [0,1,2,3])
    Iterations = norm_data[:,0]
    L2_rho = norm_data[:,1]
    L2_u = norm_data[:,2]
    L2_p = norm_data[:,3]
    ax3.cla()
    ax3.semilogy(Iterations,L2_rho, label = '$L_2$: Mass')
    ax3.semilogy(Iterations,L2_u, label = '$L_2$: Momentum')
    ax3.semilogy(Iterations,L2_p, label = '$L_2$: Energy')
    ax3.legend(loc='upper right')
  

    # ax2.xlabel('Iterations')
    # ax2.ylabel('$L_2$')
    

animation = FuncAnimation(fig, update, interval=100)

plt.show()

