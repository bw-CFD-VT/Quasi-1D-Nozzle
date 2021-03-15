# AOE 6145
# Homework 2: Quasi-1D Nozzle FVM Code
# Brendan Walsh (PID: bwalsh4)

#To be ran in conjunction with 1D code to monitor flow variables and residuals live

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec
import numpy as np

Case_Flag = 1
grid_ID = 'Check.txt' #Update to match Grid ID text in main c++ file
# color = 'b'   #Grid 1
# color = 'g'   #Grid 2 
color = 'r'   #Grid 3
# color = 'm'   #Grid 4
# color = 'c'   #Grid 5
# color = 'y'   #Grid 6

fig = plt.figure(constrained_layout=True)
spec = gridspec.GridSpec(ncols=5, nrows=2, figure=fig)
ax1 = fig.add_subplot(spec[0, 0:1])
ax2 = fig.add_subplot(spec[0, 1:2])
ax3 = fig.add_subplot(spec[1, 0:1])
ax4 = fig.add_subplot(spec[1, 1:2])
ax5 = fig.add_subplot(spec[0, 2:5])
ax6 = fig.add_subplot(spec[1, 2:5])

#NORMS
Iterations = L2_rho = L2_u = L2_p = []

#Flow Variables
x_cell_center = area = []
mach_exact = rho_exact = u_exact = press_exact = temp_exact = []
mach_CFD = rho_CFD = u_CFD = press_CFD = temp_CFD = []


#Exact Solution Data
exact_data = np.loadtxt('Exact_Solution_Isentropic'+grid_ID, delimiter=",", dtype='float',usecols = [0,1,2,3,4,5,6])
x_cell_center = exact_data[:,0]
area = exact_data [:,1]
mach_exact = exact_data[:,2]
rho_exact = exact_data[:,3]
u_exact = exact_data[:,4]
press_exact = exact_data[:,5]
temp_exact = exact_data[:,6]

#Function to extract current iteration's CFD solution
def read_data(filename,grid_ID):
    with open(filename+grid_ID, 'r') as result_CFD:
        result_CFD = result_CFD.readlines()[-1]
    
    result_CFD = np.fromstring(result_CFD, dtype=float, sep=',')
    result_CFD = result_CFD[1:]
    result_CFD = result_CFD[:-1]
    
    return result_CFD


def update(i):
    
    #-------------------- CFD Density ------------------------------------#
    rho_CFD = read_data('rho',grid_ID) 
  
    ax1.cla()
    ax1.set(xlabel = 'x (m)',ylabel = '\u03C1 (kg/$m^3$)')
    if Case_Flag == 1:
        ax1.plot(x_cell_center,rho_exact,'k', label = 'Exact')
    ax1.plot(x_cell_center,rho_CFD, color + '--', label = 'CFD')
    ax1.legend(loc = 'upper right')
    #---------------------------------------------------------------------#

    #-------------------- CFD Velocity (x-comp.) -------------------------#
    u_CFD = read_data('u',grid_ID)

    ax2.cla()
    ax2.set(xlabel = 'x (m)',ylabel = 'u (m/s)')
    if Case_Flag == 1:
        ax2.plot(x_cell_center,u_exact,'k', label = 'Exact')
    ax2.plot(x_cell_center,u_CFD, color + '--', label = 'CFD')
    ax2.legend(loc = 'upper left')
    #---------------------------------------------------------------------#

    #-------------------- CFD Pressure -----------------------------------#
    press_CFD = read_data('press',grid_ID)

    ax3.cla()
    ax3.set(xlabel = 'x (m)',ylabel = 'p (Pa)')
    if Case_Flag == 1:
        ax3.plot(x_cell_center,press_exact,'k', label = 'Exact')
    ax3.plot(x_cell_center,press_CFD, color + '--', label = 'CFD')
    ax3.legend(loc = 'upper right')
    #---------------------------------------------------------------------#

    #-------------------- CFD Temperature --------------------------------#
    temp_CFD = read_data('temp',grid_ID)

    ax4.cla()
    ax4.set(xlabel = 'x (m)',ylabel = 'T (K)')
    if Case_Flag == 1:
        ax4.plot(x_cell_center,temp_exact,'k', label = 'Exact')
        ax4.legend(loc = 'upper right')
    ax4.plot(x_cell_center,temp_CFD, color + '--', label = 'CFD')
    if Case_Flag == 2:
        ax4.legend(loc = 'lower left')
    #---------------------------------------------------------------------#

    #-------------------- CFD Mach ---------------------------------------#
    mach_CFD = read_data('mach',grid_ID)

    ax5.cla()
    ax5.set(xlabel = 'x (m)',ylabel = 'Mach')
    if Case_Flag == 1:
        ax5.plot(x_cell_center,mach_exact,'k', label = 'Exact')
    ax5.plot(x_cell_center,mach_CFD, color + '--', label = 'CFD')
    ax5.legend(loc = 'upper left')
    #---------------------------------------------------------------------#

    #-------------------- Residual Norms ---------------------------------#
    norm_data = np.loadtxt('residual_norm'+grid_ID, delimiter=",", dtype='float', usecols = [0,1,2,3])
    Iterations = norm_data[:,0]
    L2_rho = norm_data[:,1]
    L2_u = norm_data[:,2]
    L2_p = norm_data[:,3]

    ax6.cla()
    ax6.set(xlabel = 'Iterations',ylabel = '$L_2$ Norm')
    ax6.semilogy(Iterations,L2_rho, label = '$L_2$: Mass')
    ax6.semilogy(Iterations,L2_u, label = '$L_2$: Momentum')
    ax6.semilogy(Iterations,L2_p, label = '$L_2$: Energy')
    ax6.legend(loc='upper right')
    #---------------------------------------------------------------------#

animation = FuncAnimation(fig, update, interval=100)

plt.show()

