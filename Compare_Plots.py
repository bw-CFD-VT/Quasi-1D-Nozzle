# AOE 6145
# Homework 2: Quasi-1D Nozzle FVM Code
# Brendan Walsh (PID: bwalsh4)

#To be ran after obtaining solution for all grids -> compare results
# import matplotlib.gridspec as gridspec
# from matplotlib.ticker import ScalarFormatter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

p_0 = 300000.0
T_0 = 600.0

# grid_ID = "_1.txt"
color1,color2,color3,color4,color5,color6 = 'b','g','r','m','c','y'   #Grid 1-6

#Mach # Result All Grids
x_1 = x_2 = x_3 = x_4 = x_5 = x_6 = []
mach_CFD_1 = mach_CFD_2 = mach_CFD_3 = mach_CFD_4 = mach_CFD_5 = mach_CFD_6 = []
#Exact Solution Data -> Extract location of cell center for specific grid and corresponding Mach # from CFD Solution
def x_mach_cell_center(filename,grid_ID):
    x_cell_center = result_CFD = []
    exact_data = np.loadtxt('Exact_Solution_Isentropic'+grid_ID, delimiter=",", dtype='float',usecols = [0])
    x_cell_center = exact_data

    with open(filename+grid_ID, 'r') as result_CFD:
        result_CFD = result_CFD.readlines()[-1]
    result_CFD = np.fromstring(result_CFD, dtype=float, sep=',')
    result_CFD = result_CFD[1:]
    result_CFD = result_CFD[:-1]

    return x_cell_center,result_CFD

x_1,mach_CFD_1 = x_mach_cell_center('mach','_1.txt')
x_2,mach_CFD_2 = x_mach_cell_center('mach','_2.txt')
x_3,mach_CFD_3 = x_mach_cell_center('mach','_3.txt')
x_4,mach_CFD_4 = x_mach_cell_center('mach','_4.txt')
x_5,mach_CFD_5 = x_mach_cell_center('mach','_5.txt')
x_6,mach_CFD_6 = x_mach_cell_center('mach','_6.txt')

plt.figure(1)
plt.plot(x_1,mach_CFD_1,color1, label = 'h = 1')
plt.plot(x_2,mach_CFD_2,color2, label = 'h = 2')
plt.plot(x_3,mach_CFD_3,color3, label = 'h = 4')
plt.plot(x_4,mach_CFD_4,color4, label = 'h = 8')
plt.plot(x_5,mach_CFD_5,color5, label = 'h = 16')
plt.plot(x_6,mach_CFD_6,color6, label = 'h = 32')
plt.xlabel('x (m)')
plt.ylabel('Mach')
plt.legend(loc = 'upper left')


# Observed Order + Global DE error

def Norm_Error(filename,grid_ID):

    with open(filename+grid_ID, 'r') as Error_Norm:
        Error_Norm = Error_Norm.readlines()[-1]
    Error_Norm = np.fromstring(Error_Norm, dtype=float, sep=',')

    return Error_Norm

error_norm_U_1 = Norm_Error('error_norm_U','_1.txt')
error_norm_U_2 = Norm_Error('error_norm_U','_2.txt')
error_norm_U_3 = Norm_Error('error_norm_U','_3.txt')
error_norm_U_4 = Norm_Error('error_norm_U','_4.txt')
error_norm_U_5 = Norm_Error('error_norm_U','_5.txt')
error_norm_U_6 = Norm_Error('error_norm_U','_6.txt')



r = 2.0
h_global = [1,2,3,4,5,6]
h_order = [1,2,3,4,5]
order_second = [1e-4,0.0006125,0.0023029,0.0065678,0.0154316,0.0312788]

Global_DE_rho = [error_norm_U_1[0],error_norm_U_2[0],error_norm_U_3[0],error_norm_U_4[0],error_norm_U_5[0],error_norm_U_6[0]]
Global_DE_u = [error_norm_U_1[1],error_norm_U_2[1],error_norm_U_3[1],error_norm_U_4[1],error_norm_U_5[1],error_norm_U_6[1]]
Global_DE_p = [error_norm_U_1[2],error_norm_U_2[2],error_norm_U_3[2],error_norm_U_4[2],error_norm_U_5[2],error_norm_U_6[2]]

p_hat_rho = [np.log(error_norm_U_2[0]/error_norm_U_1[0])/np.log(r),np.log(error_norm_U_3[0]/error_norm_U_2[0])/np.log(r),
            np.log(error_norm_U_4[0]/error_norm_U_3[0])/np.log(r),np.log(error_norm_U_5[0]/error_norm_U_4[0])/np.log(r),
            np.log(error_norm_U_6[0]/error_norm_U_5[0])/np.log(r)]
p_hat_u = [np.log(error_norm_U_2[1]/error_norm_U_1[1])/np.log(r),np.log(error_norm_U_3[1]/error_norm_U_2[1])/np.log(r),
            np.log(error_norm_U_4[1]/error_norm_U_3[1])/np.log(r),np.log(error_norm_U_5[1]/error_norm_U_4[1])/np.log(r),
            np.log(error_norm_U_6[1]/error_norm_U_5[1])/np.log(r)]
p_hat_p = [np.log(error_norm_U_2[2]/error_norm_U_1[2])/np.log(r),np.log(error_norm_U_3[2]/error_norm_U_2[2])/np.log(r),
            np.log(error_norm_U_4[2]/error_norm_U_3[2])/np.log(r),np.log(error_norm_U_5[2]/error_norm_U_4[2])/np.log(r),
            np.log(error_norm_U_6[2]/error_norm_U_5[2])/np.log(r)]                        

plt.figure(2)
plt.loglog(h_global,Global_DE_rho,linestyle='-', marker='o',label = 'Mass')
plt.loglog(h_global,Global_DE_u,linestyle='-', marker='o',label = 'Momentum')
plt.loglog(h_global,Global_DE_p,linestyle='-', marker='o',label = 'Energy')
plt.loglog(h_global,order_second,'k',label='$2^{nd}$ Order Slope')
plt.xlabel('h')
plt.ylabel('Global DE $L_{2}$ Norm')
plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f'))
plt.gca().xaxis.set_minor_formatter(mtick.FormatStrFormatter('%.0f'))
plt.ylim(10e-10, 10e6)
plt.legend(loc = 'upper left')

plt.figure(3)
plt.semilogx(h_order,p_hat_rho,linestyle='-', marker='o',label = 'Mass')
plt.semilogx(h_order,p_hat_u,linestyle='-', marker='o',label = 'Momentum')
plt.semilogx(h_order,p_hat_p,linestyle='-', marker='o',label = 'Energy')
plt.xlabel('h')
plt.ylabel('Order of Accuracy, p')
plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f'))
plt.gca().xaxis.set_minor_formatter(mtick.FormatStrFormatter('%.0f'))
plt.ylim(0, 3)
plt.legend(loc = 'lower right')





def x_const_comparison(filename,grid_ID):
    iteration = CFD_soln = []
    CFD_soln = np.loadtxt(filename+grid_ID, delimiter=",", dtype='float',usecols = [0,50])
    iteration = CFD_soln[:,0]
    result_CFD = CFD_soln[:,1]

    return iteration,result_CFD

plt.figure(4)
x_3_1,press_3_1 = x_const_comparison('press','_3.txt')
x_3_2,press_3_2 = x_const_comparison('press','_3_2.txt')
x_3_3,press_3_3 = x_const_comparison('press','_3_3.txt')
x_3_4,press_3_4 = x_const_comparison('press','_3_4.txt')
x_3_5,press_3_5 = x_const_comparison('press','_3_5.txt')
x_3_6,press_3_6 = x_const_comparison('press','_3_6.txt')
plt.plot(x_3_1,press_3_1/p_0,color1, label = 'Cells = 160, CFL = 0.125, $K_{2}$ = 0.00, $K_{4}$ = 1/32')
plt.plot(x_3_2,press_3_2/p_0,color2, label = 'Cells = 160, CFL = 0.125, $K_{2}$ = 0.25, $K_{4}$ = 1/32')
plt.plot(x_3_3,press_3_3/p_0,color3, label = 'Cells = 160, CFL = 0.125, $K_{2}$ = 0.50, $K_{4}$ = 1/32')
plt.plot(x_3_4,press_3_4/p_0,color4, label = 'Cells = 160, CFL = 0.125, $K_{2}$ = 0.50, $K_{4}$ = 1/64')
plt.plot(x_3_5,press_3_5/p_0,color5, label = 'Cells = 160, CFL = 0.125, $K_{2}$ = 0.25, $K_{4}$ = 1/64')
plt.plot(x_3_6,press_3_6/p_0,color6, label = 'Cells = 160, CFL = 0.125, $K_{2}$ = 0.00, $K_{4}$ = 1/64')

plt.xlabel('Iterations')
plt.ylabel('p/$p_{0}$')
plt.legend(loc = 'lower right')

error_norm_U_3_1 = Norm_Error('error_norm_U','_3.txt')
error_norm_U_3_2 = Norm_Error('error_norm_U','_3_2.txt')
error_norm_U_3_3 = Norm_Error('error_norm_U','_3_3.txt')
error_norm_U_3_4 = Norm_Error('error_norm_U','_3_4.txt')
error_norm_U_3_5 = Norm_Error('error_norm_U','_3_5.txt')
error_norm_U_3_6 = Norm_Error('error_norm_U','_3_6.txt')

h_global = [1,2,3,4,5,6]

Global_DE_rho_const = [error_norm_U_3_1[0],error_norm_U_3_2[0],error_norm_U_3_3[0],error_norm_U_3_4[0],error_norm_U_3_5[0],error_norm_U_3_6[0]]
Global_DE_u_const = [error_norm_U_3_1[1],error_norm_U_3_2[1],error_norm_U_3_3[1],error_norm_U_3_4[1],error_norm_U_3_5[1],error_norm_U_3_6[1]]
Global_DE_p_const = [error_norm_U_3_1[2],error_norm_U_3_2[2],error_norm_U_3_3[2],error_norm_U_3_4[2],error_norm_U_3_5[2],error_norm_U_3_6[2]]
plt.figure(5)
plt.semilogy(h_global,Global_DE_rho_const,linestyle='-', marker='o',label = 'Mass')
plt.semilogy(h_global,Global_DE_u_const,linestyle='-', marker='o',label = 'Momentum')
plt.semilogy(h_global,Global_DE_p_const,linestyle='-', marker='o',label = 'Energy')
plt.xlabel('h')
plt.ylabel('Global DE $L_{2}$ Norm')
plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f'))
plt.gca().xaxis.set_minor_formatter(mtick.FormatStrFormatter('%.0f'))
plt.ylim(10e-10, 10e6)
plt.legend(loc = 'upper left')



plt.show()
