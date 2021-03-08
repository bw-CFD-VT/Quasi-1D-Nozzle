import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time
import csv
import numpy as np

x = []
L2_rho = []
L2_u = [] 
L2_p= []  


def update(i):
    data = np.loadtxt('norm.txt', delimiter=",", dtype='float',usecols = [0,1,2,3])

    x = data[:,0]
    L2_rho = data[:,1]
    L2_u = data[:,2]
    L2_p = data[:,3]
    plt.cla()
    plt.plot(x,L2_rho, label = 'L2_rho')
    plt.plot(x,L2_u, label = 'L2_u')
    plt.plot(x,L2_p, label = 'L2_p')

    plt.legend(loc='upper left')
    plt.tight_layout()

    plt.gca().relim()
    plt.gca().autoscale_view()
animation = FuncAnimation(plt.gcf(), update, interval=10)


plt.tight_layout()
plt.show()