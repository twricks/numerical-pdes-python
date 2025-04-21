#Finite difference method code to find numerical solution of 1-D heat equation
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.pyplot import cm
from pylab import *
import array as arr
import statistics as stat
import PyQt6
np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)
import random
matplotlib.use('qtagg')
matplotlib.pyplot.autoscale(enable=False, axis='both',tight=None)
smlltxt = 8
mdtxt = 10
lrgtext = 12
plt.rc('font', size=smlltxt)
# controls default text sizes
plt.rc('axes', titlesize=smlltxt)
# fontsize of the axes title
plt.rc('axes', labelsize=mdtxt)
# fontsize of the x and y labels
plt.rc('xtick', labelsize=smlltxt)
# fontsize of the tick labels
plt.rc('ytick', labelsize=smlltxt)
# fontsize of the tick labels
plt.rc('legend', fontsize=smlltxt)# legend fontsize
plt.rc('figure', titlesize=lrgtext)# fontsize of the figure title
#Finite difference method code to find numerical solution of 1-D heat equation
#spatial interval [0,1] divided into N+1 equally spaced sample points xn=n*xh
#time interval [0,T] divided into M+1, tk = k*th

M=200 #sets the precision in time; number of rows will be M+1
N = 9 #sets the precision in space; number of columns willbe N+1
T=1 #length of the time interval. Set to 1.
S=1 #length of the space interval
hx = round(S/N,3)
ht = round(T/M,3)
timestep = "{:.2e}".format(ht)
#pi=3.1415
alph=0.4 #diffusivity constant, set to 1; decreasing it arbitarily only results
#in greater time to solution, as shown in respective contour plots, but any
#increase results in errors.
BC0=70 #temperatures maintained at the spatial boundaries.
BC1=40
A = BC1 - BC0
#time, space, temperature meshes
t=np.zeros((M+1),dtype=float) #number of points = divisor +0
x=np.zeros((N+1),dtype=float)
ux=np.zeros((N+1),dtype=float) #list used to print all
#values in the spatial interval for each time
um=np.zeros((M+1),dtype=float)
u=np.zeros((len(t),len(x)),dtype=float) #time in rows,
#columns in space, a value at every point
#Fixing boundary temperatures
for i in arange(0,len(t)):
    u[i][0]=BC0
    ux[0]=BC0
    u[i][N]=BC1
    ux[N]=BC1
#IC: stepping through all non-boundary points at the initial time, setting f(x)
for j in arange(1,len(x)-1):
    u[0][j]= random.randint(0,15)
#f(x); initial distribution of temperatures
ux[j]=random.randint(0,100) #for display

for m in arange(0,len(t)): #loading a time mesh
    t[m]=m*ht

for q in arange(0,len(x)): #loading a space mesh
    x[q]=q*hx
#next to establish boundary & initial conditions; the boundary conditions correspond to
#temperatures maintained at both ends, u(0,t)=BC0,
#u(1,t)=BC1. the initial condition u(x,0)=f(x) corresponds
#to some initial distribution of heat values in space, f(x),
#for example, sin(x).

#Finite Difference algorithm
for k in arange(0,len(t)-1): #rows
    for n in arange(1,len(x)-1): #columns
        u[k+1][n]= u[k][n]+alph*(ht/(pow(hx,2)))*(u[k][n+1]-
        2*(u[k][n])+u[k][n-1])
        ux[n]=u[k][n] ##used for display
    um[k]=np.mean(ux)
    um[len(t)-1]=np.mean(ux) #mean of the last row.
#print("t: ",k)
print(ux)
    
#Plotting average over time
fig = plt.figure()
ax = fig.add_subplot(111)
plt.suptitle('Heat Average vs. Time')
plt.title('Time step: ' + str(timestep))
plt.plot(t,um,marker='', color='green', linestyle='dashed',
linewidth=2, label = 'Average over time')
plt.legend(loc='upper left')
plt.xlabel('Time')
plt.ylabel ('Average Temperature')
#xlim(0-1,T+1)
ylim(0,100)
#plt.show()

err=0
finerr=0
#Error in steady-state solution
for e in arange(0,len(ux)):
    err = ux[e]-(A*x[e]+BC0)
    finerr = err + finerr
    err = 0
FinalError = "{:.2e}".format(finerr)
#Plotting end state
fig = plt.figure()
ax = fig.add_subplot(111)
plt.suptitle('Final Heat Distribution')
plt.title('Time step: ' + str(timestep)+', Error:' +
str(FinalError) + '.')
plt.plot(x,ux,marker='', color='green', linestyle='dashed',
linewidth=2,label='Final Temperature')
plt.legend(loc='lower left')
plt.xlabel('Length')
plt.ylabel ('Temperature')
#xlim(0-1,S+1)
ylim(0,100) #e.g, Degrees C
plt.show()