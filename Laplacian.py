import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import *
from matplotlib.pyplot import cm
from pylab import *
import array as arr
import PyQt6
np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)
#finite difference method to find numerical solution to Laplace Equation in a
#two-dimensional region
#with the boundaries held at specified potentials
matplotlib.use('qtagg')
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
plt.rc('legend', fontsize=smlltxt)
# legend fontsize
plt.rc('figure', titlesize=lrgtext) # fontsize of the figure title

NR=30 #Number of rows, number of columns.
NC=30
V=np.zeros((NR,NC),dtype=float)
#print ("enter the value of the potential at each boundary:")
#vt = float(input("Enter potential on the top boundary:"))
vt=1
#vb = float(input("Enter potential on the bottom boundary:"))
vb =1
#vl = float(input("Enter potential on the left boundary: :"))
vl=0.3
#vr = float(input("Enter potential on the right boundary:"))
vr = 0.7

#print(vt,vb,vl,vr)
#print ("Enter the value of the relaxation parameter: ")
#w = float(input('Relaxation parameter '))
w=1.5
rlxprm = "{:.2e}".format(w)
#print ("Enter the value of the iteration parameter: ")
#dm = float(input('Iteration parameter '))
dm = 0.000001
#print (w)
#setting the initial values. we do this by stepping through every grid point
#along every boundary and loading the boundary condition into an array
for i in arange(0,NC):
    V[0][i]= vt
    V[NC-1][i]= vb
#print(V)
for j in arange(0,NR):
    V[j][0]=vl
    V[j][NR-1]=vr
#print(V)
#prints correctly, with [0][0] beginning in the top left on
#the screen
#iteration procedure; Vnext =
#(w/4)[vjplus+viplus+vjminusnext+viminusnext]
#+(1-w)V
dw=1
while dw > dm:
    dw = 0
    for I in range(1,NR-1): #rows
        for J in range(1,NC-1): #columns
            VIJ = w * ((V[I][J+1]+V[I+1][J]+V[I][J-1]+V[I-1][J])/4)+(1-w)*V[I][J]
            dc = abs(VIJ-V[I][J])
            if (dc > dw):
                dw=dc
            V[I][J] = VIJ
#print(V)

#print("range(NC)",range(NC))
#print("range(NR)",range(NR))
##Now to plot the data
x = np.linspace(0, NC, NC)
y = np.linspace(0, NR, NR)
X,Y= meshgrid(x,y)
fig = plt.figure()
ax = plt.axes(projection = '3d')
ax.plot_wireframe(X,Y,V, color = 'darkgreen',
label='Potential')
plt.suptitle('Electrostatic Potential Distribution')
plt.title('Relaxation Parameter: ' + str(rlxprm) + '.')
ax.set_ylabel('Y')
ax.set_xlabel('X')
ax.set_zlabel('V')
plt.legend(loc='lower left')
ax.view_init(32, -123)
plt.show() #print the wireframe plot