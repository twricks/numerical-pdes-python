import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.pyplot import cm
from pylab import *
import array as arr
import PyQt6
matplotlib.use('qtagg')
np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)
from numpy import linalg as LA
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
#Code using the Finite-difference method to find a numerical solution of the 
# infinite square well
#problem in quantum mechanics
#BC: Because V(x) =infinity outside the well, the wave function u(x) must 
# vanish at the boundaries, V(x)=0 abs(x) <a.
c=0
a=1 #spatial interval

N= 200 #sets precision in space
h= a/N #step size
Es = 6 #how many energy levels we wish to calculate
x=np.zeros((N+1),dtype=float) #spatial mesh= number of columns in u = N +1 to account for 0
disp=np.zeros((N+1),dtype=float)
anlyt=np.zeros((N+1),dtype=float)
v=np.zeros((N+1),dtype=float) #potential array
H=np.zeros((N+1,N+1),dtype=float) #Hu = Eu energy eigenvalue problem
u=np.zeros((Es,N+1),dtype=float) #Wavefunction mesh
#loading spatial mesh
for q in arange(0,len(x)): #loading a space mesh
    x[q]=q*h
#v[25]=600
#v[30]=600
#loading n x n H matrix, for which the eigenvalues will be found;
H[0][0]= 2/ pow(h,2)+v[0]
H[0][1]= -1/pow(h,2)

while (c < N-3):
    for i in arange(0,len(x)-2):
        H[i+1][i]= -1/pow(h,2)
        H[i+1][i+1]= 2/ pow(h,2)+v[i]
        H[i+1][i+2]= -1/pow(h,2)
    c = c+1
H[N][N-1]=-1/pow(h,2)
H[N][N]= 2/ pow(h,2)+v[N]
#print(H)
E,vecs=LA.eig(H)
#print(E) #energy eigenvalues
#print(vecs) #normalized 1-D eigenvectors, each column 
#(the way the linalg package works) corresponding to each energy

for m in arange(0,Es): #loading wave mesh, stepping through whole space (u[C] for however many energies u[R])
    for n in arange(0,len(x)):
        u[m][n]=vecs[n][m]**2 
##phi squared is the probability of finding the particle at that point, and the
#integral of phi squared over the spatial interval is always 1.

e=int(input("Enter Energy Level E=1,2,3..."))
for g in arange(0,len(x)):
    disp[g]=u[e-1][g] 
#row changed to view each energy level
    anlyt[g]=sin(x[g]*pi*e)**2
#normalizing the analytical solution for plotting
intgrlphi2 = sum(anlyt)
normsol = anlyt/ intgrlphi2
#print(sum(normsol))
#print(sum(disp)) ##checking to be sure that wavefunction is normalized, i.e. that integral of phi^2 is one.
stepsize = "{:.2e}".format(h)
plt.suptitle('Eigensolutions to TISE (E = ' + str(e) +').')
plt.title('Space step: ' + stepsize + '.')
plt.plot(x,disp,marker='', color='green',
linestyle='dashed', linewidth=2,label='Numerical Solution')
plt.plot(x,normsol,marker='', color='maroon',
linewidth=2,label='Analytical Solution')
plt.legend(loc='lower left')
plt.xlabel('X')
plt.ylabel ('$ψ^2$')
plt.show()