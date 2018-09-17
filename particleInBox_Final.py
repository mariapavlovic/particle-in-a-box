import numpy as np
import matplotlib.pyplot as plt

## physical parameters
h = 1 # planck's constant; setting to 1 for simplicity
L = 1 # length of box

k = 1 # wavenumber
E_theoretical = np.pi**2*k**2/2 # energy level (theoretical)
V = 0 # potential energy inside box

## x-axis parameters
x_0 = 0
x_max = L
x_len = 1000
dx = x_max/x_len

E = 3 # approximation of E, to start
dE = 0.01 # amount to increment E for successive approximations
threshold = 0.001 # threshold for difference between psi(L) and 0
E_vals = np.arange(0,7,dE)
m = 1

psi = np.zeros(1000)
psiDer = np.zeros(1000)

psi[0] = 0
psiDer[0] = 1


sums = 0

psi[x_len-1] = 1

xAxis = np.arange(0, L, dx)


while psi[x_len-1] > threshold:
    E += dE
    for i in np.arange(x_len-1):
        psi[i+1] = psi[i] + psiDer[i]*dx
        psiDer[i+1] = psiDer[i] + (-2*m/h**2)*E*psi[i]*dx    
       
for k in psi:
    sums += k**2
    sumdx = sums*dx

psiSum = 0
A = 1/np.sqrt(sumdx)
psi_norm = A*psi

for b in psi:
    psiSum += (b*A)**2
    psiSumdx = psiSum*dx
    

plt.plot(xAxis, psi**2)
plt.title("Probability of a Particle's Position in a Box")
plt.xlabel("dx")
plt.ylabel("Psi**2")
plt.show()

plt.plot(xAxis, psi_norm**2)
plt.title("Probability of a Particle's Position in a Box - Scaled")
plt.xlabel("dx")
plt.ylabel("(Psi*A)**2")
plt.show()
print("probability", psiSumdx)

#print("sum", psi_norm)