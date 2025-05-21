import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import os

# Parameters
# TODO adapt to what you need (folder path executable input filename)
executable = r"/Users/a-x-3/Desktop/Exercice6_2025_student/exe"  # Name of the executable (NB: .exe extension is required on Windows)
repertoire = r"/Users/a-x-3/Desktop/Exercice6_2025_student"
os.chdir(repertoire)

input_filename = 'configuration.in.example'  # Name of the input file

# ------------------------------------- Simulations ----------------------------------- #

#nsteps = np.array([2,3,4,5,6,7,8,9])*200 # valeurs utilisées pour la convergence
#nx = np.array([2,3,4,5,6,7,8,9])*20 # valeurs utilisées pour la convergence

#nsteps = np.array([50e3]) # valeurs utilisées pour le tsunami
#nx = np.array([10000]) # valeurs utilisées pour le tsunami 

Nsteps = np.array([800]) 
Nintervals = np.array([512])

paramstr = 'Nsteps'  # Parameter name to scan
param = Nsteps  # Parameter values to scan

paramstr2 = 'Nintervals'  # Parameter name to scan
param2 = Nintervals  # Parameter values to scan

nsimul = len(Nsteps)

# Simulations
outputs = []  # List to store output file names
convergence_list = []
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}_{paramstr2}={param2[i]}.out"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    cmd = f"{executable} {input_filename} {paramstr}={param[i]:.15g} {paramstr2}={param2[i]:.15g} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

lw = 1.5
fs = 16

psi2 = np.loadtxt(outputs[-1]+"_psi2.out")
obs = np.loadtxt(outputs[-1]+ "_obs.out")
pot = np.loadtxt(outputs[-1]+ "_pot.out")

t = obs[0,:]
x = pot[0,:]


def xmoy_t (obs) :

    t = obs[:,0]
    xmoy = obs[:,4]

    #for i in xmoy :
        #print(i)
    #for i in t :
        #print(i)

    #plt.scatter(1,1)
    plt.figure()
    plt.plot(t,xmoy)
    plt.xlabel("Temps [s]", fontsize = fs)
    plt.ylabel("$x_{moy}(t)$", fontsize = fs)

def pmoy_t (obs) :

    t = obs[:,0]
    pmoy = obs[:,6]

    plt.figure()
    plt.plot(t,pmoy)
    plt.xlabel("Temps [s]", fontsize = fs)
    plt.ylabel("$p_{moy}(t)$", fontsize = fs)

def Vplot(pot) :

    x = pot[:,0]
    V = pot[:,1]

    plt.figure()
    plt.plot(x,V)
    plt.xlabel("x []", fontsize = fs)
    plt.ylabel("V", fontsize = fs)

def ColorPlot (obs,psi,pot) : 

    t = obs[:,0]
    x = pot[:,0]
    abs_ = psi[:,0]

    plt.pcolor(x,t,abs_)
    

    print(psi.shape)

xmoy_t(obs)
pmoy_t(obs)
plt.show()

