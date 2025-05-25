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

def ObsPlot ( nom_obs ) :

    t = obs[:,0] # temps 
    o = np.array([]) # observable 
    ylab = "" # ylabel

    if ( nom_obs == "prob_gauche") :
        ylab = "$P_{gauche}$"
        o = obs[:,1]   
    elif ( nom_obs == "prob_droite") :
        ylab = "$P_{droite}$"
        o = obs[:,2]       
    elif ( nom_obs == "E") :
        ylab = "$E$"
        o = obs[:,3]        
    elif ( nom_obs == "xmoy") : 
        ylab = "$x_{moy}(t)$"
        o = obs[:,4]
    elif ( nom_obs == "x2moy") :
        ylab = "$x_{moy}^2(t)$"
        o = obs[:,5]
    elif ( nom_obs == "pmoy") :
        ylab = "$p_{moy}(t)$"
        o = obs[:,6]
    elif ( nom_obs == "p2moy") :
        ylab = "$p_{moy}^2(t)$"
        o = obs[:,7]
    else :
        ValueError("Le nom de l'observable doit être : E , xmoy , x2moy , pmoy , p2moy ")

    plt.figure()
    plt.plot(t,o,color="black")
    plt.xlabel("Temps [s]", fontsize = fs)
    plt.ylabel(ylab , fontsize = fs)        

def Vplot(pot) :

    x = pot[:,0]
    V = pot[:,1]

    plt.figure()
    plt.plot(x,V,color = "black")
    plt.xlabel("x []", fontsize = fs)
    plt.ylabel("V", fontsize = fs)


def ColorPlot (obs,psi,pot) : 

    t = obs[:,0]
    x = pot[:,0]

    absidx = np.arange(start = 0 , stop = psi.shape[1], step = 3)

    psi_abs = psi[:,absidx]
    


    plt.figure()
    plt.pcolor(x,t,psi_abs)
    plt.xlabel("x [m]", fontsize = fs)
    plt.ylabel("Temps [s]", fontsize = fs)
    cbar = plt.colorbar()
    cbar.set_label("$|\\psi(x,t)|$", fontsize = fs)

    

ObsPlot("p2moy")
ObsPlot("pmoy")
Vplot(pot)
ColorPlot(obs,psi2,pot)
plt.show()

