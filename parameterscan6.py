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

w0 = 100.0

def ObsPlot ( nom_obs , analytique = False ) : # plot l'observable correspondante 

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
    plt.plot(t,o,color="black", label = "Quantique")

    if (analytique and nom_obs == "xmoy") :
        plt.plot( t ,  np.cos(w0*t) +  np.sin(w0*t) , color = "red" , label = "Classique" , linestyle = "dashed" )
        plt.legend()
    elif (analytique and nom_obs == "pmoy") :
        plt.plot( t , w0*np.sin(w0*t) - w0*np.cos(w0*t) , color = "red" , label = "Classique" , linestyle = "dashed" )
        plt.legend()
        
    plt.xlabel("Temps [s]", fontsize = fs)
    plt.ylabel(ylab , fontsize = fs)        

def Vplot(pot) :

    x = pot[:,0]
    V = pot[:,1]

    plt.figure()
    plt.plot(x,V,color = "black")
    plt.xlabel("x [m]", fontsize = fs)
    plt.ylabel("V", fontsize = fs)

def Incertitude (obs) : # pour vérifier le principe d'incertitude d'Heisenberg

    t = obs[:,0]
    Dx = np.sqrt(obs[:,5] - pow(obs[:,4],2)) # Delta x 
    Dp = np.sqrt(obs[:,7] - pow(obs[:,6],2)) # Delta p 
    DxDp = Dx*Dp 

    plt.figure()
    plt.plot(t,DxDp, color = "black")
    plt.axhline(y = 0.5, color = "red", label = "$y = \\frac{\\hbar}{2}$") # ici on a pris des unités normalisées donc hbar = 1 (voir exercice p1)
    plt.xlabel("Temps [s]", fontsize = fs)
    plt.ylabel("$\\langle \\Delta x \\rangle \\langle \\Delta p \\rangle $", fontsize = fs)
    plt.legend()

def ColorPlot (obs,psi,pot, partie = "module" ) : # partie = "module" / "reelle" / "imaginaire"

    t = obs[:,0]
    x = pot[:,0]

    idx = np.array([])
    leg = ""

    if ( partie == "module" ) : 
        idx = np.arange(start = 0 , stop = psi.shape[1], step = 3) # on choisit uniquement les indexes correspondant au module (voir c++)
        leg = "$|\\psi(x,t)|$" # légende correspondante
    elif ( partie == "reelle" ) :
        idx = np.arange(start = 1 , stop = psi.shape[1], step = 3) # on choisit uniquement les indexes correspondant à la partie réelle (voir c++)
        leg = "$Re[\\psi(x,t)]$" # légende correspondante
    elif ( partie == "imaginaire" ) :
        idx = np.arange(start = 2 , stop = psi.shape[1], step = 3) # on choisit uniquement les indexes correspondant à la partie réelle (voir c++)
        leg = "$Im[\\psi(x,t)]$" # légende correspondante
    else :
        ValueError("partie = 'module' / 'reelle' / 'imaginaire'")

    psi_val = psi[:,idx]
    
    plt.figure()
    plt.pcolor(x,t,psi_val)
    plt.xlabel("x [m]", fontsize = fs)
    plt.ylabel("Temps [s]", fontsize = fs)
    cbar = plt.colorbar()
    cbar.set_label(leg, fontsize = fs)
    

ObsPlot("pmoy", True)
Incertitude(obs)
#Vplot(pot)
ColorPlot(obs,psi2,pot,"reelle")
ObsPlot("xmoy",True)
plt.show()

