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

Nintervals = np.array([512,800,1000,1500,1700,1800,2000]) # valeurs utilisées pour la conv nsteps
Nsteps = np.ones(len(Nintervals),dtype = int)*800 # valeurs utilisées pour la conv nsteps

#Nsteps = np.array([800,1000,1200,1500,1700,2000], dtype = int) # valeurs utilisées pour la conv nx
#Nintervals = np.ones(len(Nsteps) , dtype = int)*512 # valeurs utilisées pour la conv nx

#Nsteps = np.array([800]) 
#Nintervals = np.array([512])

#Nsteps = np.array([3000])
#Nintervals = np.array([512])

paramstr = 'Nsteps'  # Parameter name to scan
param = Nsteps  # Parameter values to scan

paramstr2 = 'Nintervals'  # Parameter name to scan
param2 = Nintervals  # Parameter values to scan



paramstr3 = "psi_off"
param3 = 0 # false => on écrit le fichier psi2 
if len(Nintervals) > 1 : # pour le test de convergence, on écrit pas les valeurs de psi car trop volumineux
    param3 = 1 # true => on écrit pas dans le fichier psi2 
    

nsimul = len(Nsteps)

# Simulations
outputs = []  # List to store output file names
convergence_list = []
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}_{paramstr2}={param2[i]}.out"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    cmd = f"{executable} {input_filename} {paramstr}={param[i]:.15g} {paramstr2}={param2[i]:.15g} {paramstr3}={param3:.15g} output={output_file}"
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
x0 = -0.5
p0 = 1.0 # changer

def ObsPlot ( nom_obs , classique = False ) : # plot l'observable correspondante 

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
        ylab = "$E$ [J]"
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

    Emean = np.mean(obs[:,3]) # pour l'oscillateur classique  

    if (classique and nom_obs == "xmoy") :
        plt.plot( t , x0*np.cos(w0*t) +  (p0/w0)*np.sin(w0*t) , color = "red" , label = "Classique" , linestyle = "dashed" )
        plt.legend()
    elif (classique and nom_obs == "pmoy") :
        plt.plot( t , - x0*w0*np.sin(w0*t) + p0*np.cos(w0*t) , color = "red" , label = "Classique" , linestyle = "dashed" )
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

    if param3 > 0 :
        ValueError("Nintervals.size > 1 : on écrit pas phi dans le fichier")

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

def Convergence ( order = 2 , Nsteps_fixe = False ) : # conv en ordre 2 pour Nsteps fixe et conv en ordre pour Nintervals fixe 

    xfin = []

    xlab = f"$(\\Delta t)^{order}$" # Label x
    titr = f"$n_x = {Nintervals[0]}$" # titre du graphe 

    if Nsteps_fixe :
        
        xlab = f"$(\\Delta x)^{order}$" # Label x
        titr = "$n_{steps} = $" + f"{Nsteps[0]}" # titre du graphe 

    for i in range(nsimul) :
        
        obser = np.loadtxt(outputs[i]+ "_obs.out")
        tobs = obser[:,0]
        tfin = tobs[-1]

        if ( tfin < 0.08 ) :
            ValueError(f"La simulation {outputs[i]} ne se termine pas au temps prévu tfin = {tfin} (voir mesh error), le test de convergence sera pas correct")
                
        xmoye = obser[:,4]
        xfin.append(xmoye[-1])

    plt.figure()
    plt.title(titr)
    plt.xlabel(xlab, fontsize = fs)
    plt.ylabel("$x_{moy}(t)$",fontsize = fs)
    
    if Nsteps_fixe :
        plt.plot(pow(1/Nintervals,order),xfin,"k+-")
    else :
        plt.plot(pow(1/Nsteps,order),xfin,"k+-")

    

ObsPlot("E")    
Convergence(2,True)
ObsPlot("pmoy", True)
ObsPlot("xmoy", True)
#Incertitude(obs)
#Vplot(pot)
#ColorPlot(obs,psi2,pot)
#ObsPlot("xmoy",True)
plt.show()

