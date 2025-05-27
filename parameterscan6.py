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

#Nintervals = np.array([512,800,1000,1500,1700,1800,2000]) # valeurs utilisées pour la conv nsteps
#Nsteps = np.ones(len(Nintervals),dtype = int)*800 # valeurs utilisées pour la conv nsteps

#Nsteps = np.array([800,1000,1200,1500,1700,2000], dtype = int) # valeurs utilisées pour la conv nx
#Nintervals = np.ones(len(Nsteps) , dtype = int)*512 # valeurs utilisées pour la conv nx

#Nsteps = np.array([1000,1200,1500,1700,2000,3000,4000,5000,6000,7000])
#Nintervals = np.ones(len(Nsteps) , dtype = int)*512

Nsteps = np.array([800]) 
Nintervals = np.array([512])

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

t = obs[:,0]
x = pot[:,0]

w0 = 100.0
x0 = -0.5
p0 = 0.0 # changer

def ftPlot() : # j'ajouterai ton animation après parce que ça bugait un peu 

    t = obs[:,0]

    plt.ion() # pour faire l'animation visuelle

    modidx = np.arange(start = 0 , stop = psi2.shape[1], step = 3)
    reeidx = np.arange(start = 1 , stop = psi2.shape[1], step = 3)
    imaidx = np.arange(start = 2 , stop = psi2.shape[1], step = 3)

    for i in range(t.shape[0]) :
    
        mod_a_t = psi2[i,modidx] # module au temps t
        ree_a_t = psi2[i,reeidx] # partie réelle au temps t 
        ima_a_t = psi2[i,imaidx] # partie imaginaire au temps t 
 
        plt.plot(x,mod_a_t , color = "black")
        plt.plot(x,ree_a_t , color = "blue")
        plt.plot(x,ima_a_t , color = "red")
        plt.title(f"t = {t[i]}")
        plt.draw()
        plt.pause(0.005)
        plt.close()
        
    plt.ioff() # pour arrêter


def ObsPlot ( nom_obs , classique = False , obs = obs ) : # plot l'observable correspondante 

    t = obs[:,0] # temps 
    o = np.array([]) # observable
    o2 = np.array([]) # autre observable si besoin d'en afficher deux ( prob gauche droite )
    oo2 = np.array([])
    ylab = "" # ylabel
    lab1 = "" # label pour la courbe si besoin
    couleur = "black"

    if ( nom_obs == "prob_gauche") :
        ylab = "$P_{x<0}(t)$"
        o = obs[:,1]   
    elif ( nom_obs == "prob_droite") :
        ylab = "$P_{x>0}(t)$"
        o = obs[:,2]
    elif ( nom_obs == "prob_d_g") : # plot la prob d'être à gauche et celle d'être à droite sur le même graphe ainsi que leur somme 
        o = obs[:,1] # prob gauche   
        o2 = obs[:,2] # prob droite
        oo2 = o+o2 # pour vérifier que la somme est bien 1
        lab1 = "$P_{x<0}(t)$"
        couleur = "blue"
    elif ( nom_obs == "E") :
        ylab = "$Energie$ [J]"
        o = obs[:,3]        
    elif ( nom_obs == "xmoy") : 
        ylab = "$\\langle x \\rangle (t)$"
        lab1 = "Quantique"
        o = obs[:,4]
    elif ( nom_obs == "x2moy") :
        ylab = "$\\langle x^2 \\rangle (t)$"
        o = obs[:,5]
    elif ( nom_obs == "pmoy") :
        ylab = "$\\langle p \\rangle(t)$"
        lab1 = "Quantique"
        o = obs[:,6]
    elif ( nom_obs == "p2moy") :
        ylab = "$p_{moy}^2(t)$"
        o = obs[:,7]
    elif( nom_obs == "ptot") :# probabilité totale
        ylab = "$P_{totale}(t)$"
        o = obs[:,-1]
    
    else :
        ValueError("Le nom de l'observable doit être : E , xmoy , x2moy , pmoy , p2moy ... ")

    plt.figure()
    plt.plot(t,o,color=couleur, label = lab1)

    if ( nom_obs == "prob_d_g" ) :
        plt.plot(t,o2,color="red", label = "$P_{x>0}(t)$")
        #plt.plot(t,oo2,color="black" , label = "$P_{x>0}(t) + P_{x<0}(t)$")
        #plt.plot(t,obs[:,-1], color = "magenta" , linestyle = "dashed" , label = "$P_{tot}$")
        plt.legend(fontsize = fs - 2)

    Emean = np.mean(obs[:,3]) # pour l'oscillateur classique  

    if (classique and nom_obs == "xmoy") :
        plt.plot( t , x0*np.cos(w0*t) + (p0/w0)*np.sin(w0*t) , color = "red" , label = "Classique" , linestyle = "dashed" ) # 
        plt.legend(fontsize = fs - 2)
    elif (classique and nom_obs == "pmoy") :
        plt.plot( t , - x0*w0*np.sin(w0*t) + p0*np.cos(w0*t) , color = "red" , label = "Classique" , linestyle = "dashed" )
        plt.legend(fontsize = fs - 2)
        
    plt.xlabel("Temps [s]", fontsize = fs)
    plt.ylabel(ylab , fontsize = fs)        

def Vplot(pot) :

    x = pot[:,0]
    V = pot[:,1]
    Emoy = np.mean(obs[:,3])

    plt.figure()
    plt.plot(x,V,color = "black")
    plt.axhline(y=Emoy, color = "magenta" , linestyle = "dashed" , label = "y = $\\langle E \\rangle$")
    plt.xlabel("x [m]", fontsize = fs)
    plt.ylabel("Potentiel [J]", fontsize = fs)
    plt.legend(fontsize = fs - 2)

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
    plt.legend(fontsize = fs - 2)

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

def Convergence ( order = 2 , Nsteps_fixe = False ) : # conv en ordre 2 pour Nsteps fixe et conv en ordre 2 pour Nintervals fixe 

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
            Warning(f"La simulation {outputs[i]} ne se termine pas au temps prévu tfin = {tfin} (voir mesh error), le test de convergence sera pas correct")
                
        xmoye = obser[:,4]
        xfin.append(xmoye[-1])

    plt.figure()
    plt.title(titr)
    plt.xlabel(xlab, fontsize = fs)
    plt.ylabel("$ \\langle x \\rangle (t_{fin})$",fontsize = fs)
    
    if Nsteps_fixe :
        plt.plot(pow(1/Nintervals,order),xfin,"k+-")
    else :
        plt.plot(pow(1/Nsteps,order),xfin,"k+-")


def Ptrans ( trans = 0.035 ) :

    V0 = np.array([100,500,1000,2000,3000], dtype = int)
    param4 = V0
    paramstr4 = "V0"

    nsimul2 = len(V0)

    paramstr3 = "psi_off"
    param3 = 1 
    # Simulations
    outputs2 = []  # List to store output file names
    for i in range(nsimul2):
        output_file = f"{paramstr4}={param4[i]}.out"
        outputs2.append(output_file)
        cmd = f"{executable} {input_filename} {paramstr3}={param3:.15g} {paramstr4}={param4[i]:.15g} output={output_file}"
        print(cmd)
        subprocess.run(cmd, shell=True)
        print('Done.')

    EV0s  = [] # <E>/V0
    Ptrans = []
    
    for i in range(nsimul2) :

        obsV0s = np.loadtxt(outputs2[i]+ "_obs.out")
        Emean = np.mean(obsV0s[:,3])
        EV0s.append((Emean / V0[i]))
        Pdroite = obsV0s[:,2]
        #ObsPlot("prob_d_g" , obs = obsV0s)
        #plt.title(f"$V_0 = {V0[i]}$")
        t = obsV0s[:,0]
        idxtrans = np.argmax(t >= trans)
        Ptrans.append(Pdroite[idxtrans]) # on trouve le premier indice tq t > ttrans

    plt.figure()
    plt.scatter(EV0s,Ptrans, color = "black" , s = 15)
    plt.xlabel("$\\langle E \\rangle / V_0$",fontsize = fs)
    plt.ylabel("$P_{x>0}(t_{trans})$",fontsize = fs)

##def Matteo_Pg_Pd () :     
##
##    plt.figure(figsize=(8, 4))
##    plt.plot(obs[:, 0], obs[:, 1], label="$P_{x<0}(t)$", color="blue")
##    plt.plot(obs[:, 0], obs[:, 2], label="$P_{x>0}(t)$", color="green")
##    plt.axvline(x=0.035, linestyle="dashed", color="red", alpha=0.6, label="$t_{trans}$")
##    plt.xlabel("Temps [s]", fontsize=fs)
##    plt.ylabel("Probabilité", fontsize=fs)
##    #plt.title(f"Probabilités: $V_0 = {V0}, E_i = {obs[0,3].round()}$ ", fontsize=fs)
##    plt.legend()
##    plt.grid(True, linestyle=":")
##    plt.tight_layout()    

#Matteo_Pg_Pd ()
ObsPlot("pmoy" , True)
ObsPlot("xmoy" , True)
ObsPlot("E" , True)
#ObsPlot("ptot" , True)
#ColorPlot(obs,psi2,pot, partie = "module")
#ColorPlot(obs,psi2,pot, partie = "reelle")
#Incertitude(obs) 
Vplot(pot)
#ObsPlot("prob_d_g" , True)
#Convergence( order = 2 , Nsteps_fixe = True )

plt.show()

