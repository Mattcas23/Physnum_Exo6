import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import os

# Parameters
# TODO adapt to what you need (folder path executable input filename)
executable = r"/Users/Matte/Desktop/EPFL/BA4/Physnum/Physnum_Exo4/exe" #r"/Users/a-x-3/Desktop/Ex3_2024_student/exe"  # Name of the executable (NB: .exe extension is required on Windows)
repertoire = r"/Users/Matte/Desktop/EPFL/BA4/Physnum/Physnum_Exo4"#r"/Users/a-x-3/Desktop/Ex3_2024_student"
os.chdir(repertoire)

input_filename = 'configuration.in.example'  # Name of the input file

#----------------------------------------- Valeurs du Configfile --------------------------------- # 

Values = np.genfromtxt("/Users/Matte/Desktop/EPFL/BA4/Physnum/Physnum_Exo4/configuration.in.example")

#R = Values[0,-1]
#r1 = Values[1,-1]
#epsilon_a = Values[2,-1]
#epsilon_b= Values[3,-1]
#uniform_rho_case = Values[4,-1]
#VR = Values[5,-1]
#rho0  =  Values[6,-1]
N1    = Values[9,-1]
#N2  = Values[8,-1]
#maxit  = Values[9,-1]
#output = Values[10,-1]
#sampling = Values[11,-1]
#verbose = Values[12,-1]
# ---------------------------------------------------------------
lw = 1.5
fs = 16
N1_list = [16, 32, 64, 128, 256]  # exemple
phi_center_list = []
config_path = "/Users/Matte/Desktop/EPFL/BA4/Physnum/Physnum_Exo4/configuration.in.example"
phi_output_path = "/Users/Matte/Desktop/EPFL/BA4/Physnum/Physnum_Exo4/phioutput.out"
# ------------------------------------------- Simulations ---------------------------------------------

data = np.loadtxt("/Users/Matte/Desktop/EPFL/BA4/Physnum/Physnum_Exo4/phioutput.out")
data2 = np.loadtxt("/Users/Matte/Desktop/EPFL/BA4/Physnum/Physnum_Exo4/Eoutput.out") # Load the output file of the i-th simulation
r = data[:, 0]
phi = data[:,1] # nombre de pas de temps total pour la simulation
rE= data2[:,0]
E = data2[:,1]
print(phi,N1)
def conv(): #phi(r=0) en fonction du nombre de points de maillage N1(=N2)
    for N1 in N1_list:
        # Charger et modifier le fichier de configuration
        config_data = np.genfromtxt(config_path)
        config_data[10, -1] = N1  # Modifier la ligne correspondant à N1
        np.savetxt(config_path, config_data, fmt="%.8e")  # Sauvegarder

        # Lancer l'exécutable
        subprocess.run([executable])

        # Lire la sortie (phi au centre)
        data = np.loadtxt(phi_output_path)
        phi = data[:, 1]  # colonne du potentiel
        phi_center = phi[0]  # ou un autre critère : moyenne, etc.
        phi_center_list.append(phi_center)
        print(phi)

    # Tracer la convergence
    plt.figure()
    plt.plot(N1_list, phi_center_list, 'ko-', linewidth=lw)
    plt.xlabel("Nombre de points N1", fontsize=fs)
    plt.ylabel(r"$\phi(r=0)$ (V)", fontsize=fs)
    plt.title("Convergence de $\phi(r=0)$ en fonction de N1", fontsize=fs)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def pot():
    plt.plot(r,phi)
    #plt.plot(rE,E)
    plt.xlabel("r,rE", fontsize=fs)
    plt.ylabel("phi, E", fontsize=fs)
    plt.grid()
pot()
conv() #marche pas 
plt.show()
