#include "ConfigFile.tpp"
#include <chrono>
#include <cmath>
#include <complex> // Pour les nombres complexes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
typedef vector<complex<double>> vec_cmplx;

// Fonction resolvant le systeme d'equations A * solution = rhs
// où A est une matrice tridiagonale
template<class T>
void triangular_solve(vector<T> const& diag,  vector<T> const& lower, vector<T> const& upper,
                 vector<T> const& rhs, vector<T>& solution)
{
    vector<T> new_diag = diag;
    vector<T> new_rhs = rhs;

    // forward elimination
    for (int i(1); i < diag.size(); ++i) {
        T pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i] -= pivot * new_rhs[i - 1];
    }

    solution.resize(diag.size());

    // solve last equation
    solution[diag.size() - 1] = new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    // backward substitution
    for (int i = diag.size() - 2; i >= 0; --i) {
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];
    }
}

/// TODO Potentiel V(x) :
double V(double xa, double xL, double xR, double xb, double m, double w0, double x, double V0, const double PI)
{
    if(x<=xa and x>=xL ){
        return 0.5*m*pow(w0*(x-xa)/(xL-xa),2);
    }else if (x>xa and x<xb){
        return V0*pow(sin((x-xa)*PI/(xb-xa)),2);
    }else {
        return 0.5*m*pow(w0*(x-xb)/(xR-xb),2);
    }
}


// Declaration des diagnostiques de la particule d'apres sa fonction d'onde psi :
//  - prob: calcule la probabilite de trouver la particule dans un intervalle [x_i, x_j]
//  - E:    calcule son energie,
//  - xmoy: calcule sa position moyenne,
//  - x2moy:calcule sa position au carre moyenne,
//  - pmoy: calcule sa quantite de mouvement moyenne,
//  - p2moy:calcule sa quantite de mouvement au carre moyenne.

// Attention ! : contrôler tous les indices dans les sommes 

/// TODO: calculer la probabilite de trouver la particule dans un intervalle [x_i, x_j]
double prob( vec_cmplx const & psi , double dx , size_t i , size_t j  )
{
	double integ = 0. ; 
	
	for ( size_t k(i) ; k < j - 1 ; ++k )
	{ integ += ( norm(psi[k]) + norm(psi[k+1]) ); }	
	integ *= (dx/2.) ; 
	
    return sqrt(integ); // l'intégrale donne la probabilité au carré 
}

/// TODO calculer l'energie
double E( vec_cmplx const & psi, vec_cmplx const &  dH, vec_cmplx const & aH , vec_cmplx const &  cH , double dx, size_t i, size_t j)
{
	complex<double> integ = 0. ; 
	
	for ( size_t k(i) ; k < j -	1 ; ++k )
	{ integ += ( norm(psi[k])*(dH[k]+cH[k]+aH[k]) + norm(psi[k+1])*(dH[k+1]+cH[k+1]+aH[k+1]) )  ; }
	
	return real(integ) ; 
	
}

/// TODO calculer xmoyenne
double xmoy(vec_cmplx const & psi , vector<double> const & x , double dx )
{
	double integ = 0. ; 
	// cout << "psi size : " << (psi.size()-1) << endl ;  ; 
	for ( size_t k(0) ; k < psi.size() - 1 ; ++k )	
	{ integ += ( ( norm(psi[k])*x[k] + norm(psi[k+1])*x[k+1] ) ) ; } // norm retourne le module au carré // cout << "k : " << k << endl ;
	integ *= (dx/2.) ; // on met dx/2 en évidence dans la somme 
	
    return integ;
}

/// TODO calculer x.^2 moyenne
double x2moy(vec_cmplx const & psi , vector<double> const & x , double dx )
{
	double integ = 0. ; 
	for ( size_t k(0) ; k < psi.size() - 1 ; ++k )	
	{ integ += ( ( norm(psi[k])*(pow(x[k],2) + norm(psi[k+1])*pow(x[k+1],2)) ) ) ; } // vérifier 
	integ *= (dx/2.) ; // on met dx/2 en évidence dans la somme 
	
    return integ;
}

/// TODO calculer p moyenne
double pmoy(vec_cmplx const & psi , double dx , double hbar)
{
	complex<double> complex_i = complex<double>(0, 1); // Nombre imaginaire i
	size_t N(psi.size()-1) ; // dernier indice de l'array psi
	
	complex<double> integ = 0. ; 
	for ( size_t k(1) ; k < psi.size() - 2 ; ++k ) // on ne prend pas les valeurs aux bords donc 1 et -3 
	{ integ += ( conj(psi[k]) * (psi[k+1] - psi[k-1]) / (2.*dx) + conj(psi[k+1]) * (psi[k+2] - psi[k]) / (2.*dx) ) ;}
	
	integ += conj(psi[0])*( psi[0] + psi[1] ) / dx +  conj(psi[1])*(psi[0]+psi[2]) / (2.*dx) ; // bord gauche (différences finies à gauche) page 200 
	integ += conj(psi[N])*( psi[N] + psi[N-1]) / dx +  conj(psi[N-1])*( psi[N]+ psi[N-2] ) / (2.*dx) ; // bord droit  (différences finies à droite) page 200
	
	integ *= ( - complex_i * hbar * dx/2. ) ; // on met ihdx/2 en évidence dans la somme (trapèze)
	
    return abs(integ);
}

/// TODO calculer p.^2 moyenne
double p2moy(vec_cmplx const & psi , double dx , double hbar)
{	
	complex<double> integ = 0. ; 
	for ( size_t k(1) ; k < psi.size() - 2 ; ++k ) // on ne prend pas les valeurs aux bords donc 1 et -3 
	{ integ += ( conj(psi[k]) * (psi[k+1] - 2.*psi[k] + psi[k-1]) / pow(dx,2) + conj(psi[k+1]) * (psi[k+2] - 2.*psi[k+1] + psi[k]) / pow(dx,2) ) ;}
	
	// on ne rajoute rien pour les deux bords car 0 ( voir indication ) 
	
	integ *= ( - pow(hbar,2) * dx / 2. ) ; // on met ihdx/2 en évidence dans la somme (trapèze)
	
    return real(integ);
}

/// TODO calculer la normalization
vec_cmplx normalize(vec_cmplx const& psi, double const& dx)
{
    vec_cmplx psi_norm(psi.size());
        
    psi_norm = psi ; 
    
    for ( auto & el : psi_norm ) // on passe par référence 
    { el /= prob(psi,dx,0,psi.size()) ; }
        
    return psi_norm;
}




int
main(int argc, char** argv)
{
    complex<double> complex_i = complex<double>(0, 1); // Nombre imaginaire i
    const double PI = 3.1415926535897932384626433832795028841971e0;

    string inputPath("configuration.in.example"); // Fichier d'input par defaut
    if (argc > 1) // Fichier d'input specifie par l'utilisateur ("./Exercice8 config_perso.in")
        inputPath = argv[1];

    ConfigFile configFile(
      inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for (int i(2); i < argc;
         ++i) // Input complementaires ("./Exercice8 config_perso.in input_scan=[valeur]")
        configFile.process(argv[i]);

    // Set verbosity level. Set to 0 to reduce printouts in console.
    const int verbose = configFile.get<int>("verbose");
    configFile.setVerbosity(verbose);

    // Parametres physiques :
    double hbar = 1.;
    double m = 1.;
    double tfin = configFile.get<double>("tfin");
    double xL = configFile.get<double>("xL");
    double xR = configFile.get<double>("xR");
    double xa = configFile.get<double>("xa");
    double xb = configFile.get<double>("xb");
    double V0 = configFile.get<double>("V0");
    double om0 = configFile.get<double>("om0");
    double n  = configFile.get<int>("n"); // Read mode number as integer, convert to double
    bool   psi_off = configFile.get<bool>("psi_off") ; 

    double x0 = configFile.get<double>("x0");
    double sigma0 = configFile.get<double>("sigma_norm") * (xR - xL);

    int Nsteps = configFile.get<int>("Nsteps");
    int Nintervals = configFile.get<int>("Nintervals");


	double L = xR - xL ; // ajouté 
    /// TODO: initialiser le paquet d'onde, equation (4.116) du cours
    double k0 = 2. * PI * n / L ;

    int Npoints = Nintervals + 1;
    double dx = (xR - xL) / Nintervals;
    double dt = tfin / Nsteps;

    const auto simulationStart = std::chrono::steady_clock::now();

    // Maillage :
    vector<double> x(Npoints);
    for (int i(0); i < Npoints; ++i)
        x[i] = xL + i * dx;

    // Initialisation de la fonction d'onde :
    vec_cmplx psi(Npoints);

    // initialization time and position to check Probability
    double t = 0;
    unsigned int Nx0 = floor((0 - xL)/(xR-xL)*Npoints); //chosen xR*0.5 since top of potential is at half x domain
  
    /// TODO initialize psi
    for (int i(0); i < Npoints; ++i)
    { psi[i] = exp( complex_i * k0 * x[i] ) * exp( - 0.5 * pow((x[i] - x0),2) / pow(sigma0,2) ) ; }
       
    // Modifications des valeurs aux bords :
    psi[0] = complex<double>(0., 0.);
    psi[Npoints - 1] = complex<double>(0., 0.);
    
    // Normalisation :
    psi = normalize(psi, dx);

    // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
    vec_cmplx dH(Npoints), aH(Nintervals), cH(Nintervals); // matrice Hamiltonienne
    vec_cmplx dA(Npoints), aA(Nintervals),
      cA(Nintervals); // matrice du membre de gauche de l'equation (4.100)
    vec_cmplx dB(Npoints), aB(Nintervals),
      cB(Nintervals); // matrice du membre de droite de l'equation (4.100)

    complex<double> a =
      complex_i * hbar * dt / (4.*m*dx*dx); // Coefficient complexe a de l'equation (4.100)

    /// TODO: calculer les éléments des matrices A, B et H.
    // Ces matrices sont stockées sous forme tridiagonale, d:diagonale, c et a: diagonales
    // supérieures et inférieures
    for (int i(0); i < Npoints; ++i) // Boucle sur les points de maillage
    {
		complex<double> bi = complex_i * dt * V(xa,xL,xR,xb,m,om0,x[i],V0,PI) / (2.*hbar) ; 
		
		if ( i == 0 or i == Npoints-1){ dH[i] = 1. ; }
        else { dH[i] = 2.*pow(hbar,2) / (2.*m*pow(dx,2)) + V(xa,xL,xR,xb,m,om0,x[i],V0,PI) ; }
        dA[i] = 1. + a + bi ; 
        dB[i] = 1. - a - bi ; 
        
    }
    for (int i(0); i < Nintervals; ++i) // Boucle sur les intervalles
    {
        aH[i] =  - 1.*pow(hbar,2) / (2.*m*pow(dx,2)) ;
        aA[i] = -a ; 
        aB[i] =  a ; 
        cH[i] =  - 1.*pow(hbar,2) / (2.*m*pow(dx,2)) ;
        cA[i] = -a ;  
        cB[i] =  a ;   
    }

    // Conditions aux limites: psi nulle aux deux bords
    /// TODO: Modifier les matrices A et B pour satisfaire les conditions aux limites

	// voir slides 34 Week 10 
	dH[0]   = 1. ; dA[0] = 1. ; dB[0] = 1. ;
	dH[-1]  = 1. ; dA[-1] = 1. ; dB[-1] = 1. ;
	aH[0]  = 0. ;  aA[0] = 0.  ; aB[0] = 0.;
	aH[-1] = 0. ;  aA[-1] = 0. ; aB[-1] = 0.;
	cH[0]  = 0. ;  cA[0] = 0.  ; cB[0] = 0.;
	cH[-1] = 0. ;  cA[-1] = 0. ; cB[-1] = 0.;



    // Fichiers de sortie :
    string output = configFile.get<string>("output");

    ofstream fichier_potentiel((output + "_pot.out").c_str());
    fichier_potentiel.precision(15);
    for (int i(0); i < Npoints; ++i)
        fichier_potentiel << x[i] << " " << V(xa,xL,xR,xb,m,om0,x[i],V0,PI) << endl;
    fichier_potentiel.close();

    ofstream fichier_psi((output + "_psi2.out").c_str());
    fichier_psi.precision(6);

    ofstream fichier_observables((output + "_obs.out").c_str());
    fichier_observables.precision(15);

    // t0 writing
    for (int i(0); i < Npoints; ++i){
        fichier_psi << abs(psi[i])  << " " << real(psi[i]) << " "  << imag(psi[i]) << " ";
        }
    fichier_psi << endl;

    // Ecriture des observables :
    /// TODO: introduire les arguments des fonctions prob, E, xmoy, x2moy, pmoy et p2moy
    ///       en accord avec la façon dont vous les aurez programmés plus haut
    fichier_observables << t << " " << prob(psi,dx,0,psi.size()-1) << " " << prob(psi,dx,0,psi.size()-1) // attention : contrôler tous les indices ( xa et 0 pour le premier prob et o et xb pour le deuxiüme prob ) 
                << " " << E(psi,dH,aH,cH,dx,0,x.size()-1) << " " << xmoy (psi,x,dx) << " "  
                << x2moy(psi,x,dx) << " " << pmoy (psi,dx,hbar) << " " << p2moy(psi,dx,hbar) << endl; 

    // Boucle temporelle :    
    while (t < tfin) {

        // Multiplication psi_tmp = B * psi :
        vec_cmplx psi_tmp(Npoints, 0.);
        for (int i(0); i < Npoints; ++i)
            psi_tmp[i] = dB[i] * psi[i];
        for (int i(0); i < Nintervals; ++i) {
            psi_tmp[i] += cB[i] * psi[i + 1];
            psi_tmp[i + 1] += aB[i] * psi[i];
        }

        // Resolution de A * psi = psi_tmp :
        triangular_solve(dA, aA, cA, psi_tmp, psi);
        t += dt;

        // t0 writing
        if ( not psi_off )  // pour le test de convergence on écrit pas psi car trop volumineux. 
        {
        for (int i(0); i < Npoints; ++i){
            fichier_psi << abs(psi[i])  << " " << real(psi[i]) << " "  << imag(psi[i]) << " ";
            }
        }
        fichier_psi << endl;

        // Ecriture des observables :
	/// TODO: introduire les arguments des fonctions prob, E, xmoy, x2moy, pmoy et p2moy
	///       en accord avec la façon dont vous les aurez programmés plus haut
        fichier_observables << t << " " << prob(psi,dx,0,psi.size()) << " " << prob(psi,dx,0,psi.size()) // Attenzione ! Controllare gli indici ! 
                    << " " << E(psi,dH,aH,cH,dx,0,x.size()-1) << " " << xmoy (psi,x,dx) << " "  
                    << x2moy(psi,x,dx) << " " << pmoy (psi,dx,hbar) << " " << p2moy(psi,dx,hbar) << endl; 

    } // Fin de la boucle temporelle





    fichier_observables.close();
    fichier_psi.close();

    const auto simulationEnd = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsedSeconds = simulationEnd - simulationStart;
    std::cout << "Simulation finished in " << setprecision(3) << elapsedSeconds.count()
              << " seconds" << std::endl;
}
