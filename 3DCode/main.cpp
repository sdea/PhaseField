/*
    Armadillo  project
    Based on the linear algebra C++ Library Armadillo:
    Armadillo (http://arma.sourceforge.net/)

    Created by Salvatore De Angelis on 1/9/17.
    Copyright © 2017 Salvatore De Angelis. All rights reserved.
*/

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;
typedef double T;
typedef long int  lint;


/// Function to perform derivatives
/// Laplacian operator
cube d2x_3D(cube const& C, T const& h) {


    lint end = C.n_rows -1;

    // Computation using compact stencil
    cube out1 = C(span(3,end-1),span(2,end-2),span(2,end-2)) + C(span(1,end-3),span(2,end-2),span(2,end-2));
    cube out2 = C(span(2,end-2),span(3,end-1),span(2,end-2)) + C(span(2,end-2),span(1,end-3),span(2,end-2));
    cube out3 = C(span(2,end-2),span(2,end-2),span(3,end-1)) + C(span(2,end-2),span(2,end-2),span(1,end-3));
    cube outCenter = 6*C(span(2,end-2),span(2,end-2),span(2,end-2));


    return (out1 + out2 + out3 - outCenter)/(h*h);

}

/// Function for gradient in 2D
cube gradx(cube const& C, T const& h) {

    // !!!! Only works for square matrix !!!
    lint end = C.n_rows -1;
    cube Cx = C(span(3,end-1), span(1, end-3), span(1, end-3)) - C(span(1,end-3),span(1,end-3), span(1, end-3));

    return Cx/(2*h);
}

cube grady(cube const& C, T const& h) {

    // !!!! Only works for square matrix !!!
    lint end = C.n_rows -1;
    cube Cy = C(span(1, end-3),span(3,end-1), span(1, end-3)) - C(span(1,end-3),span(1,end-3), span(1, end-3));

    return Cy/(2*h);
}

cube gradz(cube const& C, T const& h) {

    // !!!! Only works for square matrix !!!
    lint end = C.n_rows -1;
    cube Cz = C(span(1, end-3), span(1, end-3), span(3,end-1)) - C(span(1,end-3), span(1, end-3),span(1,end-3));

    return Cz/(2*h);
}


/// Covert radiant to deg
T convToDeg(T rad_angle) {

    return ((T)180/datum::pi)*rad_angle;

}

/// Print header function
void printHeader(lint& N, T& delta, std::string& def_str) {


    cout << "\n#######################################################\n" << std::endl;
    cout << "P H A S E  F I E L D  C A L C U L A T I O N\n" << std::endl;
    cout << "Contact angle calculation based on the SBM method\n" << std::endl;
    cout << "Code written by: Salvatore De Angelis & Martina Trini\n" << std::endl;

    cout << "PARAMETERS:  \nDomain size:  " << N << "  Delta:  " << delta << std::endl;
    cout << "Writing output in:  " << def_str << std::endl;
}


int main(int argc, const char * argv[]) {


    // Parameters of simulation
    if (argc != 2) {

        cout <<"Wrong number of parameters!!!" << std::endl;
        cout <<"Required: (T) delta" << std::endl;

        return -1;

    }

    cube geom;
    std::string nameGeom = "geom.bin";

    geom.load(nameGeom,arma::raw_binary);
    geom.reshape(101,101,101);
    //geom.save("geom.mat", arma_ascii,101);

    cube Cini;
    std::string nameCini = "Cini.bin";

    Cini.load(nameCini,arma::raw_binary);
    Cini.reshape(101,101,101);

    cube Psi;
    std::string namePsi = "Psi.bin";

    Psi.load(namePsi,arma::raw_binary);
    Psi.reshape(101,101,101);


    T h = 0.2;
    T dt = 1e-6;
    lint Time = 1;

    T W = 1.;
    lint bcType = 1;
    T eps = 1.;


    // Set last number of arrays
    lint end = geom.n_rows -1;
    lint dend = geom.n_rows -5;


    // Contact angle values
    T theta = ((T)100/(T)180)*datum::pi;
    T cosTheta = cos(theta);

    // Print contact angles
    cout << "\nTheta:  " << theta << "  Cos(theta):  \n" << cosTheta <<std::endl;


    // Compute the gradient of Psi and the term dPsi/Psi
    cube invPsi = 1./(Psi+1);
    cube dxPsi = gradx(Psi,h);
    cube dyPsi = grady(Psi,h);
    cube dzPsi = gradz(Psi,h);
    cube modGradPsi = sqrt(dxPsi%dxPsi + dyPsi%dyPsi + dzPsi%dzPsi);
    cube dxPsiOvPsi = invPsi(span(2,end-2),span(2,end-2),span(2,end-2))%dxPsi;
    cube dyPsiOvPsi = invPsi(span(2,end-2),span(2,end-2),span(2,end-2))%dyPsi;
    cube dzPsiOvPsi = invPsi(span(2,end-2),span(2,end-2),span(2,end-2))%dzPsi;
    invPsi.save("invPsi.dat", raw_ascii);
    dxPsi.save("dxPsi.dat", raw_ascii);
    modGradPsi.save("modGradPsi.dat", raw_ascii);
    dxPsiOvPsi.save("dxPsiOvPsi.dat", raw_ascii);

    // Impose mobility
    cube M = ones(size(Cini));

    // Solver
    lint keepFfunc = 0;
    lint keepMeasuredThea =0;
    lint showIm = 0;
    cube F(size(Cini));
    cube dF(size(Cini));
    cube C = Cini;

    // Criteria of convergence
    vec Ft_arma;
    std::vector<T> Ft;

    T numIt = 5000.;
    lint maxIter = (lint) (numIt/dt);

    for (lint iTT = 1; iTT < maxIter; ++iTT) {

        cout << "Iteration:  " << iTT << std::endl;

        // Energy functional and its derivatives
        F = (W/2)*C%C%(1-C)%(1-C);
        dF = W*(1-C)%(1-2*C)%C;

        // Derivatives of mobility
        cube dxM = gradx(M,h);
        cube dyM = grady(M,h);
        cube dzM = gradz(M,h);

        // Derivatives order parameter
        cube dxC = gradx(C,h);
        cube dyC = grady(C,h);
        cube dzC = gradz(C,h);

        // Intermediate g function
        cube g = (dxPsiOvPsi%dxC + dyPsiOvPsi%dyC + dzPsiOvPsi%dzC) +
                 d2x_3D(C,h) +
                 (modGradPsi%invPsi(span(2,end-2),span(2,end-2),span(2,end-2))%sqrt(2*F(span(2,end-2),span(2,end-2),span(2,end-2)))/eps)*cosTheta;

        // Auxiliary term
        cube k = dF(span(2,end-2),span(2,end-2),span(2,end-2)) - (eps*eps)*g;

        // Finish here the equation
        cube dxk = gradx(k,h);
        cube dyk = grady(k,h);
        cube dzk = gradz(k,h);

        cube dCdt = M(span(4,end-4),span(4,end-4),span(4,end-4))%(dxPsiOvPsi(span(2,dend-2),span(2,dend-2),span(2,dend-2))%dxk +
                                                            dyPsiOvPsi(span(2,dend-2),span(2,dend-2),span(2,dend-2))%dyk +
                                                            dzPsiOvPsi(span(2,dend-2),span(2,dend-2),span(2,dend-2))%dzk) +
                    M(span(4,end-4),span(4,end-4),span(4,end-4))%d2x_3D(k,h) +
                                                            (dxM(span(2,dend-2),span(2,dend-2),span(2,dend-2))%dxk +
                                                            dyM(span(2,dend-2),span(2,dend-2),span(2,dend-2))%dyk +
                                                            dzM(span(2,dend-2),span(2,dend-2),span(2,dend-2))%dzk);

	// Forward euler integration
	C(span(4,end-4),span(4,end-4),span(4,end-4)) = C(span(4,end-4),span(4,end-4),span(4,end-4)) + dt*dCdt;

	// Boundary conditions (No-Flux)
	if (bcType == 1) {
        C(span(1,3),span(4,end-4),span(4,end-4)) = C(span(5,7),span(4,end-4),span(4,end-4));
        C(span(end-3,end-1),span(4,end-4),span(4,end-4)) = C(span(end-7,end-5),span(4,end-4),span(4,end-4));
        C(span(4,end-4),span(5,7),span(4,end-4)) = C(span(4,end-4),span(1,3),span(4,end-4));
        C(span(4,end-4),span(end-3,end-1),span(4,end-4)) = C(span(4,end-4),span(end-7,end-5),span(4,end-4));
        C(span(4,end-4),span(4,end-4),span(5,7)) = C(span(4,end-4),span(4,end-4),span(1,3));
        C(span(4,end-4),span(4,end-4),span(end-3,end-1)) = C(span(4,end-4),span(4,end-4),span(end-7,end-5));

}

	// Check on the internal energy functional
	cube modGradC = sqrt(dxC%dxC+dyC%dyC+dzC%dzC);

	if (iTT%100 == 0){

        ++keepFfunc;
        cube Ffunc = F(span(2,end-2),span(2,end-2),span(2,end-2)) + ((eps*eps/2)*modGradC%modGradC);
        T SumFfunc = accu(Ffunc);

        Ft.push_back(SumFfunc);
        cout << "Ft:  " << Ft[keepFfunc -1] << std::endl;

}




} // End of for loop


/*cout << "Writing output...  " << std::endl;

Ft_arma = conv_to<vec>::from(Ft);
Ft_arma.save(localFolder + "Ft.dat", raw_ascii);
vec Mtheta = conv_to<vec>::from(MeasuredTheta);
Mtheta.save(localFolder + "MTheta.dat", raw_ascii);
C.save(localFolder + "Cend.dat", raw_ascii);*/



} // End of main
