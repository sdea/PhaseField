/*
    Armadillo  project
    Based on the linear algebra C++ Library Armadillo:
    Armadillo (http://arma.sourceforge.net/)
    Created by Salvatore De Angelis on 1/9/17.
    Copyright © 2017 Salvatore De Angelis. All rights reserved.
*/

#include <iostream>
#include <string> 
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

/// Surface mobility function

cube SM(T Mob, cube C, cube Psi) {

    cube M = zeros(size(Psi));
        for(int i = 0; i < Psi.n_rows; i++){
            for(int j = 0; j < Psi.n_cols; j++){
                for(int k = 0; k < Psi.n_slices; k++){
                    if (C(i,j,k) > 0.1 && C(i,j,k) < 0.9){

                M(i,j,k) = Mob*C(i,j,k)*C(i,j,k)*(1-C(i,j,k)*C(i,j,k))*pow(Psi(i,j,k),6)*((10*Psi(i,j,k)*Psi(i,j,k))*(10*Psi(i,j,k)*Psi(i,j,k))-15*Psi(i,j,k)+6);

                }// End if
            } // End for depth
        } // End for cols
    } // End for rows

    return (M);

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

 void print3Dmat(cube& matrix, std::string filename, int dim) {

        cube exp_geom = reshape(matrix, 1,1, dim*dim*dim);
        mat col_geom = exp_geom.slice(2);
        exp_geom.save(filename, arma::raw_ascii);

 }

int main(int argc, const char * argv[]) {

    // Create directory related to delta
    std::string def_str;

    // Making directory of output
    def_str.append("mkdir output");
    system(def_str.c_str());

    // Local (in loop) output folder
    std::string outputFolder;
    outputFolder.append("output/");
    // localFolder.append(str1); localFolder.append(argv[1]);
    // localFolder.append("_"); localFolder.append(argv[2]); localFolder.append("/");

    // Print important parameters
    // printHeader(N,delta,localFolder);
    // printHeader(localFolder);


    cube geom;
    std::string nameGeom = "geom.bin";

    geom.load(nameGeom,arma::raw_binary);
    geom.reshape(101,101,101);

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
    cube invPsi = 1./Psi;
    cube dxPsi = gradx(Psi,h);
    cube dyPsi = grady(Psi,h);
    cube dzPsi = gradz(Psi,h);
    cube modGradPsi = sqrt(dxPsi%dxPsi + dyPsi%dyPsi + dzPsi%dzPsi);
    cube dxPsiOvPsi = invPsi(span(2,end-2),span(2,end-2),span(2,end-2))%dxPsi;
    cube dyPsiOvPsi = invPsi(span(2,end-2),span(2,end-2),span(2,end-2))%dyPsi;
    cube dzPsiOvPsi = invPsi(span(2,end-2),span(2,end-2),span(2,end-2))%dzPsi;

    print3Dmat(invPsi, outputFolder + "invPsi.dat", 101);
    print3Dmat(dxPsi, outputFolder + "dxPsi.dat", 97);
    print3Dmat(dxPsiOvPsi, outputFolder + "dxPsiOvPsi.dat", 97);
    print3Dmat(modGradPsi, outputFolder + "modGradPsi.dat", 97);

    // Mobility
    cube M(size(Cini));
    cube dMC(size(Cini));
    T Mob = 6.;
    cube G(size(Cini));
    cube dGPsi(size(Cini));

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

        // Mobility function
        M = SM(Mob,C,Psi);

        // Derivatives of mobility
        dMC = 2*C-4*C%C%C;
        dMC.elem(find(dMC > 0.95 || dMC < 0.05)).zeros();

        // Interpolation function G
        G = pow(Psi,6)%(10*Psi%Psi-15*Psi+6);
        dGPsi = pow(Psi,5)%(80*Psi%Psi-105*Psi+36);


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
    ((Mob*(dMC(span(4,end-4),span(4,end-4),span(4,end-4)))%(dxC(span(2,dend-2),span(2,dend-2),span(2,dend-2)))%(G(span(4,end-4),span(4,end-4),span(4,end-4)))) +
    (M(span(4,end-4),span(4,end-4),span(4,end-4))%(dGPsi(span(4,end-4),span(4,end-4),span(4,end-4)))%(dxPsi(span(2,dend-2),span(2,dend-2),span(2,dend-2))))%dxk +
    Mob*(dMC(span(4,end-4),span(4,end-4),span(4,end-4)))%(dyC(span(2,dend-2),span(2,dend-2),span(2,dend-2)))%(G(span(4,end-4),span(4,end-4),span(4,end-4))) +
    (M(span(4,end-4),span(4,end-4),span(4,end-4))%(dGPsi(span(4,end-4),span(4,end-4),span(4,end-4)))%(dyPsi(span(2,dend-2),span(2,dend-2),span(2,dend-2))))%dyk +
    Mob*(dMC(span(4,end-4),span(4,end-4),span(4,end-4)))%(dzC(span(2,dend-2),span(2,dend-2),span(2,dend-2)))%(G(span(4,end-4),span(4,end-4),span(4,end-4))) +
    (M(span(4,end-4),span(4,end-4),span(4,end-4))%(dGPsi(span(4,end-4),span(4,end-4),span(4,end-4)))%(dzPsi(span(2,dend-2),span(2,dend-2),span(2,dend-2))))%dzk);

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
     	if(iTT%10 == 0) {

    	  print3Dmat(C, "output/C_" + std::to_string(iTT) + ".dat", 101);

	}


} // End of for loop





} // End of main
