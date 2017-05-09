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
mat d2x_2D(mat const& C, T const& h) {


    lint end = C.n_rows -1;

    // Computation using compact stencil
    mat out1 = C(span(3,end-1),span(2,end-2)) + C(span(1,end-3),span(2,end-2));
    mat out2 = C(span(2,end-2),span(3,end-1)) + C(span(2,end-2), span(1,end-3));
    mat outCenter = 4*C(span(2,end-2),span(2,end-2));


    return (out1 - outCenter + out2)/(h*h);

}

/// Function for gradient in 2D
mat gradx(mat const& C, T const& h) {

    // !!!! Only works for square matrix !!!
    lint end = C.n_rows -1;
    mat Cx = C(span(3,end-1), span(1, end-3)) - C(span(1,end-3),span(1,end-3));

    return Cx/(2*h);
}

mat grady(mat const& C, T const& h) {

    // !!!! Only works for square matrix !!!
    lint end = C.n_rows -1;
    mat Cy = C(span(1, end-3),span(3,end-1)) - C(span(1,end-3),span(1,end-3));

    return Cy/(2*h);
}



/// BiHarmonic operator
mat d4x_2D(mat const& C, T const& h) {

    lint end = C.n_rows -1;

    //  Computation using compact stencil
    mat out1 = C(span(4,end),span(2,end-2)) - 4*C(span(3,end-1),span(2,end-2)) - 4*C(span(1,end-3),span(2,end-2)) + C(span(0,end-4),span(2,end-2));
    mat out2 = C(span(2,end-2),span(4,end)) -4*C(span(2,end-2),span(3,end-1)) - 4*C(span(2,end-2),span(1,end-3)) + C(span(2,end-2),span(0,end-4));
    mat outCenter = 12*C(span(2,end-2),span(2,end-2));

    T h2 = h*h;

    return (out1 + outCenter + out2)/(h2*h2);

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
    if (argc != 3) {

        cout <<"Wrong number of parameters!!!" << std::endl;
        cout <<"Required: (int) N, (T) delta" << std::endl;

        return -1;

    }


    // Reading parameters from command line
    lint N = atoi(argv[1]);
    T delta = atof(argv[2]);

    // Create directory related to delta
    std::string str1 = "output_";
    std::string def_str;

    // Making directory for every delta and dimension
    def_str.append("mkdir output/"); def_str.append(str1); def_str.append(argv[1]);
    def_str.append("_"); def_str.append(argv[2]);
    system(def_str.c_str());

    // Local (in loop) output folder
    std::string localFolder;
    localFolder.append("output/");
    localFolder.append(str1); localFolder.append(argv[1]);
    localFolder.append("_"); localFolder.append(argv[2]); localFolder.append("/");





    // Print important parameters
    printHeader(N,delta,localFolder);


    // lint N = 100;
    T h = 0.2;
    T dt = 1e-5;
    lint Time = 1;

    T Q = 1.;
    lint bcType = 1;
    T eps = 1.;

    // Set geometry
    mat geom(N+4,N+4);
    geom.zeros();

    // Set last number of arrays
    lint end = geom.n_rows -1;
    lint dend = geom.n_rows -5;

    // Geometry using subslicing
    geom(span(((N+4)/2-1),end), span(0,end)).ones();
    geom(span(0,((N+4)/2-2)), span(0,((N+4)/2)-1)).fill(2.);
    geom.save(localFolder + "InitialGeometry.dat", raw_ascii);

    // Contact angle values
    T theta = ((T)100/(T)180)*datum::pi;
    T cosTheta = cos(theta);

    // Print contact angles
    cout << "\nTheta:  " << theta << "  Cos(theta):  \n" << cosTheta <<std::endl;

    // Make order parameter Psi
    Col<T> x_vec(N+4);
    x_vec = linspace<vec>(1,N+4,(N+4));
    T x0 = (T)(N+4)/2;

    // Popolate the matrix
    mat Psi(N+4,N+4);
    vec psi_x = (1-tanh((x_vec - x0)/delta))/2;
    Psi = repmat(psi_x,1,N+4);
    Psi.save(localFolder + "Psi.dat", raw_ascii);

    // Make the C matrix (order parameter)
    mat C(N+4,N+4);
    C.zeros();
    C(span(0, end), span(0,((N+4)/2)-1)).ones();
    C.save(localFolder + "Cini.dat", raw_ascii);

    // Compute the gradient of Psi and the term dPsi/Psi
    mat invPsi = 1./Psi;
    mat dxPsi = gradx(Psi, h);
    mat dyPsi = grady(Psi,h);
    mat modGradPsi = sqrt(dxPsi%dxPsi + dyPsi%dyPsi);
    mat dxPsiOvPsi = invPsi(span(2,end-2),span(2,end-2))%dxPsi;
    mat dyPsiOvPsi = invPsi(span(2,end-2),span(2,end-2))%dyPsi;
    invPsi.save(localFolder + "invPsi.dat", raw_ascii);
    dxPsi.save(localFolder + "dxPsi.dat", raw_ascii);
    modGradPsi.save(localFolder + "modGradPsi.dat", raw_ascii);
    dxPsiOvPsi.save(localFolder + "dxPsiOvPsi.dat", raw_ascii);

    // Impose mobility and tpb regions
    mat M = ones(N+4,N+4);
    mat tpbRegion(size(C(span(4,end-4),span(4,end-4))));

    // Solver
    lint keepAngle = 0;
    lint keepMeasuredThea =0;
    lint showIm = 0;
    mat F(N+4,N+4);
    mat dF(N+4,N+4);
    mat Ffunc;
    T SumFfunc;

    // Criteria of convergence
    vec Ft_arma;
    std::vector<T> Ft;
    T new_meas_theta = 1;
    T old_meas_theta = 1;
    std::vector<T> MeasuredTheta;
    rowvec Measuredtheta_arma;
    T meas_theta;
    T contact_cos;
    vec vec_contact;

    // Convergenze criteria
    T error = 1.;

    // Max iteration
    T numIt = 5000.;
    lint maxIter = (lint) (numIt/dt);
    lint afterN = 100000;

    for (lint iTT = 1; iTT < maxIter; ++iTT) {

        //cout << "Iteration:  " << iTT*dt << std::endl;

        // Energy functional and its derivatives
        F = (Q/4)*C%C%(1-C)%(1-C);
        dF = (Q/2)*(1-C)%(1-2*C)%C;

        // Derivatives of mobility
        mat dxM = gradx(M,h);
        mat dyM = grady(M,h);

        // Derivatives order parameter
        mat dxC = gradx(C,h);
        mat dyC = grady(C,h);

        // Intermediate g function
        mat g = (dxPsiOvPsi%dxC + dyPsiOvPsi%dyC) +
                 d2x_2D(C,h) +
                 (modGradPsi%invPsi(span(2,end-2),span(2,end-2))%sqrt(2*F(span(2,end-2),span(2,end-2)))/eps)*cosTheta;

        // Auxiliary term
        mat k = dF(span(2,end-2),span(2,end-2)) - (eps*eps)*g;

        // Finish here the equation
        mat dxk = gradx(k,h);
        mat dyk = grady(k,h);

        mat dCdt = M(span(4,end-4),span(4,end-4))%(dxPsiOvPsi(span(2,dend-2),span(2,dend-2))%dxk +
                                                            dyPsiOvPsi(span(2,dend-2),span(2,dend-2))%dyk)+
                   M(span(4,end-4), span(4,end-4))%d2x_2D(k,h) +
                   (dxM(span(2,dend-2),span(2,dend-2))%dxk + dyM(span(2,dend-2),span(2,dend-2))%dyk);

	// Forward euler integration
	C(span(4,end-4),span(4,end-4)) = C(span(4,end-4),span(4,end-4)) + dt*dCdt;

	// Boundary conditions (No-Flux)
	if (bcType == 1) {
        C(span(1,3),span(4,end-4)) = C(span(5,7),span(4,end-4));
        C(span(end-3,end-1),span(4,end-4)) = C(span(end-7,end-5),span(4,end-4));
        C(span(4,end-4),span(5,7)) = C(span(4,end-4),span(1,3));
        C(span(4,end-4),span(end-3,end-1)) = C(span(4,end-4),span(end-7,end-5));

}

	// Check on the internal energy functional
	mat modGradC = sqrt(dxC%dxC+dyC%dyC);




	if (iTT%10000 == 0){

        ++keepAngle;
        Ffunc = F(span(2,end-2),span(2,end-2)) + ((eps*eps/2)*modGradC%modGradC);
        SumFfunc = accu(Ffunc);

        Ft.push_back(SumFfunc);

        }



    if (iTT%1000 == 0) {

        ++keepMeasuredThea;
        tpbRegion.zeros();
        tpbRegion.elem(find(C(span(4,end-4),span(4,end-4)) > 0.2 && C(span(4,end-4),span(4,end-4)) < 0.8 &&
                            Psi(span(4,end-4),span(4,end-4)) > 0.2 && Psi(span(4,end-4),span(4,end-4)) < 0.8)).ones();

        mat den = ((1/modGradC)%(1/modGradPsi));
        mat num = (dxC%dxPsi+dyC%dyPsi);
        mat MeasuredCos = num%den;

        // Select area where to compute conctact angle
        vec_contact = MeasuredCos.elem(find(tpbRegion == 1.));
        contact_cos = mean(vec_contact);
        meas_theta = convToDeg(acos(-contact_cos));
        MeasuredTheta.push_back(meas_theta);

        // Convergence criteria
        new_meas_theta = meas_theta;
        error = abs(new_meas_theta - old_meas_theta);

        old_meas_theta = new_meas_theta;

        uvec region = find(tpbRegion == 1);
        cout << "Iteration:  " << iTT  <<"  Measured theta:  " << meas_theta
                <<"  Region:  " << region.n_elem <<"  Error:  " << error <<  std::endl;


        if (iTT > afterN) {
            if (error < 1e-06) {

                cout << "Simulation finished!!!" << std::endl;
                break;
            }
        }

	}

} // End of for loop


cout << "Writing output...  " << std::endl;

Ft_arma = conv_to<vec>::from(Ft);
Ft_arma.save(localFolder + "Ft.dat", raw_ascii);
vec Mtheta = conv_to<vec>::from(MeasuredTheta);
Mtheta.save(localFolder + "MTheta.dat", raw_ascii);
C.save(localFolder + "Cend.dat", raw_ascii);



} // End of main
