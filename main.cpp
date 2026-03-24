#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <numeric>
#include <vector>
#include <cmath>
#include "functions.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main() {
    //Task 2
    MatrixXd ctrlPts = ReadDatatoMatrix("../ctrlPoints_normalSpaces.txt");
    MatrixXd l = ReadDatatoMatrix("../distObs_2026.txt");
    VectorXd sigma = .001+2e-6*l.array();
    MatrixXd P = sigma.array().square().inverse().matrix().asDiagonal();
    MatrixXd x_hat(2,1);
    x_hat << 250,
             200;
    double stop = 1e-5;
    VectorXd delta, w;
    MatrixXd A;
    do {
        MatrixXd A = DesignMatrix(ctrlPts, x_hat);
        MatrixXd N = A.transpose()*P*A;
        MatrixXd w = Misclosure(l, ctrlPts, x_hat);
        MatrixXd u = A.transpose()*P*w;
        delta = -N.inverse()*u;//during testing cout delta to confirm it is converging
        x_hat += delta;}
    while (delta.cwiseAbs().maxCoeff() > stop);
    cout << "Last delta: " << endl << delta << endl; //should be 0,0
    cout << "Approximated Coordinates of P: " << endl;
    cout << x_hat << endl;
    A=DesignMatrix(ctrlPts,x_hat); //repopulate outside the loop with correction
    w=Misclosure(l,ctrlPts,x_hat); // as above
    MatrixXd N = A.transpose()*P*A;
    VectorXd v_hat = A*delta+w;
    double sigma0Hat_sq = (v_hat.transpose()*P*v_hat)(0,0)/(l.rows() - 2);
    cout << "a-posteriori distance: " << sigma0Hat_sq << endl;
    VectorXd l_hat = l+v_hat;
    WriteMatrixToFile(l_hat, "../l_hat.txt", 10);
    MatrixXd Cx_hat = sigma0Hat_sq*N.inverse();
    VectorXd stdX = Cx_hat.diagonal().cwiseAbs().array().sqrt().matrix();
    WriteMatrixToFile(stdX, "../Cx_hat.txt", 10);
    MatrixXd Cl_hat = A*Cx_hat*A.transpose();
    VectorXd stdL = Cl_hat.diagonal().cwiseAbs().array().sqrt().matrix();
    WriteMatrixToFile(stdL, "../Cl_hat.txt", 10);
    MatrixXd Cl = P.inverse();
    MatrixXd Cv = Cl-Cl_hat;
    VectorXd stdV = Cv.diagonal().cwiseAbs().array().sqrt().matrix();
    WriteMatrixToFile(stdV, "../Cv.txt", 10);
    MatrixXd check = Misclosure(l_hat,ctrlPts,x_hat).array().abs();
    WriteMatrixToFile(v_hat, "../v_hat_coords.txt", 10);
    cout << "Largest Misclosure: " << check.cwiseAbs().maxCoeff() << endl;
    double rho_x_hat = Cx_hat(0,1)/(sqrt(Cx_hat(0,0))*sqrt(Cx_hat(1,1)));
    cout << "rho_x_hat: " << rho_x_hat << endl;
    // Task 3
    MatrixXd laz = ReadDatatoMatrix("../az_rad.txt");
    double azStDevSec = 20.0;
    double azStDevRad = (azStDevSec/3600.0)*M_PI/180.0;
    double azWeight = 1.0/(azStDevRad*azStDevRad);
    MatrixXd Paz = MatrixXd::Identity(laz.rows(), laz.rows()) * azWeight;
    VectorXd deltaAz;
    MatrixXd x_hatAz(2,1);
    x_hatAz << 250,
               200;
    do {
        MatrixXd Aaz = DesignMatrixAz(ctrlPts, x_hatAz);
        MatrixXd Naz = Aaz.transpose()*Paz*Aaz;
        VectorXd waz = MisclosureAz(laz, ctrlPts, x_hatAz);
        MatrixXd uaz = Aaz.transpose()*Paz*waz;
        deltaAz = -Naz.inverse()*uaz;//during testing cout delta to confirm it is converging
        x_hatAz += deltaAz;
    } while (deltaAz.cwiseAbs().maxCoeff() > stop);
    cout << "Last deltaAz: " << endl << deltaAz << endl; //should be 0,0
    cout << "Approximated Coordinates of Paz: " << endl;
    cout << x_hatAz << endl;
    //
    MatrixXd Aaz = DesignMatrixAz(ctrlPts,x_hatAz); //repopulate outside the loop with correction
    VectorXd waz = MisclosureAz(laz,ctrlPts,x_hatAz); // as above
    MatrixXd Naz = Aaz.transpose()*Paz*Aaz;
    VectorXd v_hataz = Aaz*deltaAz+waz;
    double sigma0Hat_sqaz = (v_hataz.transpose()*Paz*v_hataz)(0,0)/(laz.rows() - 2);
    cout << "a-posteriori azimuth: " << sigma0Hat_sqaz << endl;
    VectorXd l_hataz = laz+v_hataz;
    VectorXd l_hat_deg = l_hataz * (180.0 / M_PI);
    WriteMatrixToFile(l_hat_deg, "../l_hat_deg.txt", 6);
    MatrixXd Cx_hataz = sigma0Hat_sqaz*Naz.inverse();
    VectorXd stdXaz = Cx_hataz.diagonal().cwiseAbs().array().sqrt().matrix();
    WriteMatrixToFile(stdX, "../Cx_hataz.txt", 10);
    MatrixXd Cl_hataz = Aaz*Cx_hataz*Aaz.transpose();
    VectorXd stdLaz = Cl_hataz.diagonal().cwiseAbs().array().sqrt().matrix();
    WriteMatrixToFile(stdL, "../Cl_hataz.txt", 10);
    MatrixXd Claz = sigma0Hat_sqaz*Paz.inverse();
    MatrixXd Cvaz = Claz-Cl_hataz;
    VectorXd stdVaz = Cvaz.diagonal().cwiseAbs().array().sqrt().matrix();
    WriteMatrixToFile(stdV, "../Cvaz.txt", 10);
    MatrixXd checkaz = MisclosureAz(l_hataz,ctrlPts,x_hatAz).array().abs();
    cout << "Largest Misclosure: " << checkaz.cwiseAbs().maxCoeff() << endl;
    //cout << checkaz << endl;
    double rho_x_hataz = Cx_hataz(0,1)/(sqrt(Cx_hataz(0,0))*sqrt(Cx_hataz(1,1)));
    cout << "rho_x_hataz: " << rho_x_hataz << endl;
    WriteMatrixToFile(v_hataz, "../v_hat_az.txt", 10);

    return 0;
}
