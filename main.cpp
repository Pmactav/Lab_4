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
        delta = -N.inverse()*u;
        cout << "delta: " << endl;
        cout << delta << endl;
        x_hat += delta;}
    while (delta.cwiseAbs().maxCoeff() > stop);
    cout << "Approximated Coordinates of P: " << endl;
    cout << x_hat << endl;
    A=DesignMatrix(ctrlPts,x_hat); //repopulate outside the loop
    w=Misclosure(l,ctrlPts,x_hat); // as above
    MatrixXd N = A.transpose()*P*A;
    VectorXd v_hat = A*delta+w; //update v_hat
    double sigma0Hat_sq = (v_hat.transpose()*P*v_hat)(0,0)/(l.rows() - 2);
    //cout << "a-posteriori: " << sigma0Hat_sq << endl;
    VectorXd l_hat = l+v_hat;
    //cout << "l-hat: " << endl << l_hat << endl;
    MatrixXd Cx_hat = sigma0Hat_sq*N.inverse();
    MatrixXd Cl_hat = A*Cx_hat*A.transpose();
    MatrixXd Cl = P.inverse();
    MatrixXd Cv = Cl-Cl_hat;
    MatrixXd check = Misclosure(l_hat,ctrlPts,x_hat).array().abs();
    //cout << "Largest Misclosure: " << check.cwiseAbs().maxCoeff() << endl;
    // Task 3
    MatrixXd az = ReadDatatoMatrix("../az_rad.txt");
    cout << "az: " << endl << az << endl;
    double azStDevSec = 20.0;
    double azStDevRad = (azStDevSec/3600.0)*M_PI/180.0;
    cout << "azStDevRad: " << azStDevRad << endl;
    double azWeight = 1.0/(azStDevRad*azStDevRad);
    cout << "azWeight: " << azWeight << endl;
    MatrixXd Paz = MatrixXd::Identity(az.rows(), az.rows()) * azWeight;
    MatrixXd Aaz = DesignMatrixAz(ctrlPts, x_hat);
    MatrixXd Naz = Aaz.transpose()*Paz*Aaz;
    VectorXd waz = MisclosureAz(az, ctrlPts, x_hat);
    MatrixXd uaz = Aaz.transpose()*Paz*waz;
    VectorXd deltaAz = -Naz.inverse()*uaz;
    cout << waz << endl;

    return 0;
}
