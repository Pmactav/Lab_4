#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <numeric>
#include <vector>
#include "functions.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main() {
    MatrixXd az = ReadDatatoMatrix("../az_ddeg.txt");
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
    A=DesignMatrix(ctrlPts,x_hat); //make sure most recent values is used
    w=Misclosure(l,ctrlPts,x_hat); // as above
    VectorXd v_hat = A*delta+w; //update v_hat
    double sigma0Hat_sq = (v_hat.transpose()*P*v_hat)(0,0)/(l.rows() - 2);
    cout << "a-posteriori: " << sigma0Hat_sq << endl;
    VectorXd l_hat = l*v_hat;
    return 0;
}
