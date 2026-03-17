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
    MatrixXd A = DesignMatrix(ctrlPts, 250, 200);
    VectorXd sigma = .001+2e-6*l.array();
    MatrixXd P = sigma.array().square().inverse().matrix().asDiagonal();
    MatrixXd N = A.transpose()*P*A;
    MatrixXd w = Misclosure(l, ctrlPts, 250, 200);
    MatrixXd u = A.transpose()*P*w;
    MatrixXd delta = -N.inverse()*u;
    cout << delta << endl;
}
