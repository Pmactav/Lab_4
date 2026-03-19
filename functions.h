//
// Created by Peter on 2/4/2026.
//

#ifndef LAB_2_FUNCTIONS_H
#define LAB_2_FUNCTIONS_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

MatrixXd ReadDatatoMatrix(const string& filename);
void WriteMatrixToFile(const MatrixXd& Mat, const string& filename, unsigned int precision);
MatrixXd DesignMatrix(const MatrixXd& ctrlPts, const MatrixXd& x_hat);
MatrixXd Misclosure(const MatrixXd& l,const MatrixXd &ctrlPts, const MatrixXd& x_hat);
MatrixXd DesignMatrixAz(const MatrixXd& ctrlPts, const MatrixXd& x_hat);
MatrixXd MisclosureAz(const MatrixXd& az,const MatrixXd &ctrlPts, const MatrixXd& x_hat);

#endif //LAB_2_FUNCTIONS_H