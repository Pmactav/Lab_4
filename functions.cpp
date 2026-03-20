//
// Created by Peter on 2/4/2026.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>
#include "functions.h"

using namespace std;
using namespace Eigen;

MatrixXd ReadDatatoMatrix(const std::string &filename) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "ERROR: Could not open file " << filename << "\n";
        return MatrixXd(0,0);}
    vector<vector<double>> data;
    string line;
    while (getline(infile, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        vector<double> row;
        double value;
        while (ss >> value) {row.push_back(value);}
        if (!row.empty())data.push_back(row);}
    if (data.empty())
        return MatrixXd(0,0);
    int rows = data.size();
    int cols = data[0].size();
    MatrixXd M(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            M(i,j) = data[i][j];
    return M;
}

void WriteMatrixToFile(const MatrixXd &Mat, const string& filename, unsigned int precision){
    ofstream out(filename, ios::out);
    if (out.fail()){
        cout << "Could not open output file " << filename << endl;
        exit(1);
    }
    out << fixed << setprecision(precision);
    for (int i = 0; i < Mat.rows(); ++i){
        for (int j = 0; j < Mat.cols(); ++j){
            out << Mat(i, j);
            if (j != Mat.cols() - 1)
                out << " ";
        }
        out << endl;
    }
    out.close();
}

MatrixXd DesignMatrix(const MatrixXd &ctrlPts, const MatrixXd& x_hat) {
    MatrixXd A(ctrlPts.rows(),2);
    double x0 = x_hat(0,0);
    double y0 = x_hat(1,0);
    for (int i = 0; i < ctrlPts.rows(); i++) {
        double deltaX = (x0 - ctrlPts(i, 0));
        double deltaY = (y0 - ctrlPts(i, 1));
        double denom = sqrt(deltaX * deltaX + deltaY * deltaY);
        A(i, 0) = deltaX / denom;
        A(i, 1) = deltaY / denom;
    }
    return A;
}

MatrixXd Misclosure(const MatrixXd &l, const MatrixXd &ctrlPts, const MatrixXd& x_hat) {
    MatrixXd W(ctrlPts.rows(), 1);
    double x0 = x_hat(0,0);
    double y0 = x_hat(1,0);
    for (int i = 0; i < ctrlPts.rows(); i++) {
        double deltaX = (x0 - ctrlPts(i, 0));
        double deltaY = (y0 - ctrlPts(i, 1));
        double diff = sqrt(deltaX * deltaX + deltaY * deltaY) - l(i, 0);
        W(i, 0) = diff;
    }
    return W;
}

MatrixXd DesignMatrixAz(const MatrixXd &ctrlPts, const MatrixXd &x_hat) {
    int n = ctrlPts.rows();
    MatrixXd Aaz(n, 2);
    double x0 = x_hat(0, 0);
    double y0 = x_hat(1, 0);
    for (int i = 0; i < n; i++) {
        double deltaX = (x0 - ctrlPts(i, 0));
        double deltaY = (y0 - ctrlPts(i, 1));
        double denom = deltaX * deltaX + deltaY * deltaY;
        Aaz(i, 0) = -deltaY / denom;
        Aaz(i, 1) = deltaX / denom;
    }
    return Aaz;
}

MatrixXd MisclosureAz(const MatrixXd &az, const MatrixXd &ctrlPts, const MatrixXd &x_hat) {
    int n = ctrlPts.rows();
    VectorXd Waz(n);
    double x0 = x_hat(0, 0);
    double y0 = x_hat(1, 0);
    for (int i = 0; i < n; i++) {
        double deltaX = (x0 - ctrlPts(i, 0));
        double deltaY = (y0 - ctrlPts(i, 1));
        double theta = atan2(-deltaX, -deltaY);
        double wi = az(i, 0)-theta;
        wi = wi - 2*M_PI*floor((wi + M_PI) / (2 * M_PI));
        Waz(i)=wi;}
    return Waz;
}
