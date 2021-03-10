#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <vector>
#include <iomanip>

using namespace Eigen;

long double start = 0.; // segment of solution
long double finish = 1.;

void preciseSolution(std::array<long double, 11>& y1_pr, std::array<long double, 11>& y2_pr, Matrix<long double, 2, 1> initCondition) {
	long double A = initCondition(0, 0);
	long double B = initCondition(1, 0);
	long double C1 = A / 10 - B/ 4;
	long double C2 = A / 10 + B / 4;
	std::cout << A << "; " << B << std::endl;
	std::cout << C1 << " ; " << C2 << std::endl;
	long double x = 0.;
	long double step = .1;
	for (int i = 0; i < 11; ++i) {
		y1_pr[i] = C1 * 5 * exp(-199. * x) + C2 * 5 * exp(x);
		y2_pr[i] = C1 * (-2) * exp(-199 * x) + C2 * 2 * exp(x);
		x += step;
	}
}

void createGrid(int numberOfKnots, std::vector<long double>& grid) {
	grid.clear();
	long double h = (finish - start) / (numberOfKnots - 1);
	for (int i = 0; i < numberOfKnots; ++i)
		grid.push_back(start + i * h);
}

long double findStep(int numberOfKnots) {
	return (finish - start) / (numberOfKnots - 1);
}

template <typename T>
void printVector(T& v) {
	for (auto it = v.begin(); it < v.end(); ++it)
		std::cout << *it << " ; ";
	std::cout << std::endl;
}

void solve(std::array<long double, 11>& y1, std::array<long double, 11>& y2,long double h, int numberOfKnots, Matrix<long double, 2, 2> A, Matrix<long double, 2, 1> y) { // y - initCondition
	int k = (numberOfKnots - 1)/10; // step to add in solution array
	//std::cout << "k = " << k << "; h = " << h << std::endl;
	y1[0] = y(0, 0);
	y2[0] = y(1, 0);
	int kCurrent = 0;
	Matrix<long double, 2, 2> M;
	M = (Matrix<long double, 2, 2>::Identity() - A*h).inverse();
	for (int i = 1; i < numberOfKnots; ++i) {
		//std::cout << y(0, 0) << " ; ";
		y = M * y;
		if (i % k == 0) {
			kCurrent++;
			y1[kCurrent] = y(0, 0);
			y2[kCurrent] = y(1, 0);
		}
	}
	//std::cout << y(0, 0) << std::endl;
}

long double findDelta(std::array<long double, 11>& y, const std::array<long double, 11>& y_pr) {
	long double delta = fabsl(y[0] - y_pr[0]);
	for (int i = 1; i < 11; ++i) {
		if (fabsl(y[i] - y_pr[i]) > delta)
			delta = fabsl(y[i] - y_pr[i]);
	}
	return delta;
}

template <typename T>
void outputGenerator(T& y1, T& y1_pr, T& y2, T& y2_pr, int numberOfKnots) {
	std::ofstream fout("output_final_res.txt");
  	//fout.setf(std::setw(10));
    	//fout.left();
	fout << "Number of knots = " << numberOfKnots << std::endl;
	fout << "A = -5*10^(-15); B = 2*10^(-15)" << std::endl;
	for (int i = 0; i < 11; ++i) {
		fout << i + 1 << "	";
		//fout << std::setw(10);
		fout << std::fixed << std::setprecision(2) << .1 * i << "	";
		//fout << std::setw(10);
		fout << std::fixed << std::scientific << std::setprecision(5) << y1[i] << "	"; 
                fout << std::fixed << std::scientific << std::setprecision(5) << y1_pr[i] << "	";
		fout << std::fixed << std::scientific << std::setprecision(5) << fabsl(y1[i] - y1_pr[i]) << "	";
               	fout << std::fixed << std::scientific << std::setprecision(5) << y2[i] << "	";
                fout << std::fixed << std::scientific << std::setprecision(5) << y2_pr[i] << "	";
                fout << std::fixed << std::scientific << std::setprecision(5) << fabsl(y2[i] - y2_pr[i]);
		fout << std::endl;		
	}
}

int main(int argc, char* argv[])
{
	std::cout.setf(std::ios::fixed);
	std::cout.precision(8);
	Matrix<long double, 2, 1> initCondition(-5*pow(10., 15), 2.*pow(10., 15));
	Matrix<long double, 2, 2> A;
	A(0, 0) = -99.;
	A(0, 1) = 250;
	A(1, 0) = 40;
	A(1, 1) = -99;
	std::array<long double, 11> y1;
	std::array<long double, 11> y2;
	std::array<long double, 11> y1_pr;
        std::array<long double, 11> y2_pr;
	preciseSolution(y1_pr, y2_pr, initCondition);
	int numberOfKnots = 11;
	long double delta = 100;
	std::ofstream fout("output_final.txt");
	fout << "Number of knots:	Delta:" << std::endl;
	numberOfKnots = 671088641;
	//while (fabsl(delta) > pow(10., -6)) {
		fout << numberOfKnots << "		        ";
		solve(y1, y2, findStep(numberOfKnots), numberOfKnots, A, initCondition);
		delta = findDelta(y1, y1_pr);
		if (findDelta(y2, y2_pr) > delta)
			delta = findDelta(y2, y2_pr);
		numberOfKnots = (numberOfKnots - 1) * 2 + 1;
		fout << std::fixed << std::scientific << std::setprecision(5) << delta << std::endl;
	//}
	std::cout << "numberOfKnots = " << (numberOfKnots - 1) / 2 + 1 << std::endl;
	std::cout << "delta = " << delta << std::endl;
	std::cout << "number solutions" << std::endl;
	printVector(y1);
	printVector(y2);
	preciseSolution(y1_pr, y2_pr, initCondition);
	std::cout << "precise solutions" << std::endl;
	printVector(y1_pr);
	printVector(y2_pr);
	outputGenerator(y1, y1_pr, y2, y2_pr, (numberOfKnots - 1) / 2 + 1);
	return 0;
}

