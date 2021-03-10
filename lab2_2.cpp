#include <iostream>
#include <cmath>
#include <utility>
#include <array>
#include <vector>
#include <fstream>
#include <iomanip>

long double x0 = 1 / sqrt(2);
long double u0 = 0;
long double u1 = 1;
std::ofstream fout("out.txt");

long double k(long double x) {
	return pow(x, 2) + .5;
}

long double q(long double x) {
	if (x < x0) return exp(-pow(x, 2));
	else return 1.;
}

long double f(long double x) {
	if (x < x0) return cos(x);
	else return 1.;
}

template <typename T>
void print(T& v) {
        for (auto it = v.begin(); it < v.end(); ++it)
                std::cout << *it << " ; ";
        std::cout << std::endl;
}


void createGrid(int numberOfKnots, std::vector<long double>& grid, int& leftBorder, int& rightBorder) {
	grid.clear();
	long double h = 1. / (numberOfKnots - 1);
	//std::cout << h << std::endl;
	for (int i = 0; i < numberOfKnots; ++i) {
		if (((i * h) <= x0) && (((i + 1) * h) >= x0)) {
		       	leftBorder = i;
	       		rightBorder = i + 1;
		}		
		grid.push_back(i * h);
	}
	//std::cout << leftBorder << " ; " << rightBorder << std::endl;
	
}

void leftCoeffsStat(std::vector<long double>& al, std::vector<long double>& bl, int numberOfKnots) {
	long double h = 1. / (numberOfKnots - 1);
	al[0] = -k(x0)/(-2*k(x0)-q(x0)*pow(h, 2));
	bl[0] = (-f(x0)*pow(h, 2) - k(x0)*u0) / (-2*k(x0) - q(x0)*pow(h, 2));
	for (int i = 1; i < al.size(); ++i) {
		al[i] = -k(x0)/(-2*k(x0) - q(x0)*pow(h, 2) + k(x0)*al[i-1]);
		bl[i] = (-f(x0)*pow(h, 2) - k(x0)*bl[i-1]) / (-2*k(x0)-q(x0)*pow(h, 2) + k(x0)*al[i-1]);
	}
/*	al[0] = -1/(-2-pow(h, 2));
        bl[0] = (-pow(h, 2)) / (-2 - pow(h, 2));
        for (int i = 1; i < al.size(); ++i) {
                al[i] = -1/(-2 - pow(h, 2) + al[i-1]);
                bl[i] = (-pow(h, 2) - bl[i-1]) / (-2-pow(h, 2) + al[i-1]);
        }*/
}

void rightCoeffsStat(std::vector<long double>& ar, std::vector<long double>& br, int numberOfKnots) {
	long double h = 1. / (numberOfKnots - 1);
	ar[ar.size()-1] = -k(x0)/(-2*k(x0)-q(x0)*pow(h, 2));
	br[br.size()-1] = (-f(x0)*pow(h, 2) - k(x0)*u1) / (-2*k(x0) - q(x0)*pow(h, 2));
        for (int i = ar.size() - 2; i > -1; --i) {
                ar[i] = -k(x0)/(-2*k(x0) - q(x0)*pow(h, 2) + k(x0)*ar[i+1]);
                br[i] = (-f(x0)*pow(h, 2) - k(x0)*br[i+1]) / (-2*k(x0)-q(x0)*pow(h, 2) + k(x0)*ar[i+1]);
        }

}

long double findDelta(std::array<long double, 11>& y, const std::array<long double, 11>& y_pr) {
	long double delta = fabsl(y[0] - y_pr[0]);
	for (int i = 1; i < 11; ++i) {
		if (fabsl(y[i] - y_pr[i]) > delta)
			delta = fabsl(y[i] - y_pr[i]);
	}
	return delta;
}

void statPrecise(std::array<long double, 11>& u) {
	long double x  = 0.;
	long double C2 = -exp(2.) / (exp(2.)-1);
	long double C1 = 1. / (exp(2.)-1);
	for (int i = 0; i < 11; ++i) {
		x = i * .1;
		u[i] = C1*exp(x) + C2*exp(-x) + 1;
	}
}

void statProgonka(int numberOfKnots, std::vector<long double>& grid, std::array<long double, 11>& u) {
	int l = (numberOfKnots - 1)/10;
	int alpha; int beta;
	createGrid(numberOfKnots, grid, alpha, beta);
	//std::cout << alpha << " ; " << beta << std::endl;
	std::vector<long double> u_(numberOfKnots);
	u_[0] = u0;	u_[numberOfKnots-1] = u1;
	//std::array<long double, alpha - 1> al;	std::array<long double, alpha - 1> bl;
	//std::array<long double, numberOfKnots - beta - 2> ar;	std::array<long double, numberOfKnots - beta - 2> br;
	std::vector<long double> al(alpha - 1);		std::vector<long double> bl(alpha - 1);
	std::vector<long double> ar(numberOfKnots - beta - 2);
	std::vector<long double> br(numberOfKnots - beta - 2);
	leftCoeffsStat(al, bl, numberOfKnots);
	rightCoeffsStat(ar, br, numberOfKnots);
	long double temp = k(x0) * (bl[bl.size()-1] + br[0]) / (k(x0) * (2 - al[bl.size()-1] - ar[0]));
	u_[alpha] = temp;	u_[beta] = temp;
	for (int i = alpha - 1; i > 0; --i)
	       u_[i] = al[i-1]*u_[i+1] + bl[i-1];
	for (int i = 0; i < ar.size(); ++i)
		u_[beta + 1 + i] = ar[i]*u_[beta + i] + br[i];
	//print(u_);
	for (int i = 0; i < u_.size(); ++i) {
		if ((i % l) == 0)
			u[i / l] = u_[i];
	}
}

void statSolve(std::array<long double, 11>& u, std::vector<long double>& grid) {
	std::array<long double, 11> u_pr;
	statPrecise(u_pr);
	std::cout << "precise solution" << std::endl;
	print(u_pr);
	int n = 11;
	statProgonka(n, grid, u);
	fout << "knots	delta\n";
	fout << n << "	";
	fout << std::fixed << std::scientific << std::setprecision(4) << findDelta(u, u_pr) << std::endl;
	int k = 1;
	while (fabsl(findDelta(u, u_pr) > .0001)) {
	//while (n < 5000) {
		n = (n - 1)*2 + 1;
		statProgonka(n, grid, u);
		fout << n << "	" << std::scientific << std::setprecision(4) << findDelta(u, u_pr) << std::endl;
	}
	fout << std::endl;
}

template <typename T>
void dinLeftCoeffs(T& al, T& bl, T& cl, T& dl, T& grid) {
	long double h = grid[1] - grid[0];
	for (int i = 0; i < al.size(); ++i) {
		al[i] = k(grid[i+1] + h/2);
		bl[i] = -(k(grid[i+1] + h/2) + k(grid[i+1] - h/2) + q(grid[i+1])*pow(h, 2.));
		cl[i] = k(grid[i+1] - h/2);
		dl[i] = -f(grid[i+1]) * pow(h, 2.);
	}
}

template <typename T>
void dinRightCoeffs(T& ar, T& br, T& cr, T& dr, T& grid, int beta) {
        long double h = grid[1] - grid[0];
        for (int i = 0; i < ar.size(); ++i) {
                ar[i] = k(grid[beta+i+1] + h/2);
                br[i] = -(k(grid[beta+i+1] + h/2) + k(grid[beta+i+2] - h/2) + q(grid[beta+1])*pow(h, 2.));
                cr[i] = k(grid[beta+i+1] - h/2);
                dr[i] = -f(grid[beta+i+1]) * pow(h, 2.);
        }
}


void dinProgonka(int numberOfKnots, std::vector<long double>& grid, std::array<long double, 11>& u) {
	int l = (numberOfKnots - 1)/10;
        int alpha; int beta;
        createGrid(numberOfKnots, grid, alpha, beta);
        //std::cout << alpha << " ; " << beta << std::endl;
        std::vector<long double> u_(numberOfKnots);
        u_[0] = u0;     u_[numberOfKnots-1] = u1;
        //std::array<long double, alpha - 1> al;        std::array<long double, alpha - 1> bl;
        //std::array<long double, numberOfKnots - beta - 2> ar; std::array<long double, numberOfKnots - beta - 2> br;
        std::vector<long double> al(alpha - 1);         
	std::vector<long double> bl(alpha - 1);
	std::vector<long double> cl(alpha - 1);
	std::vector<long double> dl(alpha - 1);
        std::vector<long double> ar(numberOfKnots - beta - 2);
        std::vector<long double> br(numberOfKnots - beta - 2);
	std::vector<long double> cr(numberOfKnots - beta - 2);
        std::vector<long double> dr(numberOfKnots - beta - 2);
        dinLeftCoeffs(al, bl, cl, dl, grid);
        dinRightCoeffs(ar, br, cr, dr, grid, beta);

	std::vector<long double> all(alpha - 1);
	std::vector<long double> bll(alpha - 1);
        std::vector<long double> arr(numberOfKnots - beta - 2);
        std::vector<long double> brr(numberOfKnots - beta - 2);
	all[0] = -al[0]/bl[0];
	bll[0] = dl[0] - (cl[0]*u0) / bl[0];
	for (int i = 1; i < al.size(); ++i) {
		all[i] = -al[i] / (bl[i] + cl[i]*all[i-1]);
		bll[i] = (dl[i] - cl[i]*bll[i-1]) / (bl[i] + cl[i]*all[i-1]);
	}
	arr[ar.size()-1] = -cr[ar.size()-1] / br[ar.size()-1];
	brr[ar.size()-1] = (dr[ar.size()-1] - cr[ar.size()-1]*u1) / br[ar.size()-1];
	for (int i = ar.size() - 2; i > -1; --i) {
		arr[i] = -cr[i] / (br[i] + ar[i]*arr[i+1]);
		brr[i] = (dl[i] - ar[i]*brr[i+1]) / (br[i] + ar[i]*arr[i+1]);
	}


        long double temp = (k(grid[alpha]) * bll[bll.size()-1] + k(grid[beta])  * brr[0]) / (k(grid[alpha]) * (1 - all[bll.size()-1]) + k(grid[beta]) * (1 - ar[0]));
        u_[alpha] = temp;       u_[beta] = temp;
        for (int i = alpha - 1; i > 0; --i)
               u_[i] = all[i-1]*u_[i+1] + bll[i-1];
        for (int i = 0; i < ar.size(); ++i)
                u_[beta + 1 + i] = arr[i]*u_[beta + i] + brr[i];
        for (int i = 0; i < u_.size(); ++i) {
                if ((i % l) == 0)
                        u[i / l] = u_[i];
        }

}

void dinSolve(std::array<long double, 11>& u, std::vector<long double>& grid) {
        int n = 11;
	int n1 = 21;
	std::array<long double, 11> u1;
        dinProgonka(n, grid, u);
	dinProgonka(n1, grid, u1);
	fout << "knots	delta\n";
	fout << n << "	";
	fout << std::fixed << std::scientific << std::setprecision(4) << findDelta(u, u1) << std::endl;
        int k = 1;
        while (fabsl(findDelta(u, u1) > .0001)) {
        //while (n < 5000) {
		n = n1;
		u = u1;
                n1 = (n1 - 1)*2 + 1;
                statProgonka(n1, grid, u1);
		fout << n << "	" << std::scientific << std::setprecision(4) << findDelta(u, u1) << std::endl;
        }
	fout << std::endl;
}


int main() {
	//std::cout.setf(std::ios::fixed);
	//std::cout.precision(8);
	std::vector<long double> grid;
	std::array<long double, 11> u;
	std::array<long double, 11> u_pr;
	fout << "STATIC RESULTS\n";
	statPrecise(u_pr);
	statSolve(u, grid);
	//std::cout << "number solution" << std::endl;
	//print(u);
	long double h = .1;
	fout << "x	precise sol	numerical sol\n";
	for (int i = 0; i < 11; ++i) {
		fout << std::fixed << std::setprecision(2) << i*h << "	";
		fout << std::scientific << std::setprecision(4) << u_pr[i] << "	" << u[i] << std::endl;
	}
	fout << "---------------------------------------------" << std::endl;
	fout << "\nDYNAMIC RESULTS\n";
	//std::cout << "---------------------------------------------" << std::endl;
	//dinProgonka(11, grid, u);
	dinSolve(u, grid);
	fout << "x	sol\n";
	        for (int i = 0; i < 11; ++i) {
                fout << std::fixed << std::setprecision(2) << i*h << "	";
                fout << std::scientific << std::setprecision(4) << u[i] << std::endl;
        }
	//std::cout << "---------------------------------------------" << std::endl;
	//print(u);

	return 0;
}
