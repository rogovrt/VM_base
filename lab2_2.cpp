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
long double kl = 1.;
long double kr = 1.;
long double ql = exp(-0.5);
long double qr = 1.;
long double fl = cos(x0);
long double fr = 1.;

long double k(long double x) {
	return pow(x, 2.) + .5;
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

long double findDelta(std::array<long double, 11>& y, const std::array<long double, 11>& y_pr) {
	long double delta = fabsl(y[0] - y_pr[0]);
	for (int i = 1; i < 11; ++i) {
		if (fabsl(y[i] - y_pr[i]) > delta)
			delta = fabsl(y[i] - y_pr[i]);
	}
	return delta;
}

void statPrecise(std::array<long double, 11>& u) {
	long double ll = sqrt(ql/kl);
	long double lr = sqrt(qr/kr);
	long double nl = fl/ql;
	long double nr = fr/qr;
	long double A11 = exp(-ll * x0) - exp(ll * x0);
	long double A12 = exp(lr * (2 - x0)) - exp(lr * x0);
	long double A21 = kl * ll * (exp(ll*x0) + exp(-ll*x0));
	long double A22 = kr*lr*(exp(lr*(2-x0)) + exp(lr*x0));
	long double B1 = nr - nl + (nl - u0)*exp(ll*x0) - (nr - u1) * exp(lr * (1 - x0));
	long double B2 = kl * ll * (u0 - nl) * exp(ll * x0) + kr * lr * (u1 - nr) * exp(lr * (1 - x0));
	long double C1 = (((u0 - nl) * A11 - B1) * A22 - ((u0 - nl) * A21 - B2) * A12) / (A11 * A22 - A12 * A21);
	long double C2 = (B1*A22 - B2*A12) / (A11*A22 - A12*A21);
	long double C3 = (B2*A11 - B1*A21) / (A11*A22 - A12*A21);
	long double C4 = (u1 - nr)*exp(lr) - C3*exp(2*lr);
	
	long double x  = 0.;
        for (int i = 0; i < 11; ++i) {
                x = i * .1;
		if (x < x0)
                	u[i] = C1*exp(ll*x) + C2*exp(-ll*x) + nl;
		else
			u[i] = C3*exp(lr*x) + C4*exp(-lr*x) + nr;
        }

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
void statLeftCoeffs(T& al, T& bl, T& cl, T& dl, T& grid) {
        long double h = grid[1] - grid[0];
        for (int i = 0; i < al.size(); ++i) {
                al[i] = kl;
                bl[i] = -(2*kl + ql*pow(h, 2.));
                cl[i] = kl;
                dl[i] = -fl * pow(h, 2.);
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

template <typename T>
void statRightCoeffs(T& ar, T& br, T& cr, T& dr, T& grid, int beta) {
        long double h = grid[1] - grid[0];
        for (int i = 0; i < ar.size(); ++i) {
                ar[i] = kr;
                br[i] = -(2*kr + qr*pow(h, 2.));
                cr[i] = kr;
                dr[i] = -fr * pow(h, 2.);
        }
}

void dinProgonka(int numberOfKnots, std::vector<long double>& grid, std::array<long double, 11>& u) {
	int l = (numberOfKnots - 1)/10;
        int alpha; int beta;
        createGrid(numberOfKnots, grid, alpha, beta);
        std::vector<long double> u_(numberOfKnots);
        u_[0] = u0;     u_[numberOfKnots-1] = u1;
        
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
	bll[0] = (dl[0] - cl[0]*u0) / bl[0];
	for (int i = 1; i < al.size(); ++i) {
		all[i] = -al[i] / (bl[i] + cl[i]*all[i-1]);
		bll[i] = (dl[i] - cl[i]*bll[i-1]) / (bl[i] + cl[i]*all[i-1]);
	}
	arr[ar.size()-1] = -cr[ar.size()-1] / br[ar.size()-1];
	brr[ar.size()-1] = (dr[ar.size()-1] - cr[ar.size()-1]*u1) / br[ar.size()-1];
	for (int i = ar.size() - 2; i > -1; --i) {
		arr[i] = -cr[i] / (br[i] + ar[i]*arr[i+1]);
		brr[i] = (dr[i] - ar[i]*brr[i+1]) / (br[i] + ar[i]*arr[i+1]);
	}


        long double temp = (k(grid[alpha]) * bll[bll.size()-1] + k(grid[beta])  * brr[0]) / (k(grid[alpha]) * (1 - all[bll.size()-1]) + k(grid[beta]) * (1 - arr[0]));
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

void statProgonka(int numberOfKnots, std::vector<long double>& grid, std::array<long double, 11>& u) {
        int l = (numberOfKnots - 1)/10;
        int alpha; int beta;
        createGrid(numberOfKnots, grid, alpha, beta);
        std::vector<long double> u_(numberOfKnots);
        u_[0] = u0;     u_[numberOfKnots-1] = u1;

        std::vector<long double> al(alpha - 1);
        std::vector<long double> bl(alpha - 1);
        std::vector<long double> cl(alpha - 1);
        std::vector<long double> dl(alpha - 1);
        std::vector<long double> ar(numberOfKnots - beta - 2);
        std::vector<long double> br(numberOfKnots - beta - 2);
        std::vector<long double> cr(numberOfKnots - beta - 2);
        std::vector<long double> dr(numberOfKnots - beta - 2);
        statLeftCoeffs(al, bl, cl, dl, grid);
        statRightCoeffs(ar, br, cr, dr, grid, beta);

        std::vector<long double> all(alpha - 1);
        std::vector<long double> bll(alpha - 1);
        std::vector<long double> arr(numberOfKnots - beta - 2);
        std::vector<long double> brr(numberOfKnots - beta - 2);
        all[0] = -al[0]/bl[0];
        bll[0] = (dl[0] - cl[0]*u0) / bl[0];
        for (int i = 1; i < al.size(); ++i) {
                all[i] = -al[i] / (bl[i] + cl[i]*all[i-1]);
                bll[i] = (dl[i] - cl[i]*bll[i-1]) / (bl[i] + cl[i]*all[i-1]);
        }
        arr[ar.size()-1] = -cr[ar.size()-1] / br[ar.size()-1];
        brr[ar.size()-1] = (dr[ar.size()-1] - cr[ar.size()-1]*u1) / br[ar.size()-1];
        for (int i = ar.size() - 2; i > -1; --i) {
                arr[i] = -cr[i] / (br[i] + ar[i]*arr[i+1]);
                brr[i] = (dr[i] - ar[i]*brr[i+1]) / (br[i] + ar[i]*arr[i+1]);
        }

	
        long double temp = (kl * bll[bll.size()-1] + kr  * brr[0]) / (kl * (1 - all[bll.size()-1]) + kr * (1 - arr[0]));
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
        //while (n < 50000) {
		n = n1;
		u = u1;
                n1 = (n1 - 1)*2 + 1;
                dinProgonka(n1, grid, u1);
		fout << n << "	" << std::scientific << std::setprecision(4) << findDelta(u, u1) << std::endl;
        }
	fout << std::endl;
}

void statSolve(std::array<long double, 11>& u, std::vector<long double>& grid) {
	std::array<long double, 11> u_pr;
	statPrecise(u_pr);
	std::cout << "precise solution" << std::endl;
	int n = 11;
	statProgonka(n, grid, u);
	fout << "knots	delta\n";
	fout << n << "	";
	fout << std::fixed << std::scientific << std::setprecision(4) << findDelta(u, u_pr) << std::endl;
	int k = 1;
	while (fabsl(findDelta(u, u_pr) > .0001)) {
	//while (n < 50000) {
		n = (n - 1)*2 + 1;
		statProgonka(n, grid, u);
		fout << n << "	" << std::scientific << std::setprecision(4) << findDelta(u, u_pr) << std::endl;
	}
	fout << std::endl;
}

int main() {
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
