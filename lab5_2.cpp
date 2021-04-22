#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <utility>

std::ofstream fout("out.txt");

template <typename T>
void printVector(T& v) {
	for (auto it = v.begin(); it < v.end(); ++it)
		std::cout << *it << " ; ";
	std::cout << std::endl;
}

void print2D(std::vector<std::vector<long double>>& u) {
	for (int i = 0; i < u.size(); ++i) {
		for (int j = 0; j < u[i].size(); ++j)
			fout << std::fixed << std::scientific << std::setprecision(7) << u[i][j] << " ; ";
		fout << std::endl;
	}	
}

void generateTn(std::vector<long double>& t_res, int n_it, long double min, long double max) {
	std::vector<int> n;
	int i = n_it;
	n.push_back(i);
	while (i != 1) {
		if ((i % 2) == 0) i = i / 2;
		else --i;
		n.push_back(i);
	}
	std::reverse(n.begin(), n.end());
	
	std::vector<int> t;
	std::vector<int> t_next;
	t_next.push_back(1);
	i = 1;
	int flag = 0;
	while (i < n.size()) {
		t = t_next;
		t_next.resize(n[i]);
		//printVector(t);
		if (flag == 1) {
			for (int j = 0; j < n[i]-1; ++j)
				t_next[j] = t[j];
			t_next[n[i]-1] = n[i];
			flag = 0;
			++i;
		}
		else { //
		if (2*n[i] == n[i+1]) {
			for (int j = 0; j < n[i]; ++j) {
				if ((j % 2) == 0) t_next[j] = t[j/2];
				else t_next[j] = 4*n[i-1] - t_next[j-1];
			}
			++i;
		}
		else {
			for (int j = 0; j < n[i]; ++j) {
				if ((j % 2) == 0) t_next[j] = t[j/2];
				else t_next[j] = 4*n[i-1] + 2 - t_next[j-1];
			}
			++i;
			flag = 1;
		}
		}		
	}
	//s = t_next;
	
	for (int i = 0; i < t_next.size(); ++i) {
		t_next[i] = (t_next[i] + 1) / 2;
	}
	
	long double tn;
	for (int i = 0; i < n_it; ++i) {
		tn = 2 / (max + min + (max - min)*cos(M_PI*(2*t_next[i]-1)/(2*n_it)));
		t_res[i] = tn;
	}
}

void nextPart(int M, int L, std::vector<std::vector<long double>>& u, long double t) {
	long double x;
	long double y;
	long double h_y = 1. / M;
        long double h_x = 1. / L;
	std::vector<std::vector<long double>> u_next(M+1);
	std::vector<long double> c(L+1);
	
	for (int i = 0; i < M+1; ++i) {
		y = h_y*i;
		if ((i == 0) || (i == M)) {
			for (int j = 0; j < L+1; ++j)
				c[j] = 0.;
		}
		else {
		for (int j = 0; j < L+1; ++j) {
			x = -.5 + h_x*j;
			if ((j == 0) || (j == L)) c[j] = 0.;
			else
			c[j] = t/pow(h_x, 2)*(u[i][j+1] - 2*u[i][j]+u[i][j-1]) + t/pow(h_y, 2)*(u[i+1][j] - 2*u[i][j] + u[i-1][j]) - 2*y*(1-y)*t-.5*(1-4*pow(x,2))*t + u[i][j];
		}
		}
		u_next[i] = c;
		//c.clear();
		//c.resize(L+1);
	}
	
	/*
	for (int i = 0; i < M+1; ++i) { 
		u_next[i][0] = 0.;
		u_next[i][L] = 0.;
	}
	for (int j = 0; j < L+1; ++j) {
		u_next[0][j] = 0.;
		u_next[M][j] = 0.;
	}
	*/
	u = u_next;
	//print2D(u);
}

void solve(std::vector<std::vector<long double>>& solution, int L, int M) {
	//long double min_nu = 2*pow(M_PI, 2);
	//long double max_nu = 4*(pow(L, 2) + pow(M, 2));
	long double h = 1. / M;
	long double min_nu = 8*pow(sin(M_PI/(2*L)) / h, 2);
	long double max_nu = 8*pow(cos(M_PI/(2*L)) / h, 2);
	//std::cout << min_nu << " ; " << max_nu << std::endl; 
	int N = ceil(log(2./pow(10., -6)) / log((sqrt(max_nu)+sqrt(min_nu)) / (sqrt(max_nu)-sqrt(min_nu))));
	fout << N << "	| ";
	std::vector<long double> t(N);
	generateTn(t, N, min_nu, max_nu);
	//printVector(t);
	// init start approximation of solution
	std::vector<std::vector<long double>> u(M+1);
	std::vector<long double> c(L+1);
        for (int j = 0; j < L+1; ++j)
        	c[j] = 0.;
	for (int i = 0; i < M+1; ++i) {
		u[i] = c;
	}
        for (int i = 0; i < M+1; ++i) { 
                u[i][0] = 0.;
                u[i][L] = 0.;
        }
        for (int j = 0; j < L+1; ++j) {
                u[0][j] = 0.;
                u[M][j] = 0.;
        }
	
	for (int k = 0; k < N; ++k) {
		nextPart(M, L, u, t[k]);
	}
	
	int step = L/5;
	for (int i = 0; i < L+1; ++i) {
		if ((i % step) == 0) {
			for (int j = 0; j < L+1; ++j) {
				if ((j % step) == 0) solution[i/step][j/step] = u[i][j];
			}
		}
	}
	
}

void generatePrecise(std::vector<std::vector<long double>>& u) {
	long double x;
	long double y;
	std::vector<long double> c(6);
	for (int i = 0; i < 6; ++i) {
		y = i * .2;
		for (int j = 0; j < 6; ++j) {
			x = -.5 + .2*j;
			c[j] = .25 * y * (y - 1) * (1 - 4 * pow(x, 2)); 
		}
		u[i] = c;
	}
}

long double delta(const std::vector<std::vector<long double>>& u, const std::vector<std::vector<long double>>& u_pr) {
	long double e = -1.;
	long double e1;
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 6; ++j) {
			e1 = fabsl ((u[i][j] - u_pr[i][j]) / u_pr[i][j]);
			if (e1 > e) e = e1;
		}
	}
	return e;
}

int main() {
	std::vector<std::vector<long double>> u(6);
	std::vector<long double> c(6);
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 6; ++j)
			c[j] = 0.;
		u[i] = c;
	}
	std::vector<std::vector<long double>> u_pr(6);
	generatePrecise(u_pr);
	int L;
	fout << "N	| L=M	| du/u\n";
	for (int i = 1; i < 10; ++i) {	
		L = 5*i;
		solve(u, L, L);
		fout << L << "	|" << delta(u, u_pr) << std::endl;
	}
	fout << std::endl << "Последнее полученное решение : \n";
	print2D(u);
	fout << std::endl << "Аналитическое решение : \n";
	print2D(u_pr);
	return 0;
}
