#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>

void print_vector(const std::vector <long double>& x) {
    for (auto it = x.begin(); it < x.end(); ++it)
        std::cout << *it << " ; ";
    std::cout << std::endl;
}

long double grid_step(long double a, long double b, int n) { // n - quantity of ties
    return (b - a) / (n + 1);
}

std::vector <long double> grid(long double a, long double b, long double h, long double n) {
    std::vector <long double> x;
    long double c = a + h;
    for (int i = 0; i < n; ++i) {
        x.push_back(c);
        c += h;
    }
    return x;
}

std::vector <long double> new_grid(std::vector <long double>& x) {
    int n = x.size();
    long double h = (x.at(1) - x.at(0))/2;
    long double x_st = x.at(0) - 2*h;
    long double x_fin = x.at(x.size() - 1);
    std::vector <long double> x_new;
    for (int i = 0; i < 2*n; ++i) {
        x_new.push_back(x_st + (i + 1) * h);
    }
    return x_new;
}

long double f(long double x, long double y) {
    return x * pow(y, 2) + 3 * x * y;
}

std::vector <long double> sol(long double xn, long double yn, const std::vector <long double>& x, const std::vector <long double>& x_canon) {
    long double fn = f(xn, yn);
    long double fn1;
    long double yn1;
    long double h = x.at(0) - xn;
    std::vector <long double> sol;
    int j = 0;
    for (int i = 0; i < x.size() ; ++i) {
        fn1 = f(xn + h/2, yn + fn*h/2);
        yn1 = yn + h*fn1;
        fn = fn1;
        xn = x[i];
        if (fabsl(x[i] - x_canon[j]) < pow(10, -6)) {
            sol.push_back(yn1);
            ++j;
        }
        yn = yn1;
    }
    return sol;
}

int main() {
    std::cout.setf(std::ios::fixed);
    std::cout.precision(6);
    std::ofstream fout("res.txt");
    fout.setf(std::ios::fixed);
    fout.precision(6);
    long double a = 0; //integration segment
    long double b = .6;
    long double y0 = 3; //initial condition
    long double eps = pow(10, -4); //accuracy
    int k = 1; //order of accuracy
    long double h = grid_step(a, b, 11);
    std::vector <long double> x_canon = grid(a, b, h, 11);
    std::vector <long double> x = x_canon;
    std::vector <long double> yh;
    std::vector <long double> y2h;
    std::vector <long double> y_help;
    std::vector <long double> div;
    long double delta = eps + 1;
    yh = sol(a, y0, x, x_canon);
    while (delta > eps) {
        x = new_grid(x);
        y2h = sol(a, y0, x, x_canon);
        for (int i = 0; i < yh.size(); ++i)
            div.push_back(fabsl(y2h[i] - yh[i]));
        long double max_div = *std::max_element(div.begin(), div.end());
        delta = max_div / (pow(2, k) - 1);
        fout << x.size() << " ; "<< delta << "\n";
        div.clear();
        y_help = yh;
        yh = y2h;
    }
    fout << "X : \n";
    for (auto it = x_canon.begin(); it < x_canon.end(); ++it)
        fout << *it << " ; ";
    fout << "\nY : \n";
    for (auto it = y_help.begin(); it < y_help.end(); ++it)
        fout << *it << " ; ";
    print_vector(yh);
    std::cout << "---------" << std::endl;
    fout.close();
    return 0;
}