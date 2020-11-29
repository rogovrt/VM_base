#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

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

long double f(long double x, long double y) {
    return x * pow(y, 2) + 3 * x * y;
}

std::vector <long double> solution(long double xn, long double yn, long double h, int part, const std::vector <long double>& x) { //a, y0
    long double fn = f(xn, yn);
    long double fn1;
    long double yn1;
    long double h1 = h / part;
    int q = 0;
    std::vector <long double> sol;
    while(xn < x.at(x.size() - 1)) {
        fn1 = f(xn + h1 / 2, yn + fn*h1 / 2);
        yn1 = yn + h1 * fn1;
        fn = fn1;
        xn += h1;
        q++;
        if (q == part) {
            sol.push_back(yn1);
            q = 0;
        }
        yn = yn1;
    }
    return sol;
}

int main() {
    std::cout.setf(std::ios::fixed);
    std::cout.precision(8);
    long double a = 1; //integration segment
    long double b = 2;
    long double y0 = .01; //initial condition
    long double eps = pow(10, -4); //accuracy
    int k = 1; //order of accuracy
    long double h = grid_step(a, b, 11);
    std::vector <long double> x = grid(a, b, h, 11);
    //print_vector(x);
    std::vector <long double> yh;
    std::vector <long double> y2h;
    std::vector <long double> div;
    long double delta = eps + 1;
    int part = 2;
    while (delta > eps) {
        yh = solution(a, y0, h, part, x);
        y2h = solution(a, y0, h, part * 2, x);
        //print_vector(yh);
        //print_vector(y2h);
        for (int i = 0; i < yh.size(); ++i)
            div.push_back(fabsl(y2h[i] - yh[i]));
        //print_vector(div);
        long double max_div = *std::max_element(div.begin(), div.end());
        delta = max_div / (pow(2, k) - 1);
        //std::cout << delta << std::endl;
        part += 1;
        div.clear();
    }
    print_vector(x);
    print_vector(yh);
    return 0;
}