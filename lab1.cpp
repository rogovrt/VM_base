#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

void print(const std::string& s, const double& v) {
    std::cout << s << " = " << v << std::endl;
}

void print_vector(const std::vector <long double>& v) {
    for (int i =0 ; i < v.size(); ++i)
        std::cout << v[i] << " ";
    std::cout << std::endl;
}

void print_equation(const std::vector <long double>& c) {
    int n = c.size();
    int j = 0;
    for (int i = n - 1; i >  -1; --i) {
        if (c[j] != 0) {
            std::cout << "(" << c[j] << ")" << " * Z^" << i;
            if (i != 0)
                std::cout << " + ";
            else
                std::cout << " = 0;" << std::endl;
        }
        ++j;
    }

}

int quantity_of_sign_change(const std::vector <long double>& c) {
    long double num = *c.begin();
    int q = 0;
    for (auto it = c.begin() + 1; it != c.end(); ++it) {
        if ( ((num > 0) && (*it < 0)) || ((num < 0) && (*it > 0)) ) {
            ++q;
            num = *it;
        }
    }
    return q;
}

long double result(const std::vector <long double>& c, long double x) { //result of stand in a polynomial function
    int n = c.size();
    long double res = 0;
    for (int i = 0; i < n; ++i) {
        res += c[i] * pow(x, n - 1 - i);
    }
    return res;
}

std::vector <long double> derivative(const std::vector <long double>& c) {
    std::vector <long double> a;
    int n = c.size();
    long double k;
    for (int i = 0; i < (n - 1); ++i) {
        k = c[i] * (n - 1 - i);
        a.push_back(k);
    }
    return a;
}

int delta(const std::vector <long double>& c, long double a, long double b) {
    std::vector <long double> f_a;
    std::vector <long double> f_b;
    f_a.push_back(result(c,a));
    f_b.push_back(result(c,b));
    std::vector <long double> k = c;
    for (int i = 1; i < c.size(); ++i) {
        k = derivative(k);
        f_a.push_back(result(k, a));
        f_b.push_back(result(k, b));
    }
    return quantity_of_sign_change(f_a) - quantity_of_sign_change(f_b);
}

int main() {
    long double g0 = 5./3;
    long double rho0 = 1.694 * pow(10., -4.); //gr/cm3
    long double p0 = 1.013 * pow(10., 6.); //din/cm2
    long double u0 = 0.;
    long double g3 = 7./5;
    long double c3 = 3.6537 * pow(10., 4.); //cm/c
    long double p3 = 1.6738 * pow(10., 6.); //din/cm2
    long double u3 = 0.;

    long double rho3 = g3  * p3 / pow(c3, 2.);
    long double a0 = (g0 + 1)/(g0 - 1);
    long double n = 2 * g3 / (g3 - 1);
    long double mu = (u3 - u0) * sqrt((g0 - 1) * rho0 / (2 * p0));
    long double nu = 2 / (g3 - 1) * sqrt(g3 * (g0 - 1) * p3 * rho0 / (2 * p0 * rho3));
    long double x = p3 / p0;
    //double z = pow(p1 / p3, 1/n);

    std::cout.setf(std::ios::fixed);
    std::cout.precision(8);/*
    print("alpha0", a0);
    print("n", n);
    print("mu", mu);
    print("nu", nu);
    print("x", x);
*/
    long double coef1 = pow(x, 2);
    long double coef2 = - a0 * pow(nu, 2) * x;
    long double coef3 = 2 * a0 * nu * (mu + nu) * x;
    long double coef4 = - (2 + pow((nu + mu), 2) * a0) * x;
    long double coef5 = - pow(nu, 2);
    long double coef6 = 2  * nu * (mu + nu);
    long double coef7 = - pow((mu + nu), 2);
    /*
    std::cout << coef1 << " * Z^(" << 2*n << ") ";
    std::cout << coef2 << " * Z^(" << n + 2 << ") +";
    std::cout << coef3 << " * Z^(" << n + 1 << ") ";
    std::cout << coef4 << " * Z^(" << n << ") ";
    std::cout << coef5 << " * Z^(2) +";
    std::cout << coef6 << " * Z ";
    std::cout << coef7 << " + 1 = 0" << std::endl;
    */
    std::vector <long double> coefs {coef1, 0, 0, 0, 0, coef2, coef3, coef4, 0, 0, 0, 0, coef5, coef6, coef7 + 1};
    print_equation(coefs);

    long double b = *std::max_element(coefs.begin(), coefs.end() - 2);
    long double a = *std::max_element(coefs.begin() + 1, coefs.end() - 1);
//    std::cout << b << " " << fabsl(*(coefs.end() - 1)) << std::endl;
//    std::cout << a << " " << fabsl(*(coefs.begin())) << std::endl;
    long double st = fabsl(*(coefs.end() - 1)) / (fabsl(*(coefs.end() - 1)) + b);
    long double fin = 1 + a / fabsl(*(coefs.begin()));
    std::cout << "solution locates from " << st << " to " << fin << std::endl;

    int q = quantity_of_sign_change(coefs); // quantity of positive solutions
    std::cout << "quantity of positive solutions = " << q << std::endl;
/*
    std::cout << result(coefs, 1.2);
    std::vector <long double> d = derivative(coefs);
    print_equation(d);
    std::cout << result(d, 1.2);
*/

    std::cout << delta(coefs, st, fin) << std::endl;
    std::cout << delta(coefs, st, 1.5) << std::endl;
    std::cout << delta(coefs, st, .9) << std::endl;
    std::cout << delta(coefs, .9, 1.) << std::endl; //root!
    std::cout << delta(coefs, st, .8) << std::endl;
    std::cout << delta(coefs, .8, .9) << std::endl; //root!
    std::cout << delta(coefs, st, .3) << std::endl; //root!
    std::cout << delta(coefs, st, .23) << std::endl; //root!

    return 0;
}