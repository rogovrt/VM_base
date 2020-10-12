#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>

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

std::vector <long double> derivative(const std::vector <long double>& c) { //derivative of polynomial function
    std::vector <long double> a;
    int n = c.size();
    long double k;
    for (int i = 0; i < (n - 1); ++i) {
        k = c[i] * (n - 1 - i);
        a.push_back(k);
    }
    return a;
}

int delta(const std::vector <long double>& c, long double a, long double b) { //using theoren budan-furie
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


std::pair <long double, long double> localization(const std::vector <long double>& c, long double a, long double b) {
    long double r = fabsl(b - a);
    std::pair <long double, long double> res;
    while (r > .000001) {
        if (delta(c, a + r/2, b) == 1)
            a = a + r / 2;
        else
            b = b - r / 2;
        r = fabsl(b - a);
    }
    res = std::make_pair(a, b);
    return res;
}

long double find_root_newton(const std::vector <long double>& c, std::pair <long double, long double> segm) {
    long double a = segm.first;
    long double b = segm.second;
    std::vector <long double> der1 = derivative(c);
    std::vector <long double> der2 = derivative(der1);

    long double x0 = a + fabsl(b - a) / 2;
    if (result(c, x0) * result(der2, x0) < 0)
        return 0;
    long double k = result(der2, a) / (2 * result(der1, a));
    long double x1 = x0 - result(c, x0)/ result(der1, x0);
    while ((k * pow(x1 - x0, 2)) > .000000001) {
        x0 = x1;
        x1 = x0 - result(c, x0)/ result(der1, x0);
    }
    return x1;
}

long double find_root(const std::vector <long double>& c, std::pair <long double, long double> segm) {
    long double a = segm.first;
    long double b = segm.second;
    long double cd;
    if (result(c, a) * result(c, b) >= 0) std::cout << "problem!" << std::endl;
    while (b - a > .000000001) {
        cd = (b + a) / 2;
        if (result(c, a) * result(c, cd) < 0) b = cd;
        else a = cd;
    }
    return (b + a) / 2;
}

void roots(const std::vector <long double>& c, long double a, long double b, std::vector<long double>& sol) {
    long double r = fabsl(b - a);
    long double q;
        if (delta(c, a, a + r / 2) == 1) {
            q = find_root(c, localization(c, a, a + r / 2));
            sol.push_back(q);
        }
        if ((delta(c, a, a + r / 2) > 1) && (a != a + r/2)) {
            roots(c, a, a + r /4, sol);
            roots(c, a + r / 4, a + r / 2, sol);
        }
        if (delta(c, a + r / 2, b) == 1) {
            q = find_root(c, localization(c, a + r / 2, b));
            sol.push_back(q);
        }
        if ((delta(c, a + r / 2, b) > 1) && (a != a + r/2)) {
            roots(c, a + r / 2, a + 3 * r /4, sol);
            roots(c, a + 3 * r / 4, b, sol);
        }
}
int main() {
    long double g0 = 5./3;
    long double rho0 = 1.694 * pow(10., -4.); //gr/cm3
    long double p0 = 1.013 * pow(10., 6.); //din/cm2
    long double u0 = 0.;
    long double g3 = 7./5;
    long double c3 = 3.6537 * pow(10., 4.); //cm/c
    long double p3 = 1.6768 * pow(10., 6.); //din/cm2
    long double u3 = 0;

    long double rho3 = g3  * p3 / pow(c3, 2.);
    long double a0 = (g0 + 1)/(g0 - 1);
    long double n = 2 * g3 / (g3 - 1);
    long double mu = (u3 - u0) * sqrt((g0 - 1) * rho0 / (2 * p0));
    long double nu = 2 / (g3 - 1) * sqrt(g3 * (g0 - 1) * p3 * rho0 / (2 * p0 * rho3));
    long double x = p3 / p0;
    //double z = pow(p1 / p3, 1/n);

    std::cout.setf(std::ios::fixed);
    std::cout.precision(8);

    long double coef1 = pow(x, 2);
    long double coef2 = - a0 * pow(nu, 2) * x;
    long double coef3 = 2 * a0 * nu * (mu + nu) * x;
    long double coef4 = - (2 + pow((nu + mu), 2) * a0) * x;
    long double coef5 = - pow(nu, 2);
    long double coef6 = 2  * nu * (nu + mu);
    long double coef7 = - pow((mu + nu), 2);

    std::vector <long double> coefs {coef1, 0, 0, 0, 0, coef2, coef3, coef4, 0, 0, 0, 0, coef5, coef6, coef7 + 1};
    print_equation(coefs);

    long double b = *std::max_element(coefs.begin(), coefs.end() - 2); //
    long double a = *std::max_element(coefs.begin() + 1, coefs.end() - 1);
    long double st = fabsl(*(coefs.end() - 1)) / (fabsl(*(coefs.end() - 1)) + b);
    long double fin = 1 + a / fabsl(*(coefs.begin()));
    std::cout << "solution locates from " << st << " to " << fin << std::endl;

    int q = quantity_of_sign_change(coefs); // quantity of positive solutions
    std::cout << "quantity of positive solutions = " << q << std::endl;

    std::vector <long double> r;
    roots(coefs, st, fin, r);

    for (auto it = r.begin(); it < r.end(); ++it)
       std::cout << *it << std::endl;

    return 0;
}