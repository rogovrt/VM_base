#include <iostream>
#include <cmath>
#include <vector>

std::vector <long double> next_order(const std::vector <long double>& previous_order, const std::vector <long double>& x) {
    int s1 = previous_order.size();
    int s2 = x.size();
    int h = s2 - s1 + 1;
    std::vector <long double> next_order;
    for (int i = 0; i < (s1 - 1); ++i) {
        next_order.push_back((previous_order[i] - previous_order[i+1])/(x[i] - x[i+h]));
    }
    return next_order;
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

std::vector <long double> coefs_splain_i(long double xi, long double xi1, long double fi, long double fi1,
                                        long double di, long double di1) {
    std::vector <long double> c;
    long double a3 = (di1*(xi1 - xi) - 2*(fi1 - fi) + di*(xi1 - xi)) / pow((xi1 - xi), 3);
    long double a2 = (-di1*(xi1 - xi)*(xi1 + 2*xi) + 3*(fi1 - fi)*(xi1 + xi) - di*(xi1 - xi)*(xi+2*xi1)) / pow((xi1 - xi), 3);
    long double a1 = (di1*xi*(2*xi1 + xi)*(xi1 - xi) - 6*(fi1 - fi)*xi*xi1 + di*xi1*(xi1 + 2*xi)*(xi1 - xi)) / pow((xi1 - xi), 3);
    long double a0 = (-di1*xi*xi*xi1*(xi1 - xi) + fi1*xi*xi*(3*xi1-xi) + fi*xi1*xi1*(xi1 - 3*xi) - di*xi*xi1*xi1*(xi1-xi)) / pow((xi1 - xi), 3);
    c.push_back(a3);
    c.push_back(a2);
    c.push_back(a1);
    c.push_back(a0);
    return c;
}

int main() {
    //std::cout.setf(std::ios::fixed);
    //std::cout.precision(8);

    std::vector<long double> y{4. * pow(10., -5), 0.00068, 0.00518, 0.02554, 0.09624, 0.30046, 0.81548};
    std::vector<long double> x{0.52360, 0.87267, 1.22173, 1.57080, 1.91986, 2.26893, 2.61799};

    //std::vector<long double> x{-3., -2., -1., 0., 1., 2., 3.};
    //std::vector<long double> y {-182., -49., -6., 1., 2., 3., -14.};

    std::vector<long double> c;
    std::vector<long double> no;
    std::vector<long double> po = y;
    c.push_back(y[0]);
    for (int i = 1; i < 7; ++i) {
        no = next_order(po, x);
        c.push_back(no[0]);
        po = no;
    }
    for (int i = 0; i < c.size(); ++i) {
        if (c[i] > 0) std::cout << "+";
        //if (c[i] < 0) std::cout << "-";
        if (c[i] == 0) break;
        std::cout << c[i];
        for (int j = 7; j > 7 - i; --j)
            std::cout << "(x-" << x[7 - j] << ")";
    }
    std::cout << std::endl;

    std::vector<long double> coefs = {0.01475509, -0.08146606, 0.22139746, -0.32930539, 0.27495844, -0.11958351,
                                      0.02080483};

    //std::vector<long double> coefs = {-1, 3, -2, 1, 1};
    std::vector <long double> der_coefs = derivative(coefs);

    std::vector<long double> ci;
    std::vector <std::vector<long double>> splains(x.size()-1);
    for (int i = 0; i < x.size()-1; ++i) {
        ci = coefs_splain_i(x[i], x[i+1], y[i], y[i+1], result(der_coefs, x[i]), result(der_coefs, x[i+1]));
        splains.at(i) = ci;
        std::cout << i << "->" << i+1 << " : ";
        for (int j = 0; j < ci.size(); ++j) std::cout << ci[j] << "; ";
        std::cout << std::endl;
    }
    /*for (int i = 0; i < x.size()-1; ++i)*/

    long double x0 = 0.95;
    int n = -1;
    for (int i = 0; i < x.size() - 1; ++i) {
        if ((x0 > x[i]) && (x[0] < x[i+1]))
            n = i;
    }
    if (n > -1)
        std::cout << "y(x0) = " << result(splains[n], x0) << " ; x0 = " << x0 << std::endl;
    std::cout << result(coefs, x0) << std::endl;
    return 0;
}