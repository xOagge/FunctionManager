#ifndef __ROOTS__
#define __ROOTS__

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <random>
#include <functional>
using namespace std;

class Roots{
    public:
        Roots() = default;
        ~Roots() = default;

        void AddFunction(string name, function<double(double)> funct);
        //use when there is only one 0
        double Bisection(string name, double xi, double xf, double e = 1e-7);
        double RegulaFalsi(string name, double xi, double xf, double e = 1e-7);
        double NewtonRaphson(string name, double xi, double var, double e = 1e-7);
        double Secant(string name, double xi, double var, double e = 1e-7);


        //MONTE-CARLO
        vector<double> Evaluation(string name, double xi, double xf, double yi, double yf, int N);


        bool IsValid(string name, double x);

    private:
        map<string, function<double(double)>> functions; //{F_name : F}
};

#endif