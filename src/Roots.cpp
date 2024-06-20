#include "Roots.h"
using namespace std;

void Roots::AddFunction(string name, function<double(double)> funct){
    functions[name] = funct;
}

double Roots::Bisection(string name, double xi, double xf, double e){
    if (functions.find(name) == functions.end()) {
        cerr << "Function not found!" << endl;
        return 0.0; // or some error value
    }

    while((xf-xi)>e){
        double f1 = functions[name](xi);
        double f2 = functions[name]((xi+xf)/2);
        double f3 = functions[name](xf);
        if(f2*f1 <= 0.0){
            xf = (xi+xf)/2;
        }
        else{
            xi = (xi+xf)/2;
        }
    }
    return (xi+xf)/2;
}

//THE FORWARD ROOT METHODS HAVE A CONTINGENCY PLAN FOR WHEN THE USED FUNCTION HAS A DOMINIUM THAT IS CLOSED
//WHEN A nan VALUE IS FOUND WE MAKE A SMALL STEP IN THAT DIRECTION

double Roots::RegulaFalsi(string name, double xi, double xf, double e){
    if (functions.find(name) == functions.end()) {
        cerr << "Function not found!" << endl;
        return 0.0; // or some error value
    }

    while(true){
        double f1 = functions[name](xi);
        double f3 = functions[name](xf);
        //can make this in one like by making an algebric simplification of the -b/a as b is a function of a in this case
        double a = (f3-f1)/(xf-xi);
        double b = f1-a*xi;
        if(IsValid(name, -b/a)){
            double f2 = functions[name](-b/a);
            if(abs(f2) < e) return -b/a; //verify if point is found, if f2 is 0 it breaks the reducing of intervals
            if(f2*f1 <= 0.0){xf = -b/a;}
            else{ xi = -b/a;}
        }
        
        else if( -b/a > xi){xi = xi+e;}
        else if( -b/a < xi){xi = xi-e;}
    }
}

double Roots::NewtonRaphson(string name, double xi, double var, double e){
    if (functions.find(name) == functions.end()) {
        cerr << "Function not found!" << endl;
        return 0.0; // or some error value
    }

    while(true){
        if(abs(functions[name](xi)) < e){return xi;}
        double f1 = functions[name](xi);
        //can make this a method of first derivative;
        double f2 = functions[name](xi-var);
        double f3 = functions[name](xi+var);
        double a = (f3-f2)/(2*var);
        double b = f1-a*xi;
        if(IsValid(name, -b/a)){
            double f4 = functions[name](-b/a);
            xi = -b/a;
        }
        //Emergency plan, if function is non continuous let it walk in the direction of the line
        else if( -b/a > xi){xi = xi+var;}
        else if( -b/a < xi){xi = xi-var;}
    }
}

double Roots::Secant(string name, double xi, double var, double e){
    if (functions.find(name) == functions.end()) {
        cerr << "Function not found!" << endl;
        return 0.0; // or some error value
    }

    while(true){
        if(abs(functions[name](xi)) < e){return xi;}
        double f1 = functions[name](xi);
        //can make this a method of first derivative;
        double f2 = functions[name](xi-var);
        double a = (f1-f2)/(var);
        double b = f1-a*xi;
        if(IsValid(name, -b/a)){
            double f4 = functions[name](-b/a);
            xi = -b/a;
        }
        //Emergency plan, if function is non continuous let it walk in the direction of the line
        else if( -b/a > xi){xi = xi+var;}
        else if( -b/a < xi){xi = xi-var;}
    }
}


//can be used to calculate Pi using a function sqrt(1-x*x) with x[0,1] and y[0,1] -> pi = 4*return/(r*r)
vector<double> Roots::Evaluation(string name, double xi, double xf, double yi, double yf, int N){
    if (functions.find(name) == functions.end()) {
        cerr << "Function not found!" << endl;
        return {}; // or some error value
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> disX(xi, xf);
    uniform_real_distribution<> disY(yi, yf);
    int MinCount = 0;
    vector<double> values = {};
    for(int i = 1; i <= N; i++){
        double x = disX(gen);
        double y = disY(gen);
        if(abs(functions[name](x)) >= abs(y)){MinCount++;}
        values.push_back((double)(MinCount)/(double)(i));
        if(i%1000 == 0){cout << (double)(i)/(double)(N) * 100.0 << "%" << endl;}
    }
    return values;
}

bool Roots::IsValid(string name, double x) {
    double result = functions[name](x);
    return !isnan(result) && isfinite(result);
}