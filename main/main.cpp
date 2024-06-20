#include <vector>
#include "PlotsMaker.h"
#include "Roots.h"
#include <cmath>
#include <chrono>
#include <algorithm>
using namespace std;

int main(){
    Roots rootFinder;
    rootFinder.AddFunction("quadratic", [](double x){return 2*x*x+7*x-12;});
    rootFinder.AddFunction("circle", [](double x){return sqrt(1-x*x);});

    /*Both these methods will work fine as logn as the entry xi and xf represent a closed domain, as these methods work by choosing smaller inner domains*/
    //works fine as long as there is no discontinuity
    double rootBisection = rootFinder.Bisection("circle", 0, 100);
    cout << "1111\n";
    //works fine as long as there is no discontinuity
    double rootRegulaFalsi = rootFinder.RegulaFalsi("circle", 0, 1);
    cout << "1111\n";
    //works fine, both have a contingency plan used for when the fucntion has a border domain.
    /*
    Continuous domain -> more precision = -(3rt(betther derivative), 4th(less tolerance))
    Border Domain-> the Y= 0 of the tangent can be out of domain, in this case we swith the iteration to step-Xchange. In order to evade infinte loop, if we -4rth
    we need to -3rd, as in contingency the size of the steps is 3rd, and in casses such as "circle" where the derivative goes to inf in x=1 we need to consider small
    steps for a certain tolerance.
    */
    double rootNewtonRaphson = rootFinder.NewtonRaphson("circle",0.3, 1e-4, 0.10);
    cout << "1111\n";
    /*
    Continuous domain -> more precision = -(3rt(betther derivative), 4th(less tolerance))
    Border Domain-> the Y= 0 of the tangent can be out of domain, in this case we swith the iteration to step-Xchange. In order to evade infinte loop, if we -4rth
    we need to -3rd, as in contingency the size of the steps is 3rd, and in casses such as "circle" where the derivative goes to inf in x=1 we need to consider small
    steps for a certain tolerance.
    */
    double rootSecant = rootFinder.Secant("circle",0.3, 1e-4, 0.10);

    // Output results
    cout << "Bisection Method:" << endl;
    cout << "Root found: " << rootBisection << endl;

    cout << "\nRegula Falsi Method:" << endl;
    cout << "Root found: " << rootRegulaFalsi << endl;

    cout << "\nNewtonRaphson Method:" << endl;
    cout << "Root found: " << rootNewtonRaphson << endl;

    cout << "\nSecant Method:" << endl;
    cout << "Root found: " << rootSecant << endl;
    

    // Obtain Probs from rootFinder.Evaluation("linear", 0, 1, 0, 1, 10000);
    vector<double> Probs = rootFinder.Evaluation("circle", 0, 1, 0, 1, 100000);

    // Transform Probs into PiValues
    vector<double> PiValues(Probs.size());
    vector<double> x(Probs.size());

    for (int i = 0; i < Probs.size(); ++i) {
        PiValues[i] = Probs[i] * 4;  // Adjust transformation as needed
        x[i] = i;  // Assuming x[i] represents some index or iteration count
    }

    // Print sizes for verification
    cout << "Probs.size() = " << Probs.size() << endl;
    cout << "PiValues.size() = " << PiValues.size() << endl;
    cout << "x.size() = " << x.size() << endl;

    // Create PlotsMaker instance and make histogram
    PlotsMaker P;
    P.MakeHistogram(PiValues, "PiValues.pdf", "Pi", "Number of Iterations", "Approximated value");

    return 0;
}