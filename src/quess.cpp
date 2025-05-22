#include "quess.hpp"
#include <iostream>

Quess::Quess(std::initializer_list<std::pair<std::string, double>> init) {
    values["l"] = 1.0;
    values["n"] = 4.0;
    values["c"] = 1.0;
    values["A"] = 0.15;
    values["epsilon"] = 0.1;

    for (const auto &p : init) {
        if (values.count(p.first) == 0)
            throw std::invalid_argument("Invalid parameter: " + p.first);
        values[p.first] = p.second;
    }


	//  double lambda=(4*l)/(2*n+1);
	//  double gamma=(2*M_PI)/lambda;
	//  double omega= gamma*c;

    values["lambda"] = (4 * values["l"]) / (2 * values["n"] + 1);
    values["gamma"]  = (2.0 * M_PI) / values["lambda"];
    values["omega"]  = values["gamma"] * values["c"];
}

double Quess::operator[](const std::string &key) const {
    auto it = values.find(key);
    if (it == values.end())
        throw std::out_of_range("Invalid key: " + key);
    return it->second;
}

double Quess::f(double x) const {
    double eps = values.at("epsilon");
    if (x <= eps) {
        return 1.0 + (2.0 * x * x * x) / (eps * eps * eps)
               - (3.0 * x * x) / (eps * eps);
    }
    return 0.0;
}

double Quess::mollifier(double x) const {
    double eps = values.at("epsilon");
    if (std::abs(x) >= eps) return 0.0;
    double r = x / eps;
    return std::exp(-1.0 / (1.0 - r * r)) * std::exp(1.0);
}

void Quess::print() const {
    std::cout << "Quess parameters:\n";
	std::cout << "  l        = " << values.at("l")        << "\n";
    std::cout << "  n        = " << values.at("n")        << "\n";
    std::cout << "  c        = " << values.at("c")        << "\n";
	std::cout << "  omega        = " << values.at("omega")        << "\n";
	std::cout << "  lambda        = " << values.at("lambda")        << "\n";
	std::cout << "  gamma        = " << values.at("gamma")        << "\n";
    std::cout << "  A        = " << values.at("A")        << "\n";
    std::cout << "  epsilon  = " << values.at("epsilon")  << "\n";
    
}