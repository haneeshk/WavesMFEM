
#ifndef QUESS_HPP
#define QUESS_HPP

#include <string>
#include <map>
#include <initializer_list>
#include <stdexcept>
#include <cmath>

class Quess {
public:
    Quess(std::initializer_list<std::pair<std::string, double>> init);
    double operator[](const std::string &key) const;
    double f(double x) const;
    double mollifier(double x) const;
	void print() const;  // ← add this

private:
    std::map<std::string, double> values;
};

#endif // QUESS_HPP