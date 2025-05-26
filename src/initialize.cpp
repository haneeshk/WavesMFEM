#include <nlohmann/json.hpp>
#include <fstream>
#include <cmath>  // if you want to use math functions like sin, cos

using json = nlohmann::json;

json initialize_config() {
    // Build the JSON structure programmatically
    double L = 1.0;
    int n = 4;
    double lambda = 2 * L / n;
    double epsilon = 0.1;

    json config;

    config["testName"] = "test11";
    config["Physical Parameters"] = {
        {"Bar Length", L},
        {"n", n},
        {"lambda", lambda},
        {"epsilon", epsilon}
    };

    config["initial conditions"] = {
        {"A", 0.15},
        {"epsilon", epsilon}
    };

    return config;
}
