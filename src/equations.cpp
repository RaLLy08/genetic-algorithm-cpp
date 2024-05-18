#include "equations.h"

double rosensbrock(double x, double y) {
    return 100 * pow(y - x * x, 2) + pow(1 - x, 2);
}

double sphere(double a, double b) {
    return a * a + b * b;
}

double deJong(double a, double b) {
    return 3905.93 - 100 * (a * a - b) * (a * a - b) - (1 - a) * (1 - a);
}

double venter(double a, double b) {
    return 100 * pow(b - a * a, 2) + pow(1 - a, 2);
}

double ackley(double a, double b) {
    return -20 * exp(-0.4 * sqrt(0.5 * (a * a + b * b))) - exp(0.5 * (cos(2 * M_PI * a) + cos(2 * M_PI * b))) + 20 + M_E;
}

double schwefel(double a, double b) {
    return 418.9829 * 2 - (a * sin(sqrt(abs(a))) + b * sin(sqrt(abs(b))));
}

double griewank(double a, double b) {
    return 1 + (a * a + b * b) / 4000 - cos(a) * cos(b / sqrt(2));
}

double styblinski(double a, double b) {
    return 0.5 * (a * a * a * a - 16 * a * a + 5 * a + b * b * b * b - 16 * b * b + 5 * b);
}

double rastrigin(double a, double b) {
    return 20 + a * a - 10 * cos(2 * M_PI * a) + b * b - 10 * cos(2 * M_PI * b);
}