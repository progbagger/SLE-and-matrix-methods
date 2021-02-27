#pragma once
#include "SLE.h"

Vector Lagrange_interpolation_builder(const Vector&, const Vector&); // создание полинома Лагранжа
void printPolynomial(const Vector&); // вывод полинома на экран
double calculatePolynomial(const Vector&, const double); // вычислить значение полинома в точке
