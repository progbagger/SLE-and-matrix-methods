#pragma once
#include "SLE.h"

class Polynomial
{
	size_t deg; // степень многочлена
	Vector p; // сам многочлен

public:

	Polynomial();
	Polynomial(const Polynomial&);
	~Polynomial();

	Polynomial& operator = (const Polynomial&);

	Polynomial(const Vector&);
	Polynomial(const size_t&);

	size_t getSize() const;
	double get(const size_t&) const;
	double& operator [] (const size_t&);

	size_t calculateDeg() const;
	void shrink();

	// операторы для работы с полиномами
	Polynomial operator + (const Polynomial&) const;
	Polynomial& operator += (const Polynomial&);
	Polynomial operator - (const Polynomial&) const;
	Polynomial& operator -= (const Polynomial&);
	Polynomial operator * (const Polynomial&) const;
	Polynomial& operator *= (const Polynomial&);

	Polynomial df() const; // производный многочлен

	double operator () (const double&) const; // значение полинома в точке

	// операторы ввода/вывода полинома
	friend istream& operator >> (istream&, Polynomial&);
	friend ostream& operator << (ostream&, const Polynomial&);

	// операторы умножения и деления многочлена на константу
	friend Polynomial operator * (const double&, const Polynomial&);
	Polynomial& operator *= (const double&);
	Polynomial operator / (const double&);
	Polynomial& operator /= (const double&);
};

// функции создания интерполяционного многочлена
Polynomial int_L(const Vector&, const Vector&);
Polynomial int_N(const Vector&, const Vector&);

Vector Cheb(const size_t&);

// функции численного интегрирования

Vector steady_grid(const size_t&, const double&, const double&);

double mid_rect(const size_t&, const double&, const double&, double (*)(const double&));
double trapecia(const size_t&, const double&, const double&, double (*)(const double&));
double Simpson(const size_t&, const double&, const double&, double (*)(const double&));

double fu1(const double&); // для задания 2

double fu2(const double&); // для задания 3
double dfu2(const double&);

void Gaussian(const double&, const size_t&, const double&, const double&, double (*)(const double&));

Polynomial Legendre(const size_t&); // построение полинома Лежандра степени n

pair<size_t, double> dichotomy(const double&, const double&, const double&, double (*)(const double&));
pair<size_t, double> newt(const double&, const double&, const double&, double (*)(const double&), double (*)(const double&));
pair<size_t, double> newt(const double&, const double&, const double&, double (*)(const double&), double (*)(const double&), const double&); // версия с начальным приближением
pair<size_t, double> newt(const double&, const double&, const double&, Polynomial&, Polynomial&, const double&);
