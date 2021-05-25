#include <fstream>
#include "Polynomial.h"
#define EPS 0.000000001

Polynomial::Polynomial() { deg = 0; }
Polynomial::Polynomial(const Polynomial& that) { deg = that.deg; p = that.p; }
Polynomial::Polynomial(const Vector& pol) { deg = pol.getSize() - 1; p = pol; }
Polynomial::Polynomial(const size_t& d) { deg = d; p = Vector(deg + 1, 0); }
Polynomial::~Polynomial() {}

Polynomial& Polynomial::operator = (const Polynomial& that) { deg = that.deg; p = that.p; return *this; }

size_t Polynomial::getSize() const { return deg; }
double Polynomial::get(const size_t& d) const { return p.get(d); }
double& Polynomial::operator [] (const size_t& r) { return p[r]; }

/*
* Функция, вычисляющая степень полинома, получившегося
* в результате арифметических операций
*/
size_t Polynomial::calculateDeg() const
{
	size_t result = 0;
	for (size_t i = deg; i > 0; i--)
		if (this->get(i))
		{
			result = i;
			break;
		}
	return result;
}

/*
* Функция, удаляющая нулевые высшие степени из полинома
*/
void Polynomial::shrink()
{
	size_t endDeg = this->calculateDeg();
	if (deg != endDeg)
	{
		Polynomial result(endDeg);
		for (size_t i = 0; i <= endDeg; i++)
			result[i] = p.get(i);
		*this = result;
	}
}

Polynomial Polynomial::operator + (const Polynomial& that) const
{
	size_t endDeg = 0;
	if (deg <= that.deg)
		endDeg = that.deg;
	else
		endDeg = deg;
	Polynomial result(endDeg);
	for (size_t i = 0; i <= this->getSize(); i++)
		result[i] += this->get(i);
	for (size_t i = 0; i <= that.getSize(); i++)
		result[i] += that.get(i);
	result.shrink();
	return result;
}

Polynomial& Polynomial::operator += (const Polynomial& that)
{
	size_t endDeg = 0;
	if (deg <= that.deg)
		endDeg = that.deg;
	else
		endDeg = deg;
	Polynomial result(endDeg);
	for (size_t i = 0; i <= this->getSize(); i++)
		result[i] += this->get(i);
	for (size_t i = 0; i <= that.getSize(); i++)
		result[i] += that.get(i);
	result.shrink();
	*this = result;
	return *this;
}

Polynomial Polynomial::operator - (const Polynomial& that) const
{
	size_t endDeg = 0;
	if (deg <= that.deg)
		endDeg = that.deg;
	else
		endDeg = deg;
	Polynomial result(endDeg);
	for (size_t i = 0; i <= this->getSize(); i++)
		result[i] += this->get(i);
	for (size_t i = 0; i <= that.getSize(); i++)
		result[i] -= that.get(i);
	result.shrink();
	return result;
}

Polynomial& Polynomial::operator -= (const Polynomial& that)
{
	size_t endDeg = 0;
	if (deg <= that.deg)
		endDeg = that.deg;
	else
		endDeg = deg;
	Polynomial result(endDeg);
	for (size_t i = 0; i <= this->getSize(); i++)
		result[i] += this->get(i);
	for (size_t i = 0; i <= that.getSize(); i++)
		result[i] -= that.get(i);
	result.shrink();
	*this = result;
	return *this;
}

Polynomial Polynomial::operator * (const Polynomial& that) const
{
	Polynomial result(that.getSize() + this->getSize());
	for (size_t i = 0; i <= this->getSize(); i++)
		for (size_t j = 0; j <= that.getSize(); j++)
			result[i + j] += this->get(i) * that.get(j);
	result.shrink();
	return result;
}

Polynomial& Polynomial::operator *= (const Polynomial& that)
{
	Polynomial result(that.getSize() + this->getSize());
	for (size_t i = 0; i <= this->getSize(); i++)
		for (size_t j = 0; j <= that.getSize(); j++)
			result[i + j] += this->get(i) * that.get(j);
	result.shrink();
	*this = result;
	return *this;
}

Polynomial Polynomial::df() const
{
	Polynomial result(this->getSize() - 1);
	for (size_t i = 1; i <= this->getSize(); i++)
		result[i - 1] = i * this->get(i);
	result.shrink();
	return result;
}
Polynomial Polynomial::Df() const
{
	Polynomial result(this->getSize() + 1);
	for (size_t i = 1; i <= this->getSize() + 1; i++)
		result[i] = this->get(i - 1) / i;
	result.shrink();
	return result;
}

double Polynomial::operator () (const double& x) const
{
	double result = this->get(0);
	for (size_t i = 1; i <= this->getSize(); i++)
		result += this->get(i) * pow(x, i);
	return result;
}

Polynomial operator * (const double& a, const Polynomial& pol)
{
	Polynomial result = pol;
	for (size_t i = 0; i <= result.getSize(); i++)
		result[i] *= a;
	result.shrink();
	return result;
}

Polynomial& Polynomial::operator *= (const double& a)
{
	Polynomial result(this->getSize());
	for (size_t i = 0; i <= result.getSize(); i++)
		result[i] = this->get(i) * a;
	result.shrink();
	*this = result;
	return *this;
}

Polynomial Polynomial::operator / (const double& a)
{
	Polynomial result = *this;
	for (size_t i = 0; i <= this->getSize(); i++)
		result[i] /= a;
	result.shrink();
	return result;
}

Polynomial& Polynomial::operator /= (const double& a)
{
	Polynomial result(this->getSize());
	for (size_t i = 0; i <= result.getSize(); i++)
		result[i] = this->get(i) / a;
	result.shrink();
	*this = result;
	return *this;
}

istream& operator >> (istream& in, Polynomial& pol)
{
	for (size_t i = 0; i <= pol.getSize(); i++)
		in >> pol[i];
	return in;
}

ostream& operator << (ostream& out, const Polynomial& pol)
{
	bool check = false; // определить, напечатано ли первое ненулевое слагаемое
	for (int i = pol.getSize(); i >= 0; i--)
	{
		if (pol.get(i))
		{
			if (!check)
			{
				out << pol.get(i);
				check = true;
			}
			else
			{
				if (pol.get(i) < 0)
					out << " - ";
				else
					out << " + ";
				out << abs(pol.get(i));
			}
			if (i)
				out << "x^" << i;
		}
	}
	out << endl;
	return out;
}

Polynomial int_L(const Vector& x, const Vector& y)
{
	Polynomial result(x.getSize() - 1);
	for (size_t i = 0; i <= result.getSize(); i++)
	{
		// построение полинома w(x)
		Polynomial w(Vector(1, 1));
		double wxk = 1;
		for (size_t j = 0; j <= result.getSize(); j++)
			if (i != j)
			{
				w *= Polynomial({ -1 * x.get(j), 1 });
				wxk *= x.get(i) - x.get(j);
			}
		result += y.get(i) * (w / wxk);
	}
	return result;
}

Polynomial int_N(const Vector& x, const Vector& y)
{
	Polynomial result(x.getSize() - 1);
	/*
	* Составляем вектор разделённых разностей,
	* необходимых для вычисления полинома
	*/
	Vector calculations(y.getSize() + 1, 0);
	calculations[0] = y.get(0);
	for (size_t i = 1; i < calculations.getSize(); i++)
		for (size_t k = 0; k <= i; k++)
		{
			double denominator = 1;
			for (size_t j = 0; j <= i; j++)
				if (k != j)
					denominator *= x.get(k) - x.get(j);
			calculations[i] += y.get(k) / denominator;
		}
	/*
	* Составление самого полинома
	*/
	for (size_t i = 0; i <= x.getSize() - 1; i++)
	{
		Polynomial w(Vector(1, 1));
		for (size_t j = 0; j < i; j++)
			w *= Polynomial({ -1 * x.get(j), 1 });
		result += calculations.get(i) * w;
	}
	return result;
}

Vector Cheb(const size_t& n) // построение узлов Чебышёва для задания 1.3
{
	Vector result(n);
	for (size_t i = 1; i <= result.getSize(); i++)
		result[i - 1] = cos((((2.0 * i) - 1.0) / (2.0 * n)) * acos(-1.0));
	return result;
}

double fu1(const double& x)
{
	double result = 0;
	result = 1 / (x * x + 1);
	return result;
}

double fu2(const double& x)
{
	double result = 0;
	result = (1 + x * x) * (exp(x) - exp(-x)) - 1;
	return result;
}

double dfu2(const double& x)
{
	double result = 0;
	result = exp(x) + exp(-x) + 2 * x * exp(x) + x * x * exp(x) - 2 * x * exp(-x) + x * x * exp(-x);
	return result;
}

// функция создания равномерной сетки на n отрезков
Vector steady_grid(const size_t& n, const double& a, const double& b)
{
	Vector result(n + 1);
	result[0] = a; result[n] = b;
	for (size_t i = 1; i < n; i++)
		result[i] = result[i - 1] + (b - a) / n;
	return result;
}

// формула центральных прямоугольников
double mid_rect(const size_t& n, const double& a, const double& b, double (*function)(const double&))
{
	double result = 0;
	// создадим сетку размера n + 1, чтобы на ней было n отрезков разбиения
	Vector x = steady_grid(n, a, b);
	// формула центральных прямоугольников собственной персоной
	for (size_t i = 0; i < n; i++)
		result += function((x[i] + x[i + 1]) / 2) * (x[i + 1] - x[i]);
	return result;
}

// формула трапеций
double trapecia(const size_t& n, const double& a, const double& b, double (*function)(const double&))
{
	double result = 0;
	// создадим сетку размера n + 1, чтобы на ней было n отрезков разбиения
	Vector x = steady_grid(n, a, b);
	// формула трапеций...
	for (size_t i = 0; i < n; i++)
		result += ((function(x[i]) + function(x[i + 1])) / 2) * (x[i + 1] - x[i]);
	return result;
}

// формула Симпсона
double Simpson(const size_t& n, const double& a, const double& b, double (*function)(const double&))
{
	double result = 0;
	// создадим сетку размера n + 1, чтобы на ней было n отрезков разбиения
	Vector x = steady_grid(n, a, b);
	// формула Симпсона...
	for (size_t i = 0; i < n; i++)
		result += (((function(x[i]) + 4 * function((x[i] + x[i + 1]) / 2) + function(x[i + 1])) / 6) * (x[i + 1] - x[i]));
	return result;
}

Polynomial Legendre(const size_t& n)
{
	Polynomial result(Vector(1, 1));
	for (size_t i = 0; i < n; i++)
		result *= Polynomial({ -1, 0, 1 });
	for (size_t i = 0; i < n; i++)
		result = result.df();
	// вычисление факториала
	long int factorial = 1;
	for (long int i = 1; i <= n; i++)
		factorial *= i;
	// получение итогового многочлена домножением на 1 / (2^n * n!)
	result *= 1 / (pow(2, n) * factorial);
	return result;
}

// метод дихотомии решения неоинейного уравнения f(x) = 0
pair<size_t, double> dichotomy(const double& eps, const double& a, const double& b, double (*function)(const double&))
{
	double result = a, presult = b, ac = a, bc = b;
	size_t counter = 0; // счетчик итераций
	do // цикл должен сработать хотя бы один раз
	{
		++counter;
		presult = result; // для выяснения достижения указанной точности
		result = (ac + bc) / 2;
		if (!function(result))
			break;
		else if (function(ac) * function(result) < 0)
			bc = result;
		else if (function(result) * function(bc) < 0)
			ac = result;
	} while (abs(result - presult) >= eps);
	pair<size_t, double> ans{counter, result}; // возвращаем пару - количество итераций и найденный корень
	return ans;
}

// метод Ньютона решения нелинейного уравнения f(x) = 0
pair<size_t, double> newt(const double& eps, const double& a, const double& b, double (*function)(const double&), double (*dfunction)(const double&))
{
	double result = (a + b) / 2, presult = 0;
	size_t counter = 0;
	if (!function(result))
		goto ret;
	do
	{
		++counter;
		presult = result;
		result = presult - (function(presult) / dfunction(presult));
	} while (abs(result - presult) >= eps);
ret:
	pair<size_t, double> ans{ counter, result };
	return ans;
}

pair<size_t, double> newt(const double& eps, double (*function)(const double&), double (*dfunction)(const double&), const double& x0)
{
	double result = x0, presult = 0;
	size_t counter = 0;
	if (!function(result))
		goto ret;
	do
	{
		++counter;
		presult = result;
		result = presult - (function(presult) / dfunction(presult));
	} while (abs(result - presult) >= eps);
ret:
	pair<size_t, double> ans{ counter, result };
	return ans;
}

pair<size_t, double> newt(const double& eps, Polynomial& function, Polynomial& dfunction, const double& x0)
{
	double result = x0, presult = 0;
	size_t counter = 0;
	if (!function(result))
		goto ret;
	do
	{
		++counter;
		presult = result;
		result = presult - (function(presult) / dfunction(presult));
	} while (abs(result - presult) >= eps);
ret:
	pair<size_t, double> ans{ counter, result };
	return ans;
}

// квадратурная формула Гаусса
void Gaussian(const double& eps, const size_t& n, const double& a, const double& b, double (*function)(const double&))
{
	fstream fout("output.txt", ios::app);
	fout << fixed << setprecision(8);
	double result = 0;
	Polynomial result_pol;
	// найдём узлы - корни многочлена Лежандра n-й степени
	Vector x(n, 0), c(n, 0);
	Polynomial L = Legendre(n), dL = L.df();
	for (size_t i = 0; i < n; i++)
	{
		// начальное приближение для i-го корня
		x[i] = cos((acos(-1) * (4 * (i + 1) - 1)) / (4 * n + 2));
		// поиск самого корня
		x[i] = newt(eps, L, dL, x[i]).second;
	}
	// перерасчёт узлов с отрезка [-1, 1] на нужный отрезок [a, b]
	if (a != -1 || b != 1)
		for (size_t i = 0; i < n; i++)
			x[i] = ((a + b) / 2) + (x[i] * ((b - a) / 2));
	// подсчёт значений функции в узлах
	Vector y(n, 0);
	for (size_t i = 0; i < n; i++)
		y[i] = function(x[i]);
	// подсчёт весов квадратурной формулы
	for (size_t i = 0; i < n; i++) {
		Polynomial weight(Vector(1, 1));
		double denominator = 1;
		for (size_t j = 0; j < n; j++)
			if (i != j) {
				weight *= Polynomial(Vector{ -1 * x[j], 1 });
				denominator *= x[i] - x[j];
			}
		weight /= denominator;
		c[i] = (weight.Df())(b) - (weight.Df())(a);
	}
	// сама квадратурная формула
	for (size_t i = 0; i < n; i++)
		result += c[i] * function(x[i]);
	//result *= ((b - a) / 2);
	fout << "Узлы:\n" << x << "Значения функции в узлах:\n" << y << "Веса квадратурной формулы:\n" << c << "Результат = " << result << endl << endl;
	fout.close();
}
