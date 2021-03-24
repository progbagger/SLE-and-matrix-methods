#include "Polynomial.h"
#include <fstream>

Polynomial::Polynomial() { deg = 0; }
Polynomial::Polynomial(const Polynomial& that) { deg = that.deg; p = that.p; }
Polynomial::Polynomial(const size_t& d, const Vector& pol) { deg = d; p = pol; }
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
	for (size_t i = 0; i <= result.getSize(); i++)
		result[i] = (i + 1) * this->get(i + 1);
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
	for (size_t i = 0; i <= this->getSize(); i++)
		this[i] *= a;
	this->shrink();
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
	for (size_t i = 0; i <= this->getSize(); i++)
		this[i] /= a;
	this->shrink();
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
		Polynomial w(0, { 1 });
		double wxk = 1;
		for (size_t j = 0; j <= result.getSize(); j++)
			if (i != j)
			{
				w *= Polynomial(1, { -1 * x.get(j), 1 });
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
		Polynomial w(0, { 1 });
		for (size_t j = 0; j < i; j++)
			w *= Polynomial(1, { -1 * x.get(j), 1 });
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
