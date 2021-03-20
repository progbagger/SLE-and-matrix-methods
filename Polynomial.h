#pragma once
#include "SLE.h"

class Polynomial
{
	size_t deg; // ������� ����������
	Vector p; // ��� ���������

public:

	Polynomial();
	Polynomial(const Polynomial&);
	~Polynomial();

	Polynomial& operator = (const Polynomial&);

	Polynomial(const size_t&, const Vector&);
	Polynomial(const size_t&);

	size_t getSize() const;
	double get(const size_t&) const;
	double& operator [] (const size_t&);

	size_t calculateDeg() const;
	void shrink();

	// ��������� ��� ������ � ����������
	Polynomial operator + (const Polynomial&) const;
	Polynomial& operator += (const Polynomial&);
	Polynomial operator - (const Polynomial&) const;
	Polynomial& operator -= (const Polynomial&);
	Polynomial operator * (const Polynomial&) const;
	Polynomial& operator *= (const Polynomial&);

	Polynomial df() const; // ����������� ���������

	double operator () (const double&) const; // �������� �������� � �����

	// ��������� �����/������ ��������
	friend istream& operator >> (istream&, Polynomial&);
	friend ostream& operator << (ostream&, const Polynomial&);

	// ��������� ��������� � ������� ���������� �� ���������
	friend Polynomial operator * (const double&, const Polynomial&);
	Polynomial& operator *= (const double&);
	Polynomial operator / (const double&);
	Polynomial& operator /= (const double&);
};

// ������� �������� ����������������� ����������
Polynomial int_L(const Vector&, const Vector&);
Polynomial int_N(const Vector&, const Vector&);

double error(const Polynomial&, const Polynomial);