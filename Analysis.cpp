#include "Analysis.h"
#include <fstream>

void printPolynomial(const Vector& polynomial)
{
	bool check = false; // определить, напечатано ли первое ненулевое слагаемое
	std::ofstream fout("output.txt", ios::app);
	fout << "f(x) = ";
	for (int i = polynomial.getSize() - 1; i >= 0; i--)
	{
		if (polynomial.get(i))
		{
			if (!check)
			{
				fout << polynomial.get(i);
				check = true;
			}
			else
			{
				if (polynomial.get(i) < 0)
					fout << " - ";
				else
					fout << " + ";
				fout << polynomial.get(i);
			}
			if (i)
				fout << "x^" << i;
		}
	}
	fout << endl << endl;
}

double calculatePolynomial(const Vector& polynomial, const double x0)
{
	double result = polynomial.get(0);
	for (size_t i = 1; i < polynomial.getSize(); i++)
		result += polynomial.get(i) * pow(x0, i);
	return result;
}

Vector Lagrange_interpolation_builder(const Vector& nodes_x, const Vector& nodes_y) // создание полинома
{
	/*
	* Создание интерполяционного полинома Лагранжа с помощью решения системы
	* линейных уравнений Ax = b, где в качестве x выступает вектор коэффициентов
	* полинома, а "иксы" полинома известны (подставляем их сами)
	*/
	size_t size = nodes_x.getSize();
	Vector result(size);
	Matrix A(size); Vector b(size);
	for (size_t i = 0; i < size; i++)
	{
		for (size_t j = 0; j < size; j++)
			A[i][j] = pow(nodes_x.get(i), j); // заполняем матрицу коэфф-тами x (подставляем узлы)
		b[i] = nodes_y.get(i); // заполняем вектор b известными значениями в узлах
	}
	SLE system(A, b);
	result = system.HR();
	return result;
}
