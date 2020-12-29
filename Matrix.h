#pragma once
#include "Vector.h"
#include <fstream>
#include <initializer_list>

using namespace std;

class Matrix // класс для работы с квадратными матрицами
{
    double** M;
    size_t size;

    void rawClean();
    void rawCopy(const Matrix&);

public:

    Matrix(const Matrix&);
    Matrix();
    ~Matrix();
    Matrix(const size_t&);
    Matrix(const initializer_list<Vector>&);
    Matrix& operator = (const Matrix&);
    Matrix& operator = (const initializer_list<Vector>&);

    size_t getSize() const;
    double*& operator [] (const int&); // с возможностью изменения
    double*& get(const int&) const;

    // стандартные операторы для работы с матрицами
    Matrix operator + (const Matrix&);
    Matrix operator - (const Matrix&);
    Matrix operator * (const Matrix&);
    Matrix operator * (const double&);
    Matrix operator !();
    bool operator == (const Matrix&) const;
    bool operator != (const Matrix&) const;

    void swap(const int&, const int&); // смена строк местами
    double det() const; // вычисление определителя матрицы

    Matrix reflect(); // нахождение обратной матрицы методом Жордана-Гаусса
    Matrix H(); // матрица отражений
    Matrix diag(); // диагональная матрица
    Matrix lowerTriangle(); // нижняя треугольная матрица
    Matrix upperTriangle(); // верхняя треугольная матрица

    // нахождение собственных чисел и собственных векторов
    void QR();

    // операторы ввода/вывода матрицы
    friend istream& operator >> (istream&, Matrix&);
    friend ostream& operator << (ostream&, const Matrix&);

    // коммутативные операторы
    friend Matrix operator * (const double&, const Matrix&); // умножение матрицы на константу
    friend Vector operator * (const Matrix&, const Vector&); // умножение матрицы на вектор
};
