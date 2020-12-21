#pragma once
#include "Vector.h"
#include <fstream>

using namespace std;

class Matrix // класс для работы с квадратными матрицами
{
    double** M;
    int size;

    void rawClean()
    {
        for (int i = 0; i < this->size; i++)
            delete[] this->M[i];
        delete[] this->M;
    }

    void rawCopy(const Matrix& that)
    {
        this->size = that.size;
        this->M = new double* [this->size];
        for (int i = 0; i < this->size; i++)
            this->M[i] = new double[this->size];
        for (int i = 0; i < this->size; i++)
            for (int j = 0; j < this->size; j++)
                this->M[i][j] = that.M[i][j];
    }

public:

    Matrix(const Matrix&);
    Matrix();
    ~Matrix();
    Matrix(const int&);
    Matrix& operator = (const Matrix&);

    int getSize() const;
    double*& operator [] (const int&); // с возможностью изменения
    double*& get(const int&) const;

    // стандартные операторы для работы с матрицами
    Matrix operator + (const Matrix&);
    Matrix operator - (const Matrix&);
    Matrix operator * (const Matrix&);
    Matrix operator * (const double&);
    Matrix operator !();

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
