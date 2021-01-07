#pragma once
#include "Vector.h"
#include <fstream>
#include <initializer_list>
#include <iomanip>

using namespace std;

class Matrix // класс для работы с квадратными матрицами
{
    double** M;
    size_t size; // размер матрицы

    void rawClean();
    void rawCopy(const Matrix&);

public:

    Matrix(const Matrix&);
    Matrix();
    ~Matrix();
    Matrix(const size_t&);
    Matrix(const initializer_list<Vector>&);
    Matrix(const initializer_list<double>&);
    Matrix& operator = (const Matrix&);
    Matrix& operator = (const initializer_list<Vector>&);
    Matrix& operator = (const initializer_list<double>&);

    size_t getSize() const; // размер матрицы
    double*& operator [] (const size_t&); // с возможностью изменения
    double*& get(const size_t&) const; // константный

    // стандартные операторы для работы с матрицами
    Matrix operator + (const Matrix&) const; // сложение матриц
    Matrix& operator += (const Matrix&);
    Matrix operator - (const Matrix&) const; // вычитание матриц
    Matrix& operator -= (const Matrix&);
    Matrix operator * (const Matrix&) const; // перемножение матриц
    Matrix& operator *= (const Matrix&);
    Matrix operator * (const double&) const; // умножение матрицы на константу, также описана коммутативность
    Matrix& operator *= (const double&);
    Matrix operator ! () const; // транспонирование

    // сравнение матриц
    bool operator == (const Matrix&) const;
    bool operator != (const Matrix&) const;

    void swap(const size_t&, const size_t&); // смена строк местами
    double det() const; // вычисление определителя матрицы

    Matrix reflect() const; // нахождение обратной матрицы методом Жордана-Гаусса
    Matrix H() const; // матрица отражений
    Matrix diag() const; // диагональная матрица
    Matrix lowerTriangle() const; // нижняя треугольная матрица
    Matrix upperTriangle() const; // верхняя треугольная матрица

    // возвращает вектор с диагонали
    Vector diagV() const;

    // спектральная норма матрицы
    double sNorm() const;

    // нахождение собственных чисел и собственных векторов
    Matrix Hausholder() const; // преобразование Хаусхолдера
    Matrix Givens() const; // преобразование Гивенса
    void QR(const double&) const; // QR-алгоритм

    void RQI(const double&, const double&) const; // Обратные итерации со сдвигом с соотношением Рэлея

    // операторы ввода/вывода матрицы
    friend istream& operator >> (istream&, Matrix&);
    friend ostream& operator << (ostream&, const Matrix&);

    // коммутативные операторы
    friend Matrix operator * (const double&, const Matrix&); // умножение матрицы на константу
    friend Vector operator * (const Matrix&, const Vector&); // умножение матрицы на вектор
};
