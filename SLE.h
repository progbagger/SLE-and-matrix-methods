#pragma once
#include "Matrix.h"
#include "Vector.h"
#include <fstream>
#include <initializer_list>
#include <iomanip>

using namespace std;

class SLE // класс для работы с системами линейных алгебраических уравнений
{
private:

    size_t size; // размер системы
    Matrix M; // матрица системы
    Vector b; // вектор свободных членов

    void rawCopy(const SLE& that);

public:

    SLE();
    SLE(const SLE& that);
    SLE(Matrix m, Vector v);
    SLE& operator = (const SLE& that);
    SLE(size_t s);
    SLE(const pair<Matrix, Vector>&);
    SLE& operator = (const pair<Matrix, Vector>&);
    ~SLE();

    // операторы сравнения
    bool operator == (const SLE&) const;
    bool operator != (const SLE&) const;

    Vector c_getb() const; // константный
    Vector& getb(); // получение вектора
    Matrix c_getM() const; // константный
    Matrix& getM(); // получение матрицы
    size_t getSize() const; // размер системы

    // прямые методы
    Vector Gauss() const; // метод Гаусса
    Vector HR() const; // метод отражений

    // итерационные методы
    void iView(); // приведение системы к виду, пригодному для итерации
    void HZ(const double&, const Vector&) const; // метод Гаусса-Зейделя
    void Jacobi(const double&, const Vector&) const; // метод Якоби
    void SGrd(const double&, const Vector&) const; // метод сопряжённых градиентов
    void Rchd3(const double&, const Vector&, const double&, const double&) const; // трёхчленная формула реализации метода Ричардсона с чебышёвскими параметрами

    // операторы ввода/вывода СЛАУ
    friend istream& operator >> (istream&, SLE&);
    friend ostream& operator << (ostream&, const SLE&);
};