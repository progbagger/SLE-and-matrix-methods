#pragma once
#include "Matrix.h"
#include "Vector.h"
#include <fstream>
#include <initializer_list>

using namespace std;

class SLE // класс для работы с системами линейных алгебраических уравнений
{
private:

    size_t size;
    Matrix M;
    Vector b;

    void rawCopy(const SLE& that);

public:

    class Iterator;

    Iterator begin();
    Iterator end();

    SLE();
    SLE(const SLE& that);
    SLE(Matrix m, Vector v);
    SLE& operator = (const SLE& that);
    SLE(size_t s);
    SLE(const pair<Matrix, Vector>&);
    SLE& operator = (const pair<Matrix, Vector>&);
    ~SLE();

    // операторы сравнения
    bool operator == (const SLE&);
    bool operator != (const SLE&);

    Vector c_getb() const;
    Vector& getb();
    Matrix c_getM() const;
    Matrix& getM();
    size_t getSize() const;

    // прямые методы
    Vector Gauss();
    Vector HR();

    // итерационные методы
    void iView();
    void HZ(const double&, const Vector&);
    void Jacobi(const double&, const Vector&);
    void SGrd(const double&, const Vector&);
    void Rchd3(const double&, const Vector&, const double&, const double&);

    // операторы ввода/вывода СЛАУ
    friend istream& operator >> (istream&, SLE&);
    friend ostream& operator << (ostream&, const SLE&);
};
