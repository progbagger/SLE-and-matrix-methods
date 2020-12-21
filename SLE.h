#pragma once
#include "Matrix.h"
#include "Vector.h"
#include <fstream>

using namespace std;

class SLE // класс для работы с системами линейных алгебраических уравнений
{
private:

    int size;
    Matrix M;
    Vector b;

    void rawCopy(const SLE& that) { this->size = that.size; M = that.M; b = that.b; }

public:

    SLE();
    SLE(const SLE& that);
    SLE(Matrix m, Vector v);
    SLE& operator = (const SLE& that);
    SLE(int s);
    ~SLE();

    Vector c_getb() const;
    Vector& getb();
    Matrix c_getM() const;
    Matrix& getM();
    int getSize() const;

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