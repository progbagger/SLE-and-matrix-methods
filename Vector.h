#pragma once
#include <fstream>

using namespace std;

class Vector // класс для работы с векторами
{
    double* V;
    int size;

    void rawClean()
    {
        delete[] this->V;
    }

    void rawCopy(const Vector& that)
    {
        this->size = that.size;
        this->V = new double[this->size];
        for (int i = 0; i < this->size; i++)
            this->V[i] = that.V[i];
    }

public:

    Vector(const Vector&);
    Vector();
    ~Vector();
    Vector& operator = (const Vector&);
    Vector(const int&);
    Vector(const int&, double*);

    int getSize() const;
    double& operator [](const int&);
    double& get(const int&) const;

    double operator *(const Vector&);
    Vector operator /(const double&);
    Vector operator -(const Vector&);
    Vector operator +(const Vector&);
    bool operator !() const;

    double infNorm() const;
    double abs() const;
    void swap(const int& i1, const int& i2);

    friend Vector operator * (const double&, const Vector&);
    friend Vector operator * (const Vector&, const double&);

    // операторы ввода/вывода вектора
    friend istream& operator >> (std::istream&, Vector&);
    friend ostream& operator << (std::ostream&, const Vector&);
};
