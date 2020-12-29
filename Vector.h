#pragma once
#include <fstream>
#include <initializer_list>

using namespace std;

class Vector // класс для работы с векторами
{
    double* V;
    size_t size;

    void rawClean();
    void rawCopy(const Vector&);

public:

    double* next(double*);

    class Iterator
    {
        double* el;
        Vector* belong;

    public:

        Iterator();
        Iterator(double&, Vector&);
        Iterator(const Iterator&);
        ~Iterator();
        Iterator& operator = (const Iterator&);

        Iterator operator ++ (int);
        Iterator operator ++ ();

        bool operator == (const Iterator&);
        bool operator != (const Iterator&);

        double& operator * ();
    };

    Iterator begin();
    Iterator end();

    Vector(const Vector&);
    Vector();
    ~Vector();
    Vector& operator = (const Vector&);
    Vector(const size_t&);
    Vector(const size_t&, double*);
    Vector(const initializer_list<double>&);
    Vector& operator = (const initializer_list<double>&);

    size_t getSize() const;
    double& operator [](const int&);
    double& get(const int&) const;

    double operator *(const Vector&);
    Vector operator /(const double&);
    Vector operator -(const Vector&);
    Vector operator +(const Vector&);
    bool operator !() const;
    bool operator == (const Vector&);
    bool operator != (const Vector&);

    double infNorm() const;
    double abs() const;
    void swap(const int& i1, const int& i2);

    friend Vector operator * (const double&, const Vector&);
    friend Vector operator * (const Vector&, const double&);

    // операторы ввода/вывода вектора
    friend istream& operator >> (std::istream&, Vector&);
    friend ostream& operator << (std::ostream&, const Vector&);
};
