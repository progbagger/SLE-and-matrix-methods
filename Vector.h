#pragma once
#include <fstream>
#include <initializer_list>

using namespace std;

class Vector // класс для работы с векторами
{
    double* V;
    size_t size; // размер вектора

    void rawClean();
    void rawCopy(const Vector&);

public:

    Vector(const Vector&);
    Vector();
    ~Vector();
    Vector& operator = (const Vector&);
    Vector(const size_t&);
    Vector(const size_t&, const double&);
    Vector(const size_t&, const Vector&); // размер и вектор-массив
    Vector(const initializer_list<double>&);
    Vector& operator = (const initializer_list<double>&);

    size_t getSize() const; // размер вектора
    double& operator [] (const size_t&); // с возможностью изменения
    double& get(const size_t&) const; // константный

    double operator * (const Vector&) const; // скалярное произведение
    Vector operator / (const double&) const; // деление на константу
    Vector& operator /= (const double&);
    Vector& operator *= (const double&);
    Vector operator - (const Vector&) const; // вычитание векторов
    Vector& operator -= (const Vector&);
    Vector operator + (const Vector&) const; // сложение векторов
    Vector& operator += (const Vector&);
    bool operator ! () const; // ненулевой вектор
    bool operator == (const Vector&) const;
    bool operator != (const Vector&) const;

    double infNorm() const; // норма на бесконечности
    double euclidNorm() const; // евклидова норма
    double abs() const; // модуль вектора
    void swap(const size_t& i1, const size_t& i2); // поменять элементы местами

    friend Vector operator * (const double&, const Vector&); // коммутативное умножение на константу
    friend Vector operator * (const Vector&, const double&); // коммунативное умножение на константу

    // операторы ввода/вывода вектора
    friend istream& operator >> (std::istream&, Vector&);
    friend ostream& operator << (std::ostream&, const Vector&);
};