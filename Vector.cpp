#include "Vector.h"
#include <fstream>

using namespace std;

void Vector::rawClean()
{
    delete[] this->V;
}

void Vector::rawCopy(const Vector& that)
{
    this->size = that.size;
    this->V = new double[this->size];
    for (size_t i = 0; i < this->size; i++)
        this->V[i] = that.V[i];
}

size_t Vector::getSize() const { return size; }
double& Vector::operator [] (const size_t& r) { return V[r]; }
double& Vector::get(const size_t& r) const { return V[r]; }
Vector::Vector(const Vector& that) { rawCopy(that); }
Vector::Vector() { size = 0; V = new double; }
Vector::~Vector() { rawClean(); }

Vector& Vector::operator = (const Vector& that)
{
    if (this != &that)
    {
        rawClean();
        rawCopy(that);
    }
    return *this;
}

Vector::Vector(const size_t& s)
{
    size = s;
    V = new double[size];
    for (size_t i = 0; i < size; i++)
        V[i] = 1;
}

Vector::Vector(const size_t& s, const double& el) : size(s)
{
    V = new double[size];
    for (size_t i = 0; i < size; i++)
        V[i] = el;
}

Vector::Vector(const size_t& s, const Vector& m)
{
    size = s;
    V = new double[size];
    for (size_t i = 0; i < size; i++)
        V[i] = m.get(i);
}

Vector::Vector(const initializer_list<double>& list) : size(list.size())
{
    V = new double[size];
    size_t i = 0;
    for (const auto& el : list)
    {
        V[i] = el;
        i++;
    }
}

Vector& Vector::operator = (const initializer_list<double>& list)
{
    *this = Vector(list);
    return *this;
}

// ниже приведены стандартные операторы для работы с векторами
double Vector::operator * (const Vector& that) const
{
    double result = 0;
    for (size_t i = 0; i < size; i++)
        result += V[i] * that.V[i];
    return result;
}

// коммутативность операторов
Vector operator * (const double& a, const Vector& that)
{
    Vector result(that.getSize());
    for (size_t i = 0; i < result.getSize(); i++)
        result[i] = a * that.get(i);
    return result;
}

Vector operator * (const Vector& that, const double& a)
{
    Vector result(that.getSize());
    for (size_t i = 0; i < result.getSize(); i++)
        result[i] = that.get(i) * a;
    return result;
}

Vector Vector::operator / (const double& a) const
{
    Vector result(size);
    for (size_t i = 0; i < result.size; i++)
        result[i] = V[i] / a;
    return result;
}

Vector& Vector::operator /= (const double& a)
{
    for (size_t i = 0; i < size; i++)
        V[i] /= a;
    return *this;
}

Vector& Vector::operator *= (const double& a)
{
    for (size_t i = 0; i < size; i++)
        V[i] *= a;
    return *this;
}

Vector Vector::operator - (const Vector& that) const
{
    Vector result(size);
    for (size_t i = 0; i < result.size; i++)
        result[i] = V[i] - that.get(i);
    return result;
}

Vector& Vector::operator -= (const Vector& that)
{
    for (size_t i = 0; i < size; i++)
        V[i] -= that.get(i);
    return *this;

}

Vector Vector::operator + (const Vector& that) const
{
    Vector result(size);
    for (size_t i = 0; i < result.size; i++)
        result.V[i] = V[i] + that.V[i];
    return result;
}

Vector& Vector::operator += (const Vector& that)
{
    for (size_t i = 0; i < size; i++)
        V[i] += that.V[i];
    return *this;
}

bool Vector::operator ! () const
{
    bool result = true;
    for (size_t i = 0; i < size; i++)
        if (V[i] != 0)
        {
            result = false;
            break;
        }
    return result;
}

bool Vector::operator != (const Vector& that) const
{
    bool result = false;
    for (size_t i = 0; i < size; i++)
        if (V[i] != that.V[i])
        {
            result = true;
            break;
        }
    return result;
}

bool Vector::operator == (const Vector& that) const
{
    bool result = true;
    for (size_t i = 0; i < size; i++)
        if (V[i] != that.V[i])
        {
            result = false;
            break;
        }
    return result;
}

// норма вектора на бесконечности
double Vector::infNorm() const
{
    double result = std::abs(V[0]);
    for (size_t i = 1; i < size; i++)
        if (std::abs(V[i]) > result)
            result = std::abs(V[i]);
    return result;
}

// евклидова норма вектора
double Vector::euclidNorm() const
{
    double result = 0;
    for (size_t i = 0; i < size; i++)
        result += V[i] * V[i];
    return sqrt(result);
}

// операторы ввода и вывода вектора
istream& operator >> (istream& in, Vector& that)
{
    for (size_t i = 0; i < that.getSize(); i++)
        in >> that[i];
    return in;
}

ostream& operator << (ostream& out, const Vector& that)
{
    out << "(";
    for (size_t i = 0; i < that.getSize() - 1; i++)
        out << that.get(i) << ", ";
    out << that.get(that.getSize() - 1) << ')' << endl << endl;
    return out;
}

double Vector::abs() const
{
    double result = 0;
    for (size_t i = 0; i < size; i++)
        result += V[i] * V[i];
    result = sqrt(result);
    return result;
}

void Vector::swap(const size_t& i1, const size_t& i2)
{
    double temp = V[i1];
    V[i1] = V[i2];
    V[i2] = temp;
}
