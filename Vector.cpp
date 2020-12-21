#include "Vector.h"
#include <fstream>

using namespace std;

int Vector::getSize() const { return size; }
double& Vector::operator [] (const int& r) { return V[r]; }
double& Vector::get(const int& r) const { return V[r]; }
Vector::Vector(const Vector& that) { rawCopy(that); }
Vector::Vector() { size = 0; V = nullptr; }
Vector::~Vector() { rawClean(); }

Vector& Vector::operator = (const Vector& that) {
    if (this != &that)
    {
        rawClean();
        rawCopy(that);
    }
    return *this;
}

Vector::Vector(const int& s) {
    size = s;
    V = new double[size];
    for (int i = 0; i < size; i++)
        V[i] = 1;
}

Vector::Vector(const int& s, double* m) {
    size = s;
    V = new double[size];
    for (int i = 0; i < size; i++)
        V[i] = m[i];
}

// ниже приведены стандартные операторы для работы с векторами
double Vector::operator *(const Vector& that) {
    double result = 0;
    for (int i = 0; i < size; i++)
        result += V[i] * that.V[i];
    return result;
}

// коммутативность операторов
Vector operator *(const double& a, const Vector& that) {
    Vector result(that.getSize());
    for (int i = 0; i < result.getSize(); i++)
        result[i] = a * that.get(i);
    return result;
}

Vector operator * (const Vector& that, const double& a) {
    Vector result(that.getSize());
    for (int i = 0; i < result.getSize(); i++)
        result[i] = that.get(i) * a;
    return result;
}

Vector Vector::operator /(const double& a) {
    Vector result(size);
    for (int i = 0; i < result.size; i++)
        result[i] = V[i] / a;
    return result;
}

Vector Vector::operator -(const Vector& that) {
    Vector result(size);
    for (int i = 0; i < result.size; i++)
        result[i] = V[i] - that.get(i);
    return result;
}

Vector Vector::operator +(const Vector& that) {
    Vector result(size);
    for (int i = 0; i < result.size; i++)
        result.V[i] = V[i] + that.V[i];
    return result;
}

bool Vector::operator ! () const {
    bool result = true;
    for (int i = 0; i < size; i++)
        if (V[i] != 0)
        {
            result = false;
            break;
        }
    return result;
}

// норма вектора
double Vector::infNorm() const{
    double result = std::abs(V[0]);
    for (int i = 1; i < size; i++)
        if (std::abs(V[i]) > result)
            result = std::abs(V[i]);
    return result;
}

// операторы ввода и вывода вектора
istream& operator >> (istream& in, Vector& that) {
    for (int i = 0; i < that.getSize(); i++)
        in >> that[i];
    return in;
}

ostream& operator << (ostream& out, const Vector& that) {
    out << "(";
    for (int i = 0; i < that.getSize() - 1; i++)
        out << that.get(i) << ", ";
    out << that.get(that.getSize() - 1) << ')' << endl << endl;
    return out;
}

double Vector::abs() const{
    double result = 0;
    for (int i = 0; i < size; i++)
        result += V[i] * V[i];
    result = sqrt(result);
    return result;
}

void Vector::swap(const int& i1, const int& i2) {
    double temp = V[i1];
    V[i1] = V[i2];
    V[i2] = temp;
}