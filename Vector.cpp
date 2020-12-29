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
    for (int i = 0; i <= this->size; i++)
        this->V[i] = that.V[i];
}

size_t Vector::getSize() const { return size; }
double& Vector::operator [] (const int& r) { return V[r]; }
double& Vector::get(const int& r) const { return V[r]; }
Vector::Vector(const Vector& that) { rawCopy(that); }
Vector::Vector() { size = 0; V = new double; }
Vector::~Vector() { rawClean(); }

double* Vector::next(double* el)
{
    el++;
    return el;
}

Vector::Iterator::Iterator() : el(nullptr), belong(nullptr) {};
Vector::Iterator::Iterator(double& c_el, Vector& c_belong) : el(&c_el), belong(&c_belong) {};
Vector::Iterator::Iterator(const Iterator& that)
{
    el = that.el;
    belong = that.belong;
}
Vector::Iterator::~Iterator() {};

Vector::Iterator& Vector::Iterator::operator = (const Iterator& that)
{
    if (this != &that)
    {
        el = that.el;
        belong = that.belong;
    }
    return *this;
}

Vector::Iterator Vector::Iterator::operator ++ (int)
{
    Iterator tmp(*this);
    el = belong->next(el);
    return tmp;
}

Vector::Iterator Vector::Iterator::operator ++ ()
{
    el = belong->next(el);
    return *this;
}

bool Vector::Iterator::operator == (const Iterator& that)
{
    return (*el == *that.el && *belong == *that.belong);
}

bool Vector::Iterator::operator != (const Iterator& that)
{
    return (*el != *that.el || *belong != *that.belong);
}

double& Vector::Iterator::operator * ()
{
    return *el;
}

Vector::Iterator Vector::begin()
{
    return Iterator(V[0], *this);
}

Vector::Iterator Vector::end()
{
    return Iterator(V[size], *this);
}

Vector& Vector::operator = (const Vector& that) {
    if (this != &that)
    {
        rawClean();
        rawCopy(that);
    }
    return *this;
}

Vector::Vector(const size_t& s) {
    size = s;
    V = new double[size + 1];
    for (int i = 0; i < size; i++)
        V[i] = 1;
    V[size] = 0;
}

Vector::Vector(const size_t& s, double* m) {
    size = s;
    V = new double[size + 1];
    for (int i = 0; i < size; i++)
        V[i] = m[i];
    V[size] = 0;
}

Vector::Vector(const initializer_list<double>& list) : size(list.size())
{
    V = new double[size + 1];
    size_t i = 0;
    for (const auto& el : list)
    {
        V[i] = el;
        i++;
    }
    V[size] = 0;
}

Vector& Vector::operator = (const initializer_list<double>& list)
{
    *this = Vector(list);
    return *this;
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

bool Vector::operator != (const Vector& that)
{
    bool result = false;
    for (int i = 0; i < size; i++)
        if (V[i] != that.V[i])
        {
            result = true;
            break;
        }
    return result;
}

bool Vector::operator == (const Vector& that)
{
    bool result = true;
    for (int i = 0; i < size; i++)
        if (V[i] != that.V[i])
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
