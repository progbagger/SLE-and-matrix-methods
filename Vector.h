#pragma once
#include <fstream>
#include <initializer_list>

using namespace std;

class Vector // ����� ��� ������ � ���������
{
    double* V;
    size_t size; // ������ �������

    void rawClean();
    void rawCopy(const Vector&);

public:

    Vector(const Vector&);
    Vector();
    ~Vector();
    Vector& operator = (const Vector&);
    Vector(const size_t&);
    Vector(const size_t&, const double&);
    Vector(const size_t&, const Vector&); // ������ � ������-������
    Vector(const initializer_list<double>&);
    Vector& operator = (const initializer_list<double>&);

    size_t getSize() const; // ������ �������
    double& operator [] (const size_t&); // � ������������ ���������
    double& get(const size_t&) const; // �����������

    double operator * (const Vector&) const; // ��������� ������������
    Vector operator / (const double&) const; // ������� �� ���������
    Vector& operator /= (const double&);
    Vector& operator *= (const double&);
    Vector operator - (const Vector&) const; // ��������� ��������
    Vector& operator -= (const Vector&);
    Vector operator + (const Vector&) const; // �������� ��������
    Vector& operator += (const Vector&);
    bool operator ! () const; // ��������� ������
    bool operator == (const Vector&) const;
    bool operator != (const Vector&) const;

    double infNorm() const; // ����� �� �������������
    double euclidNorm() const; // ��������� �����
    double abs() const; // ������ �������
    void swap(const size_t& i1, const size_t& i2); // �������� �������� �������

    friend Vector operator * (const double&, const Vector&); // ������������� ��������� �� ���������
    friend Vector operator * (const Vector&, const double&); // ������������� ��������� �� ���������

    // ��������� �����/������ �������
    friend istream& operator >> (std::istream&, Vector&);
    friend ostream& operator << (std::ostream&, const Vector&);
};