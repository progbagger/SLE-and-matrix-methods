#pragma once
#include "Matrix.h"
#include "Vector.h"
#include <fstream>
#include <initializer_list>
#include <iomanip>

using namespace std;

class SLE // ����� ��� ������ � ��������� �������� �������������� ���������
{
private:

    size_t size; // ������ �������
    Matrix M; // ������� �������
    Vector b; // ������ ��������� ������

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

    // ��������� ���������
    bool operator == (const SLE&) const;
    bool operator != (const SLE&) const;

    Vector c_getb() const; // �����������
    Vector& getb(); // ��������� �������
    Matrix c_getM() const; // �����������
    Matrix& getM(); // ��������� �������
    size_t getSize() const; // ������ �������

    // ������ ������
    Vector Gauss() const; // ����� ������
    Vector HR() const; // ����� ���������

    // ������������ ������
    void iView(); // ���������� ������� � ����, ���������� ��� ��������
    void HZ(const double&, const Vector&) const; // ����� ������-�������
    void Jacobi(const double&, const Vector&) const; // ����� �����
    void SGrd(const double&, const Vector&) const; // ����� ���������� ����������
    void Rchd3(const double&, const Vector&, const double&, const double&) const; // ���������� ������� ���������� ������ ���������� � ������������ �����������

    // ��������� �����/������ ����
    friend istream& operator >> (istream&, SLE&);
    friend ostream& operator << (ostream&, const SLE&);
};