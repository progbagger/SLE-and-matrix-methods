#pragma once
#include "Vector.h"
#include <fstream>
#include <initializer_list>
#include <iomanip>

using namespace std;

class Matrix // ����� ��� ������ � ����������� ���������
{
    double** M;
    size_t size; // ������ �������

    void rawClean();
    void rawCopy(const Matrix&);

public:

    Matrix(const Matrix&);
    Matrix();
    ~Matrix();
    Matrix(const size_t&);
    Matrix(const initializer_list<Vector>&);
    Matrix(const initializer_list<double>&);
    Matrix& operator = (const Matrix&);
    Matrix& operator = (const initializer_list<Vector>&);
    Matrix& operator = (const initializer_list<double>&);

    size_t getSize() const; // ������ �������
    double*& operator [] (const size_t&); // � ������������ ���������
    double*& get(const size_t&) const; // �����������

    // ����������� ��������� ��� ������ � ���������
    Matrix operator + (const Matrix&) const; // �������� ������
    Matrix& operator += (const Matrix&);
    Matrix operator - (const Matrix&) const; // ��������� ������
    Matrix& operator -= (const Matrix&);
    Matrix operator * (const Matrix&) const; // ������������ ������
    Matrix& operator *= (const Matrix&);
    Matrix operator * (const double&) const; // ��������� ������� �� ���������, ����� ������� ���������������
    Matrix& operator *= (const double&);
    Matrix operator ! () const; // ����������������

    // ��������� ������
    bool operator == (const Matrix&) const;
    bool operator != (const Matrix&) const;

    void swap(const size_t&, const size_t&); // ����� ����� �������
    double det() const; // ���������� ������������ �������

    Matrix reflect() const; // ���������� �������� ������� ������� �������-������
    Matrix H() const; // ������� ���������
    Matrix diag() const; // ������������ �������
    Matrix lowerTriangle() const; // ������ ����������� �������
    Matrix upperTriangle() const; // ������� ����������� �������

    // ���������� ������ � ���������
    Vector diagV() const;

    // ������������ ����� �������
    double sNorm() const;

    // ���������� ����������� ����� � ����������� ��������
    Matrix Hausholder() const; // �������������� �����������
    Matrix QR_reflection() const; // ���������� ������� ��������� ��� ������� ���� QR-���������
    void QR(const double&) const; // QR-��������

    void RQI(const double&, const double&) const; // �������� �������� �� ������� � ������������ �����

    // ��������� �����/������ �������
    friend istream& operator >> (istream&, Matrix&);
    friend ostream& operator << (ostream&, const Matrix&);

    // ������������� ���������
    friend Matrix operator * (const double&, const Matrix&); // ��������� ������� �� ���������
    friend Vector operator * (const Matrix&, const Vector&); // ��������� ������� �� ������
};