#pragma once
#include "Vector.h"
#include <fstream>

using namespace std;

class Matrix // ����� ��� ������ � ����������� ���������
{
    double** M;
    int size;

    void rawClean()
    {
        for (int i = 0; i < this->size; i++)
            delete[] this->M[i];
        delete[] this->M;
    }

    void rawCopy(const Matrix& that)
    {
        this->size = that.size;
        this->M = new double* [this->size];
        for (int i = 0; i < this->size; i++)
            this->M[i] = new double[this->size];
        for (int i = 0; i < this->size; i++)
            for (int j = 0; j < this->size; j++)
                this->M[i][j] = that.M[i][j];
    }

public:

    Matrix(const Matrix&);
    Matrix();
    ~Matrix();
    Matrix(const int&);
    Matrix& operator = (const Matrix&);

    int getSize() const;
    double*& operator [] (const int&); // � ������������ ���������
    double*& get(const int&) const;

    // ����������� ��������� ��� ������ � ���������
    Matrix operator + (const Matrix&);
    Matrix operator - (const Matrix&);
    Matrix operator * (const Matrix&);
    Matrix operator * (const double&);
    Matrix operator !();

    void swap(const int&, const int&); // ����� ����� �������
    double det() const; // ���������� ������������ �������

    Matrix reflect(); // ���������� �������� ������� ������� �������-������
    Matrix H(); // ������� ���������
    Matrix diag(); // ������������ �������
    Matrix lowerTriangle(); // ������ ����������� �������
    Matrix upperTriangle(); // ������� ����������� �������

    // ���������� ����������� ����� � ����������� ��������
    void QR();

    // ��������� �����/������ �������
    friend istream& operator >> (istream&, Matrix&);
    friend ostream& operator << (ostream&, const Matrix&);

    // ������������� ���������
    friend Matrix operator * (const double&, const Matrix&); // ��������� ������� �� ���������
    friend Vector operator * (const Matrix&, const Vector&); // ��������� ������� �� ������
};