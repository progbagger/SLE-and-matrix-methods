#include "Matrix.h"
#include "SLE.h"

using namespace std;

void Matrix::rawClean()
{
    for (size_t i = 0; i < this->size; i++)
        delete[] this->M[i];
    delete[] this->M;
}

void Matrix::rawCopy(const Matrix& that)
{
    this->size = that.size;
    this->M = new double* [this->size];
    for (size_t i = 0; i < this->size; i++)
        this->M[i] = new double[this->size];
    for (size_t i = 0; i < this->size; i++)
        for (size_t j = 0; j < this->size; j++)
            this->M[i][j] = that.M[i][j];
}


Matrix::Matrix(const Matrix& that) { rawCopy(that); }
Matrix::Matrix() { size = 0; this->M = nullptr; }
Matrix::~Matrix() { rawClean(); }

Matrix& Matrix::operator = (const Matrix& that) 
{
    if (this != &that) {
        rawClean();
        rawCopy(that);
    }
    return *this;
}

size_t Matrix::getSize() const { return size; }
double*& Matrix::operator [] (const size_t& r) { return M[r]; } // с изменением
double*& Matrix::get(const size_t& r) const { return M[r]; }


Matrix::Matrix(const size_t& s) : size(s) 
{
    M = new double* [size];
    for (size_t i = 0; i < size; i++)
        M[i] = new double[size];
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++)
            if (i == j)
                M[i][j] = 1;
            else
                M[i][j] = 0;
}

Matrix::Matrix(const initializer_list<Vector>& list) : size(list.begin()->getSize())
{
    M = new double* [size];
    for (size_t i = 0; i < size; i++)
        M[i] = new double[size];
    size_t i = 0;
    for (const auto& v : list)
    {
        for (size_t j = 0; j < size; j++)
            M[i][j] = v.get(j);
        i++;
    }
    for (i; i < size; i++)
        for (size_t j = 0; j < size; j++)
            M[i][j] = 0;
}

Matrix::Matrix(const initializer_list<double>& list)
{
    size = ceil(sqrt(list.size()));
    M = new double* [size];
    for (size_t i = 0; i < size; i++)
        M[i] = new double[size];
    size_t ii = 0;
    for (const auto& el : list)
    {
        M[ii / size][ii % size] = el;
        ii++;
    }
    for (size_t i = ii / size; i < size; i++)
        for (size_t j = ii % size; j < size; j++)
            M[i][j] = 0;
}

Matrix& Matrix::operator = (const initializer_list<Vector>& list)
{
    *this = Matrix(list);
    return *this;
}

Matrix Matrix::H() const // нахождение матрицы отражений
{
    Matrix result(size);
    Matrix R = *this;
    for (size_t i = 0; i < size - 1; i++) {
        Vector t(size - i);
        for (size_t j = 0; j < size - i; j++)
            t[j] = R.get(j + i)[i];
        Vector b(this->size - i, t); // вектор b - первый столбец нашей матрицы размера n - i, где i - итерация главного цикла for
        for (size_t j = 0; j < size - i; j++)
            if (j == 0)
                t[j] = 1;
            else
                t[j] = 0;
        Vector c(size - i, t); // вектор c - единично-координатный вектор
        Vector temp = b - (b.abs() * c);
        Vector w = (temp / sqrt(2 * (b * temp))); // ключевой вектор, из которого строится матрица отражений для каждой итерации
        Matrix Omega(size - i);
        for (size_t j = 0; j < size - i; j++)
            for (size_t k = 0; k < size - i; k++)
                Omega[j][k] = w.get(j) * w.get(k);
        Matrix E(size - i);
        Matrix EU = E - (2 * Omega); // матрица отражений для каждой итерации
        Matrix U(size); // обобщённая матрица отражений (размера n)
        for (size_t j = 0; j < size; j++)
            for (size_t k = 0; k < size; k++)
                if (j >= i && k >= i) {
                    U[j][k] = EU.get(j - i)[k - i];
                }
                else {
                    if (j == k)
                        U[j][k] = 1;
                    else
                        U[j][k] = 0;
                }
        result = U * result; // формирование матрицы отражений в каждой итерации цикла
        R = U * R;
    }
    return result;
}

// стандартные операторы для работы с матрицами
Matrix Matrix::operator + (const Matrix& that) const
{
    Matrix result(size);
    for (size_t i = 0; i < result.size; i++)
        for (size_t j = 0; j < result.size; j++)
            result.M[i][j] = M[i][j] + that.M[i][j];
    return result;
}

Matrix& Matrix::operator += (const Matrix& that)
{
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++)
            M[i][j] += that.M[i][j];
    return *this;
}

Matrix Matrix::operator - (const Matrix& that) const
{
    Matrix result(size);
    for (size_t i = 0; i < result.size; i++)
        for (size_t j = 0; j < result.size; j++)
            result.M[i][j] = M[i][j] - that.M[i][j];
    return result;
}

Matrix& Matrix::operator -= (const Matrix& that)
{
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++)
            M[i][j] -= that.M[i][j];
    return *this;
}

Matrix Matrix::operator * (const Matrix& that) const
{
    Matrix result(size);
    for (size_t i = 0; i < result.size; i++) {
        result.M[i][i] = 0;
        for (size_t j = 0; j < result.size; j++)
            for (size_t k = 0; k < result.size; k++)
                result.M[i][j] += M[i][k] * that.M[k][j];
    }
    return result;
}

Matrix& Matrix::operator *= (const Matrix& that)
{
    Matrix result(size);
    for (size_t i = 0; i < result.size; i++) {
        result.M[i][i] = 0;
        for (size_t j = 0; j < result.size; j++)
            for (size_t k = 0; k < result.size; k++)
                result.M[i][j] += M[i][k] * that.M[k][j];
    }
    *this = result;
    return *this;
}

Matrix Matrix::operator * (const double& coeff) const
{
    Matrix result(size);
    for (size_t i = 0; i < result.size; i++)
        for (size_t j = 0; j < result.size; j++)
            result.M[i][j] = M[i][j] * coeff;
    return result;
}

Matrix& Matrix::operator *= (const double& coeff)
{
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++)
            M[i][j] *= coeff;
    return *this;
}

Matrix operator * (const double& coeff, const Matrix& that)
{
    Matrix result(that.getSize());
    for (size_t i = 0; i < result.getSize(); i++)
        for (size_t j = 0; j < result.getSize(); j++)
            result[i][j] = that.get(i)[j] * coeff;
    return result;
}

Matrix Matrix::operator !() const
{
    Matrix result(size);
    for (size_t i = 0; i < result.size; i++)
        for (size_t j = i; j < result.size; j++)
        {
            if (i == j)
                result.M[i][j] = M[i][j];
            else
            {
                double temp = M[i][j];
                result.M[i][j] = M[j][i];
                result.M[j][i] = temp;
            }
        }
    return result;
}

bool Matrix::operator == (const Matrix& that) const
{
    bool result = true;
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++)
            if (M[i][j] != that.M[i][j])
            {
                result = false;
                break;
            }
    return result;
}

bool Matrix::operator != (const Matrix& that) const
{
    bool result = false;
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++)
            if (M[i][j] != that.M[i][j])
            {
                result = true;
                break;
            }
    return result;
}

void Matrix::swap(const size_t& i1, const size_t& i2) { // смена строк местами
    for (size_t i = 0; i < size; i++) {
        double temp = M[i1][i];
        M[i1][i] = M[i2][i];
        M[i2][i] = temp;
    }
}

// вычисление определителя матрицы
double Matrix::det() const {
    Matrix T = *this;
    double result = 1;
    for (size_t i = 0; i < T.size - 1; i++) {
        if (!T.M[i][i])
            for (size_t j = i + 1; j < T.size; j++) {
                if (T.M[j][i]) {
                    T.swap(j, i);
                    break;
                }
            }
        for (size_t j = i + 1; j < T.size; j++) {
            double temp = T.M[j][i] / T.M[i][i];
            for (size_t k = 0; k < T.size; k++)
                T.M[j][k] -= T.M[i][k] * temp;
        }
    }
    for (size_t i = 0; i < T.size; i++)
        result *= T.M[i][i];
    return result;
}

// нахождение обратной матрицы методом Жордана-Гаусса
Matrix Matrix::reflect() const
{
    Matrix E(this->size);
    Matrix T = *this;
    // прямой ход метода Гаусса
    for (size_t i = 0; i < T.size - 1; i++) {
        if (!T.M[i][i])
            for (size_t j = i + 1; j < T.size; j++) {
                if (T.M[j][i]) {
                    T.swap(j, i);
                    E.swap(j, i);
                    break;
                }
            }
        for (size_t j = i + 1; j < T.size; j++) {
            double temp = T.M[j][i] / T.M[i][i];
            for (size_t k = 0; k < T.size; k++) {
                E.M[j][k] -= E.M[i][k] * temp;
                T.M[j][k] -= T.M[i][k] * temp;
            }
        }
    }
    // обратный ход метода Гаусса
    for (int i = T.size - 1; i > 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            double temp = T.M[j][i] / T.M[i][i];
            for (int k = T.size - 1; k >= 0; k--) {
                E.M[j][k] -= E.M[i][k] * temp;
                T.M[j][k] -= T.M[i][k] * temp;
            }
        }
    }
    for (size_t i = 0; i < T.size; i++) {
        for (size_t j = 0; j < this->size; j++)
            E.M[i][j] /= T.M[i][i];
        T.M[i][i] = 1;
    }
    return E;
}

Matrix Matrix::diag() const
{
    Matrix result(size);
    for (size_t i = 0; i < size; i++)
        result[i][i] = M[i][i];
    return result;
}

Matrix Matrix::lowerTriangle() const
{
    Matrix result(size);
    for (size_t i = 0; i < size; i++)
        for (size_t j = i; j < size; j++)
            if (i == j)
                result[j][i] = 0;
            else
                result[j][i] = M[j][i];
    return result;
}

Matrix Matrix::upperTriangle() const
{
    Matrix result(size);
    for (size_t i = 0; i < size; i++)
        for (size_t j = i; j < size; j++)
            if (i == j)
                result[i][j] = 0;
            else
                result[i][j] = M[i][j];
    return result;
}

Vector Matrix::diagV() const
{
    Vector result(size);
    for (size_t i = 0; i < size; i++)
        result[i] = M[i][i];
    return result;
}

double Matrix::sNorm() const
{
    double result = 0;
    for (size_t i = 0; i < size; i++)
        result += M[i][i] * M[i][i];
    return sqrt(result);
}

// операторы ввода/вывода матрицы
istream& operator >> (istream& in, Matrix& that)
{
    for (size_t i = 0; i < that.getSize(); i++)
        for (size_t j = 0; j < that.getSize(); j++)
            in >> that[i][j];
    return in;
}

ostream& operator << (ostream& out, const Matrix& that)
{
    for (size_t i = 0; i < that.getSize(); i++)
    {
        for (size_t j = 0; j < that.getSize(); j++)
            out << that.get(i)[j] << '\t';
        out << endl;
    }
    out << endl;
    return out;
}

Vector operator * (const Matrix& that, const Vector& v)
{
    Vector result(v.getSize());
    for (size_t i = 0; i < result.getSize(); i++)
    {
        result[i] = 0;
        for (size_t j = 0; j < result.getSize(); j++)
            result[i] += that.get(i)[j] * v.get(j);
    }
    return result;
}

// преобразование Хаусхолдера для первого шага QR-метода
Matrix Matrix::Hausholder() const
{
    Matrix result(size);
    double sum = 0; // для коэффициента s
    for (size_t i = 1; i < size; i++)
        sum += M[i][0] * M[i][0];
    sum = sqrt(sum);
    double s = -1 * (M[1][0] / abs(M[1][0])) * sum;
    double mu = 1 / sqrt(2 * s * (s - M[1][0]));
    Vector w(size);
    w[0] = 0;
    w[1] = M[1][0] - s;
    for (size_t i = 2; i < size; i++)
        w[i] = M[i][0];
    w *= mu;
    Matrix W(size);
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++)
            W[i][j] = w[i] * w[j];
    Matrix H = Matrix(size) - 2 * W;
    result = H * *this * H;
    return result;
}

// Построение матрицы отражений для каждого шага QR-алгоритма
Matrix Matrix::QR_reflection() const
{
    ofstream fout("output.txt", ios::app);
    fout << fixed << setprecision(8);
    Matrix Q(size), R = *this;
    // выполняется n - 1 раз
    for (size_t i = 0; i < size - 1; i++)
    {
        double s = 0;
        for (size_t j = i; j < size; j++)
            s += R[j][i] * R[j][i];
        s = -1 * (R[i][i] / abs(R[i][i])) * sqrt(s);
        double mu = 1 / sqrt(2 * s * (s - R[i][i]));
        Vector w(size, 0);
        w[i] = R[i][i] - s;
        for (size_t j = i + 1; j < size; j++)
            w[j] = R[j][i];
        w *= mu;
        Matrix W(size);
        for (size_t j = 0; j < size; j++)
            for (size_t k = 0; k < size; k++)
                W[j][k] = w[j] * w[k];
        Matrix H = Matrix(size) - 2 * W;
        Q *= H;
        R = H * R;
    }
    return R * Q;
}

// QR метод нахождения собственных чисел и векторов матрицы
void Matrix::QR(const double& e) const
{
    ofstream fout("output.txt", ios::app);
    fout << fixed << setprecision(8);
    /*
    * Первый шаг - матрица Хессинберга, получаемая в результате преобразования Хессенберга
    */
    Matrix B = (*this).Hausholder();
    Vector result = B.diagV() /* запись в результат текущих приближений собственных значений */, resn1(size);
    /*
    * Следующие шаги - матрицы отражений
    */
    size_t iterationCounter = 0;
    do
    {
        resn1 = result;
        B = B.QR_reflection();
        result = B.diagV();
        iterationCounter++;
    } while ((result - resn1).infNorm() > e); // условие конца процесса
    // печтаь результата
    fout << "Iteration: " << iterationCounter << endl;
    fout << "Eigenvalues: " << endl;
    for (size_t i = 0; i < size; i++)
        fout << i + 1 << ") " << result[i] << endl;
    fout << endl;
}

// Обратные итерации со сдвигом с соотношением Рэлея
void Matrix::RQI(const double& e, const double& lambda_e) const
{
    ofstream fout("output.txt", ios::app);
    fout << fixed << setprecision(8);
    /*
    * Шаг 0: подбор вектора, евклидома норма которого равна нулю
    */
    Vector xkn1(size, 0);
    xkn1[0] = 1;
    /*
    * Шаг 1: для k = 1, 2, ... :
    * 1.1: вычисление ro = (A * xk, xk), в случае первого шага это xkn1;
    * 1.2: вычисление yk из уравнения (A - ro * E) * yk = xkn1;
    * 1.3: нормирование yk: xk = yk / ||yk||;
    * 1.4: проверка ro на сходимость и проверка конца процесса.
    * (из Вержбицкого)
    * На деле алгоритм пришлось слегка изменить...
    */
    Vector yk(size), xk(size);
    double ro = lambda_e, difference = e + 1;
    size_t iterationCount = 0;
    while (difference > e) // условие конца - норма соседних приближений собственных векторов
    {
        Matrix M = (*this - ro * Matrix(size));
        // если матрица вырождена, то найдено точное собственное значение
        if (!M.det())
            break;
        SLE sle(M, xkn1);
        yk = sle.HR();
        xk = yk / yk.euclidNorm();
        difference = (xk - xkn1).infNorm();
        xkn1 = xk;
        ro = (*this * xk) * xk;
        iterationCount++;
    }
    // печать результата
    fout << "Iteration: " << iterationCount;
    if (!iterationCount)
        fout << " (The eigenvalue was guessed by the approximation)";
    fout << "\nEigenvalue: " << ro << endl;
    if (iterationCount)
    {
        fout << "Eigenvector:\n";
        for (size_t i = 0; i < size; i++)
            fout << "[" << i + 1 << "] = " << xk[i] << endl;
    }
    fout << endl;
}