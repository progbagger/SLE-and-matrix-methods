#include "SLE.h"
#include <fstream>
#include <iomanip>

using namespace std;

SLE::SLE() { size = 0; M = Matrix(); b = Vector(); }
SLE::SLE(const SLE& that) { this->rawCopy(that); }
SLE::SLE(Matrix m, Vector v) { size = m.getSize(); M = m; b = v; }
SLE::SLE(int s) { size = s; M = Matrix(s); b = Vector(s); }
SLE::~SLE() {}

SLE& SLE::operator = (const SLE& that) {
    if (this != &that)
        this->rawCopy(that);
    return *this;
}

// геттеры
Vector SLE::c_getb() const { return b; }
Vector& SLE::getb() { return b; }
Matrix SLE::c_getM() const { return M; }
Matrix& SLE::getM() { return M; }
int SLE::getSize() const { return size; }

// операторы ввода/вывода СЛАУ
istream& operator >> (istream& in, SLE& that)
{
    for (int i = 0; i < that.getSize(); i++)
    {
        for (int j = 0; j < that.getSize(); j++)
            in >> that.getM()[i][j];
        in >> that.getb()[i];
    }
    return in;
}

ostream& operator << (ostream& out, const SLE& that)
{
    for (int i = 0; i < that.getSize(); i++)
    {
        for (int j = 0; j < that.getSize(); j++)
            out << that.c_getM().get(i)[j] << '\t';
        out << that.c_getb().get(i) << endl;
    }
    out << endl;
    return out;
}

void SLE::iView()
{
    for (int i = 0; i < size; i++)
    {
        double t = M[i][i];
        M[i][i] = 1;
        if (i)
        {
            double temp = M[i][i];
            M[i][i] = M[i][0];
            M[i][0] = temp;
        }
        for (int j = 1; j < this->size; j++)
           M[i][j] /= -t;
        b[i] /= t;
    }
}

/* ПРЯМЫЕ МЕТОДЫ */

// решение системы методом Гаусса
Vector SLE::Gauss()
{
    Matrix Mt = M; Vector bt = b; // чтобы не испортить исходные
    Vector result(size);
    // прямой ход
    for (int i = 0; i < size - 1; i++)
    {
        if (!Mt[i][i])
            for (int j = i + 1; j < size; j++)
            {
                if (Mt[i][i])
                {
                    Mt.swap(j, i);
                    bt.swap(j, i);
                    break;
                }
            }
        for (int j = i + 1; j < size; j++)
        {
            double temp = Mt[j][i] / Mt[i][i];
            for (int k = 0; k < size; k++)
            {
                Mt[j][k] -= Mt[i][k] * temp;
            }
            bt[j] -= bt[i] * temp;
        }
    }
    // нахождение иксов
    for (int i = size - 1; i >= 0; i--)
    {
        result[i] = bt[i];
        for (int j = i + 1; j < size; j++)
            result[i] -= result[j] * Mt[i][j];
        result[i] /= Mt[i][i];
    }
    return result;
}

Vector SLE::HR() // решение СЛАУ методом отражений
{
    Matrix MM = M; Vector bb = b; // делаем доп. систему чтобы не испортить исходную
    Vector result(size);
    for (int i = 0; i < result.getSize(); i++)
        result[i] = 0;
    Matrix H = MM.H();
    // преобразовываем систему, умножая матрицу системы и вектор свободных членов на матрицу отражений слева
    MM = H * MM; bb = H * bb;
    // поиск "иксов"
    for (int i = this->size - 1; i >= 0; i--)
    {
        result[i] = bb[i];
        for (int j = i + 1; j < this->size; j++)
            result[i] -= result[j] * MM[i][j];
        result[i] /= MM[i][i];
    }
    return result;
}

/* ИТЕРАЦИОННЫЕ МЕТОДЫ */

// Метод Гаусса-Зейделя
void SLE::HZ(const double& e, const Vector& ee)
{
    //e - точность вычислений
    Vector result(size);
    SLE t = *this;
    Matrix D = t.M.diag(), L = -1 * t.M.lowerTriangle();
    Matrix H = (D - L).reflect(); // вычисление (D - L)^-1
    Vector x0 = ee; // начальное приближение, равное вектору свободных членов в преобразованной системе
    result = x0 - H * ((getM() * x0) - getb()); // первая итерация
    int m = 1;
    while ((result - x0).infNorm() > e) // условие конца - норма "соседних" вычисленных решений на бесконечности должна быть меньше или равна точности
    {
        x0 = result;
        result = x0 - (H * ((getM() * x0) - getb()));
        m++;
    }
    ofstream fout("output.txt", ios::app);
    fout << "m = " << m << endl << "x = " << result;
    fout.close();
}

// Метод Якоби
void SLE::Jacobi(const double& e, const Vector& ee) {
    ofstream fout("output.txt", ios::app);
    fout << fixed << setprecision(8);
    Vector result(size);
    Matrix H = (getM().diag()).reflect(); // матрица D^-1
    Vector x0 = ee; // начальный вектор
    result = x0 - H * ((getM() * x0) - getb());
    int m = 1;
    while ((result - x0).infNorm() > e) { // условие конца
        x0 = result;
        result = x0 - H * (getM() * x0 - getb());
        m++;
    }
    fout << "m = " << m << endl << "x = " << result;
    fout.close();
}

// Метод сопряжённых градиентов
void SLE::SGrd(const double& e, const Vector& ee)
{
    ofstream fout("output.txt", ios::app);
    fout << fixed << setprecision(8);
    Vector result = ee, x0 = ee;
    int m = 0;
    SLE t = *this;
    Vector r = (t.getM() * x0) - t.getb(), g = r; // задаём начальные значения вектора невязки и градиента
    if (!g) // проверка на не базисность вектора градиента
    {
        fout << "m = " << m << endl << "x = " << result;
        fout.close();
        return;
    }
    double a = (r * g) / ((t.getM() * g) * g); // коэфф. для вычисления вектора x
    result = x0 - (a * g); r = r - (a * (t.getM() * g)); // первые значения
    ++m;
    if (!r) // проверка...
    {
        fout << "m = " << m << endl << "x = " << result;
        fout.close();
        return;
    }
    double gamma; // коэфф. для вычисления градиента
    while ((result - x0).infNorm() > e) // основной цикл программы
    {
        x0 = result;
        gamma = ((t.getM() * r) * g) / ((t.getM() * g) * g);
        g = r - (gamma * g);
        a = (r * g) / ((t.getM() * g) * g);
        result = x0 - (a * g);
        r = r - (a * (t.getM() * g));
        ++m;
    }
    fout << "m = " << m << endl << "x = " << result;
    fout.close();
}

// метод Ричардсона с чебышёвскими параметрами (трёхчленная формула)
void SLE::Rchd3(const double& e, const Vector& ee, const double& alpha, const double& beta)
{
    ofstream fout("output.txt", ios::app);
    fout << fixed << setprecision(8);
    Vector result(size);
    try {
        if (beta + alpha == 0)
            throw - 1;
    }
    catch (int) {
        return;
    }
    /* Задаём начальные векторы и коэффициенты */
    Vector xkn1 = ee; // вектор x(k - 1)
    double w1 = -1 * (beta - alpha) / (beta + alpha); // начальный коэффициент w
    int m = 0; // номер итерации
    /* первая итерация */
    Vector xk = xkn1 - ((2 / (beta + alpha)) * ((getM() * xkn1) - getb())); // вектор xk, первая итерация
    ++m;
    double wk = w1, wkp1 = w1; // вычисление каждого нового w требует w1 и предыдущего w
    /* Рабочий цикл */
    while ((xk - xkn1).infNorm() > e)
    {
        wkp1 = 1 / ((2 * (1 / w1)) - wk);
        result = xk + (wk * wkp1) * (xk - xkn1) - (2 / (beta + alpha)) * (1 + wk * wkp1) * (getM() * xk - getb());
        xkn1 = xk;
        xk = result;
        wk = wkp1;
        ++m;
    }
    /* печать результата */
    fout << "m = " << m << endl << "x = " << result;
    fout.close();
}
