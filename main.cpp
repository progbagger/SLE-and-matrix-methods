#include <fstream>
#include <iomanip>
#include "SLE.h"
#include "Polynomial.h"
#define max(a, b) (a > b) ? a : b
using namespace std;

ifstream fin("input.txt");
ofstream fout("output.txt");

/*
* --------------------------------------------------------Список методов--------------------------------------------------------
* Решение СЛАУ и поиск собственных значений и векторов:
* --- Прямые методы решения СЛАУ (применяются к объектам типа SLE)
* - метод Гаусса: Gauss()
* - метод отражений: HR()
* --- Итерационные методы решения СЛАУ (применяются к объектам типа SLE)
* - метод Зейделя: HZ(const double& - точность, const Vector& - начальное приближение)
* - метод Якоби: Jacobi(const double& - точность, const Vector& - начальное приближение)
* - метод сопряжённых градиентов: SGrd(const double& - точность, const Vector& - начальное приближение)
* - трёхчленная формула реализации метода Ричардсона с чебышёвскими коэффициентами: Rchd3(const double& - точность,
* const Vector& - начальное приближение, const double& - левая граница спектра, const double& - правая граница спектра)
* --- Методы поиска собственных значений (применя.тся к объектам типа Matrix)
* - QR-метод: QR(const double& - точность)
* - метод обратных итераций со сдвигом с соотношением Рэлея: RQI(const double& - точность, const double& - приближение с. ч.)
* Численный анализ:
* --- Интерполирование:
* - Интерполяционный многочлен Лагранжа: int_L(const Vector& - иксы узлов, const Vector& - игреки узлов)
* - Интерполяционный многочлен Ньютона: int_N(const Vector& - иксы узлов, const Vector& - игреки узлов)
* - Построение узлов Чебышёва: Vector Cheb(const size_t& - необходимое количество узлов)
*/

int main()
{
    fout << fixed << setprecision(8); // установка точности на 8 знаков после запятой

    /*//////////////////////////////////////////////////////////
    * //////////////////////////////////////////////////////////
    * //////////////////////////////////////////////////////////
    * Раскомментировать нужное или написать своё, если требуется
    * //////////////////////////////////////////////////////////
    * //////////////////////////////////////////////////////////
    *///////////////////////////////////////////////////////////

    /*
    * Ввод начальных данных
    */

    //int N; fin >> N; // размерность матрицы или системы

    ////////////////////////////////////////////////////////////

    //double e; fin >> e; // точность вычислений (не для прямых методов)
    //Vector v0(N); fin >> v0; // начальное приближение (для итерационных методов)
    //double alpha, beta; fin >> alpha >> beta; // границы собственных чисел (для метода Ричардсона)
    //double lambda_e; fin >> lambda_e; // начальное приближение с. ч. (для метода обратных итераций)
    //SLE system(N); fin >> system; // для решения систем
    //Matrix A(N); fin >> A; // для поиска собственных значений

    /*
    * Вызов методов и печать их результата
    */

    // Прямые методы
    //fout << system.Gauss(); // метод Гаусса
    //fout << system.HR(); // метод отражений

    // Итерационные методы
    //system.HZ(e, v0); // метод Зейделя
    //system.Jacobi(e, v0); // метод Якоби
    //system.SGrd(e, v0); // метод сопряжённых градиентов
    //system.Rchd3(e, v0, alpha, beta); // метод Ричардсона

    // Методы поиска собственных значений и векторов
    //A.QR(e); // QR-метод
    //A.RQI(e, lambda_e); // метод обратных итераций

    ////////////////////////////////////////////////////////////

    // Интерполяционный многочлен Лагранжа

    /*
    double a, b;
    fin >> a >> b;
    double x0;
    fin >> x0;
    size_t n;
    fin >> n;
    Vector x(n), y(n);
    for (size_t i = 0; i < n; i++)
        fin >> x[i] >> y[i];
    Polynomial Lagr = int_L(x, y);
    fout << "f(" << x0 << ") = " << Lagr(x0) << endl << "f(x) = " << Lagr;
    */

    // Интерполяционный многочлен Ньютона

    /*
    double a, b;
    fin >> a >> b;
    double x0;
    fin >> x0;
    size_t n;
    fin >> n;
    Vector x(n), y(n);
    for (size_t i = 0; i < n; i++)
        fin >> x[i] >> y[i];
    Polynomial Newt = int_N(x, y);
    fout << "f(" << x0 << ") = " << Newt(x0) << endl << "f(x) = " << Newt;
    */

    // Узлы Чебышёва

    /*
    double a = -1, b = 1;
    size_t n = 4;
    Vector x{ -1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0 };
    Vector y(n);
    for (size_t i = 0; i < x.getSize(); i++)
        y[i] = exp(x.get(i));
    Vector xCh = Cheb(4);
    Vector yCh(n);
    for (size_t i = 0; i < xCh.getSize(); i++)
        yCh[i] = exp(xCh.get(i));
    // полином Лагранжа по равноотстоящим узлам x
    Polynomial EL = int_L(x, y);
    fout << "Equidistant nodes:\nx = " << x << "y = " << y << "L(x) = " << EL;
    // полином Лагранжа по узлам Чебышёва
    Polynomial CL = int_L(xCh, yCh);
    fout << "Chebyshev nodes:\nx = " << xCh << "y = " << yCh << "L(x) = " << CL;
    // поиск экстремума ошибки
    // Так как в программе нельзя задать явно многочлен, одним из слагаемых которого
    // является экспонента, придётся отдельно считать производную как сумму
    Polynomial dfEL = EL.df(), dfCL = CL.df();
    double ac = 0.6, bc = 0.8, eps = 0.000000001;
    // точность 10^-9 удобна для нас, поскольку требуемая точность вывода 8 знаков после запятой
    // интервал найден по графику
    double result = 0, result_prev = 0;
    // для равноотстоящих узлов
    do
    {
        result_prev = result;
        result = (ac + bc) / 2;
        if ((dfEL(result) - exp(result)) * (dfEL(ac) - exp(ac)) <= 0)
            bc = result;
        else if ((dfEL(result) - exp(result)) * (dfEL(bc) - exp(bc)) <= 0)
            ac = result;
    } while (abs(result - result_prev) >= eps);
    fout << "x = " << result << endl;
    fout << "Maximum for equidistant nodes: " << abs(exp(result) - EL(result)) << endl;
    // для узлов Чебышёва
    ac = 0.6, bc = 0.8;
    result = result_prev = 0;
    do
    {
        result_prev = result;
        result = (ac + bc) / 2;
        if ((dfCL(result) - exp(result)) * (dfCL(ac) - exp(ac)) <= 0)
            bc = result;
        else if ((dfCL(result) - exp(result)) * (dfCL(bc) - exp(bc)) <= 0)
            ac = result;
    } while (abs(result - result_prev) >= eps);
    fout << "x = " << result << endl;
    fout << "Maximum for Chebyshev nodes: " << abs(exp(result) - CL(result)) << endl;
    */
    
    fin.close();
    fout.close();
    return 0;
}
