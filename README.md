# Прямые методы решения СЛАУ - в `SLE.h` и `SLE.cpp`
* Метод Гаусса: `Vector SLE::Gauss()`
* Метод отражений: `Vector SLE::HR()`
# Итерационные методы решения СЛАУ - в `SLE.h` и `SLE.cpp`
* Метод Зейделя: `void SLE::HZ(const double&, const Vector&)`
* Метод Якоби: `void SLE::Jacobi(const double&, const Vector&)`
* Метод сопряжённых градиентов: `void SLE::SGrd(const double&, const Vector&)`
* Трёхчленная формула реализации метода Ричардсона с чебышёвскими параметрами: `void SLE::Rchd3(const double&, const Vector&, const double&, const double&)`
# Методы поиска собственных значений - в `Matrix.h` и `Matrix.cpp`
* QR-метод: `void Matrix::QR(const double&)`
* Метод обратных итераций со сдвигом с соотношением Рэлея: `void Matrix::RQI(const double&, const double&)`
# Методы численного анализа - в `Polynomial.h` и `Polynomial.cpp`
## Интерполирование
* Вывод и ввод полинома осуществляются через операторы `<<` и `>>` соответственно
* Значение полинома в точке вычисляется с помощью оператора `()`
* Построение полинома Лагранжа в данных узлах и их значениях: `Polynomial Polynomial::int_L(const Vector& - узлы, const Vector& - значения в узлах)`
* Построение полинома Ньютона в данных узлах и их значениях: `Polynomial Polynomial::int_N(const Vector& - узлы, const Vector& - значения в узлах)`
* Построение узлов Чебышёва: `Vector Cheb(const size_t& - количество необходимых узлов)`
## Численное интегрирование
* `Vector Vector steady_grid(const size_t&, const double&, const double&)` - создание равномерной сетки
* `double mid_rect(const size_t& - размерность равномерной сетки, const double& - левая граница отрезка, const double& - правая граница отрезка, double (*)(const double&) - данная функция)` - формула центральных прямоугольников
* `double trapecia(const size_t&, const double&, const double&, double (*)(const double&))` - формула трапеций
* `double Simpson(const size_t&, const double&, const double&, double (*)(const double&))` - формула Симпсона
* `void Gaussian(const double& - точность вычислений n, const size_t&, const double&, const double&, double (*)(const double&))` - квадратурная формула Гаусса на n узлах
  * `Polynomial Legendre(const size_t&)` - построение полинома Лежандра
## Решение нелинейных уравнений
* `pair<size_t, double> dichotomy(const double&, const double&, const double&, double (*)(const double&))` - метод дихотомии (используйте `dichotomy(...).first` для использования количества итераций и `dichotomy(...).second` для использования результата метода
* pair<size_t, double> newt(const double&, const double&, const double&, double (*)(const double&), double (*)(const double&)) - метод Ньютона (используйте `newt(...).first` для использования количества итераций и `newt(...).second` для использования результата метода
  * Также имеются перегруженные версии данной функции, которые принимают начальное приближение и вместо указателей на функции используют в качестве аргументов полиномы
