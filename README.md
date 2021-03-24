# Прямые методы решения СЛАУ - в SLE.h и SLE.cpp
1) Метод Гаусса: Vector SLE::Gauss()
2) Метод отражений: Vector SLE::HR()
# Итерационные методы решения СЛАУ - в SLE.h и SLE.cpp
1) Метод Зейделя: void SLE::HZ(const double&, const Vector&)
2) Метод Якоби: void SLE::Jacobi(const double&, const Vector&)
3) Метод сопряжённых градиентов: void SLE::SGrd(const double&, const Vector&);
4) Трёхчленная формула реализации метода Ричардсона с чебышёвскими параметрами: void SLE::Rchd3(const double&, const Vector&, const double&, const double&)
# Методы поиска собственных значений - в Matrix.h и Matrix.cpp
1) QR-метод: void Matrix::QR(const double&)
2) Метод обратных итераций со сдвигом с соотношением Рэлея: void Matrix::RQI(const double&, const double&)
# Методы численного анализа - в Polynomial.h и Polynomial.cpp
1) Вывод и ввод полинома осуществляются через операторы << и >> соответственно
2) Значение полинома в точке вычисляется с помощью оператора ()
3) Построение полинома Лагранжа в данных узлах и их значениях: Polynomial Polynomial::int_L(const Vector& - узлы, const Vector& - значения в узлах)
4) Построение полинома Ньютона в данных узлах и их значениях: Polynomial Polynomial::int_N(const Vector& - узлы, const Vector& - значения в узлах)
5) Построение узлов Чебышёва: Vector Cheb(const size_t& - количество необходимых узлов)
