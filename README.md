# Прямые методы решения СЛАУ
1) Метод Гаусса: Vector SLE::Gauss()
2) Метод отражений: Vector SLE::HR()
# Итерационные методы решения СЛАУ
1) Метод Зейделя: void SLE::HZ(const double&, const Vector&)
2) Метод Якоби: void SLE::Jacobi(const double&, const Vector&)
3) Метод сопряжённых градиентов: void SLE::SGrd(const double&, const Vector&);
4) Трёхчленная формула реализации метода Ричардсона с чебышёвскими параметрами: void SLE::Rchd3(const double&, const Vector&, const double&, const double&)
# Методы поиска собственных значений
1) QR-метод: void Matrix::QR(const double&)
2) Метод обратных итераций со сдвигом с соотношением Рэлея: void Matrix::INVIT(const double&, const double&)
