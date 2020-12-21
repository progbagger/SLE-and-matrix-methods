#include <fstream>
#include <iomanip>
#include "SLE.h"
using namespace std;

ifstream fin("input.txt");
ofstream fout("output.txt");

int main()
{
    fout << fixed << setprecision(8); // установка точности на 8 знаков после запятой
    int n;
    fin >> n;
    double e, a, b;
    Vector x0(n);
    fin >> e >> x0 >> a >> b;
    SLE System(n);
    fin >> System;
    System.Rchd3(e, x0, a, b);
    fin.close();
    fout.close();
    return 0;
}
