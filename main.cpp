#include <fstream>
#include <iomanip>
#include "SLE.h"
#include "Polynomial.h"
using namespace std;

ifstream fin("input.txt");
ofstream fout("output.txt");

/*
* --------------------------------------------------------������ �������--------------------------------------------------------
* ������� ���� � ����� ����������� �������� � ��������:
* --- ������ ������ ������� ���� (����������� � �������� ���� SLE) ---
* - ����� ������: Gauss()
* - ����� ���������: HR()
* --- ������������ ������ ������� ���� (����������� � �������� ���� SLE) ---
* - ����� �������: HZ(const double& - ��������, const Vector& - ��������� �����������)
* - ����� �����: Jacobi(const double& - ��������, const Vector& - ��������� �����������)
* - ����� ���������� ����������: SGrd(const double& - ��������, const Vector& - ��������� �����������)
* - ���������� ������� ���������� ������ ���������� � ������������ ��������������: Rchd3(const double& - ��������,
* const Vector& - ��������� �����������, const double& - ����� ������� �������, const double& - ������ ������� �������)
* --- ������ ������ ����������� �������� (�������.��� � �������� ���� Matrix)
* - QR-�����: QR(const double& - ��������)
* - ����� �������� �������� �� ������� � ������������ �����: RQI(const double& - ��������, const double& - ����������� �. �.)
* ��������� ������:
* --- ����������������:
* - ���������������� ��������� ��������: int_L(const Vector& - ���� �����, const Vector& - ������ �����)
*/

int main()
{
    fout << fixed << setprecision(8); // ��������� �������� �� 8 ������ ����� �������

    /*//////////////////////////////////////////////////////////
    * //////////////////////////////////////////////////////////
    * //////////////////////////////////////////////////////////
    * ����������������� ������ ��� �������� ���, ���� ���������
    * //////////////////////////////////////////////////////////
    * //////////////////////////////////////////////////////////
    *///////////////////////////////////////////////////////////

    /*
    * ���� ��������� ������
    */

    //int N; fin >> N; // ����������� ������� ��� �������

    ////////////////////////////////////////////////////////////

    //double e; fin >> e; // �������� ���������� (�� ��� ������ �������)
    //Vector v0(N); fin >> v0; // ��������� ����������� (��� ������������ �������)
    //double alpha, beta; fin >> alpha >> beta; // ������� ����������� ����� (��� ������ ����������)
    //double lambda_e; fin >> lambda_e; // ��������� ����������� �. �. (��� ������ �������� ��������)
    //SLE system(N); fin >> system; // ��� ������� ������
    //Matrix A(N); fin >> A; // ��� ������ ����������� ��������

    /*
    * ����� ������� � ������ �� ����������
    */

    // ������ ������
    //fout << system.Gauss(); // ����� ������
    //fout << system.HR(); // ����� ���������

    // ������������ ������
    //system.HZ(e, v0); // ����� �������
    //system.Jacobi(e, v0); // ����� �����
    //system.SGrd(e, v0); // ����� ���������� ����������
    //system.Rchd3(e, v0, alpha, beta); // ����� ����������

    // ������ ������ ����������� �������� � ��������
    //A.QR(e); // QR-�����
    //A.RQI(e, lambda_e); // ����� �������� ��������

    ////////////////////////////////////////////////////////////

    // ���������������� ��������� ��������

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
    fout << "f(" << x0 << ") = " << Lagr(x0) << endl << endl << "f(x) = " << Lagr;
    
    fin.close();
    fout.close();
    return 0;
}