#include <iostream>
#include <vector>
#include "main.hpp"
using namespace std;

int main()
{
    double b;
    size_t n,m;
    cin >> n >> m >> b;
    Matrix<double> A(n,m);
    cin >> A;
    //>> n >> m;
    //Matrix<double> B(n,m);
    //cin >> B;

    try
    {
        //B.transpose();
        A = A*b;
        //cout << A+(C*B*k);
        //cout << 1 << endl <<A.det() << endl << A.col.size() << " " << A.col.at(0).size() << endl <<A;
        cout << 1 << endl << A << endl;
    }
    catch(char const* i)
    {
        cerr << -1 << endl << i << endl;
    }
    return 0;
}