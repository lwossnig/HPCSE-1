#include <iostream>
#include <complex>

using namespace std;

int main()
{
	complex<double> a(2,1);
	complex<double> i(0,1);
	auto b = i * a;
	cout << a << " " << b << endl;
}

