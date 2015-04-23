#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{
	ifstream in(argv[1]);
	ofstream out(argv[2], ios::out | ios::binary);
	
	int c=0;
	while (in.good())
	{
		float val;
		in >> val;
		if (in.good())
		{
			out.write((char*)&val, sizeof(float));
			c++;
		}
	}
	
	in.close();
	out.close();
	
	cout << "Changed " << c << " floats" << endl;
	
	return 0;
}