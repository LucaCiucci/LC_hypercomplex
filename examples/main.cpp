

#include <iostream>

#include <LC/hyperxomplex/Quaternion.hpp>

using namespace lc;


int main()
{
	std::cout << "Hello World!" << std::endl;

	Quaternion<double> q1(1, 2, 3, 4);

	Quaternion<double> q2(0, 0.01, 0.02, 0.03);

	return 0;
}