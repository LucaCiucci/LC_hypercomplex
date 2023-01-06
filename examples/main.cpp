

#include <iostream>

#include <LC/hyperxomplex/Quaternion.hpp>

using namespace lc;

int main()
{
	std::cout << "Hello World!" << std::endl;

	Quaternion<double> q1(1, 2, 3, 4);

	std::cout << q1 << std::endl;
	std::cout << q1 * q1 << std::endl;

	return 0;
}