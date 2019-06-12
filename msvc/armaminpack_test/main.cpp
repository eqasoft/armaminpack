#include "config.hpp"
#include "armaminpack.hpp"
#include "armaminpack_test.hpp"


using namespace std;
using namespace arma;

int main(int argc, char **argv) {

	cout << "      " << endl;
	cout << " ------------------------------------------------------ " << endl;
	cout << " ---------------PERFORMING MINPACK TESTS -------------- " << endl;
	cout << " ------------------------------------------------------ " << endl;
	cout << "      " << endl;
	cout << "      " << endl;

	test01_arma();
	test02_arma();
	test03_arma();
	test04_arma();
	test05_arma();
	test06_arma();
	test07_arma();
	test08_arma();
	test09_arma();

	return 0;
}