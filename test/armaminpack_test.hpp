#ifndef ARMAMINPACK_TEST_HPP
#define ARMAMINPACK_TEST_HPP
#include "armaminpack.hpp"

int f01_arma(Col<double> &x, Col<double> &fvec, double *fjac, int *ldfjac,			int &iflag);


int f02_arma(Col<double> &x, Col<double> &fvec,						int &iflag);
int f03_arma(Col<double> &x, Col<double> &fvec, Mat<double> &fjac1, int &iflag);
int f04_arma(Col<double> &x, Col<double> &fvec, Mat<double> &fjac1, int &iflag);
int f05_arma(Col<double> &x, Col<double> &fvec, Mat<double> &fjac1, int &iflag);
int f06_arma(Col<double> &x, Col<double> &fvec,						int &iflag);
int f07_arma(Col<double> &x, Col<double> &fvec,						int &iflag);
int f08_arma(Col<double> &x, Col<double> &fvec, Col<double> &fjrow,	int &iflag);
int f09_arma(Col<double> &x, Col<double> &fvec, Col<double> &fjrow,	int &iflag);


int test01_arma(void);
int test02_arma(void);
int test03_arma(void);
int test04_arma(void);
int test05_arma(void);
int test06_arma(void);
int test07_arma(void);
int test08_arma(void);
int test09_arma(void);

#endif // !MINPACK_TEST_HPP