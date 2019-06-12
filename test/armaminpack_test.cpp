#include "armaminpack_test.hpp"

// FUNCTIONS

int f01_arma(Col<double> &x1, Col<double> &fvec1, double *fjac, int *ldfjac, int &iflag)
{
	double *x = x1.memptr();
	double *fvec = fvec1.memptr();
	uword n = x1.n_elem;
	int fjac_dim1, fjac_offset;

	static double x_sum__, x_prod__;

	/*  F01 is a function/jacobian routine.

	//  Parameters:

		(1) Inputs -

			int N - the number of variables.

			double X the variable values [N x 1].

			int LDFJAC, the leading dimension of FJAC, which must be at least N.

			int IFLAG:
				1, please compute F(I) (X).
				2, please compute FJAC(I,J) (X).

		(2) Output -

			double FVEC the function values at X if IFLAG = 1.  [N x 1]

			double precision FJAC(LDFJAC,N), the N by N jacobian at X, if IFLAG = 2.  [N x N]

	*/

	--fvec;
	--x;
	fjac_dim1 = *ldfjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;

	/* Function Body */
	x_prod__ = 1.;
	x_sum__ = 0.;

	for (uword i__ = 1; i__ <= n; ++i__) {
		x_prod__ *= x[i__];
		x_sum__ += x[i__];
	}

	/*  If IFLAG is 1, we are supposed to evaluate F(X). */

	if (iflag == 1) {

		for (uword i__ = 1; i__ <= n-1; ++i__) {
			fvec[i__] = x[i__] - (double)(n + 1) + x_sum__;
		}
		fvec[n] = x_prod__ - 1.;

		/*  If IFLAG is 2, we are supposed to evaluate FJAC(I,J) = d F(I)/d X(J) */

	}
	else if (iflag == 2) {

		for (uword j = 1; j <= n; ++j) {
			for (uword i__ = 1; i__ <= n-1; ++i__) {
				fjac[i__ + j * fjac_dim1] = 1.;
			}
		}
		for (uword i__ = 1; i__ <= n-1; ++i__) {
			fjac[i__ + i__ * fjac_dim1] = 2.;
		}
		for (uword j = 1; j <= n; ++j) {
			fjac[n + j * fjac_dim1] = x_prod__ / x[j];
		}
	}
	return 0;
}

int f02_arma(Col<double> &x, Col<double> &fvec, int &iflag)
{

	fvec(0) = x(0) * x(0) - x(0) * 10.0 + x(1) * x(1) + 8.0;
	fvec(1) = x(0) * x(1) * x(1) + x(0) - x(1) * 10.0 + 8.0;
	return 0;
}

int f03_arma(Col<double> &x, Col<double> &fvec, Mat<double> &fjac1, int &iflag)
{
	//double * x = x1.memptr();
	//double * fvec = fvec1.memptr();

	uword ldfjac = fjac1.n_rows;

	double * fjac = fjac1.memptr();

	int fjac_dim1, fjac_offset;

	fjac_dim1 = ldfjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;

	if (iflag == 1) {
		fvec(0) = x(0) * x(0) - x(0) * 10.0 + x(1) * x(1) + 8.0;
		fvec(1) = x(0) * x(1) * x(1) + x(0) - x(1) * 10.0 + 8.0;
	}

	else if (iflag == 2) {
		fjac[fjac_dim1 + 1] = x(0) * 2. - 10.;
		fjac[(fjac_dim1 << 1) + 1] = x(1) * 2.;
		fjac[fjac_dim1 + 2] = x[2] * x(1) + 1.;
		fjac[(fjac_dim1 << 1) + 2] = x(0) * 2. * x(1) - 10.;
	}
	return 0;
}

int f04_arma(Col<double> &x1, Col<double> &fvec1, Mat<double> &fjac1, int &iflag)
{
	uword ldfjac = fjac1.n_rows;
	Col<double> xdat = { 2.,4.,6.,8. };
	Col<double> ydat = { 2.,11.,28.,40. };

	double * fjac = fjac1.memptr();

	int fjac_dim1;
	int fjac_offset;


	fjac_dim1 = ldfjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;

	if (iflag == 1) {	fvec1 = x1(0) * xdat + x1(1) - ydat;}

	else if (iflag == 2) {
		for (uword id = 1; id <= fvec1.n_elem; ++id) {
			fjac[id + fjac_dim1] = xdat(id - 1);
			fjac[id + (fjac_dim1 << 1)] = 1.;
		}
	}
	return 0;
}

int f05_arma(Col<double> &x1, Col<double> &fvec1, Mat<double> &fjac1, int &iflag)
{
	Col<double> xdat = { 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. };
	Col<double> ydat = { 4.,13.,28.,49.,76.,109.,148.,193.,244.,301. };

	uword ldfjac = fjac1.n_rows;

	double * fjac = fjac1.memptr();
	int fjac_dim1, fjac_offset;


	fjac_dim1 = ldfjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;

	if (iflag == 1) {fvec1 = x1(0) + x1(1) * pow(xdat, x1(2)) - ydat;}

	else if (iflag == 2) {
		for (uword id = 1; id <= fvec1.n_elem; ++id) {
			fjac[id + fjac_dim1] = 1.;
			fjac[id + (fjac_dim1 << 1)] = pow(xdat(id - 1), x1(2));
			fjac[id + fjac_dim1 * 3] = x1(1) * log(xdat(id - 1)) * pow(xdat(id - 1), x1(2));
		}
	}
	return 0;
}

int f06_arma(Col<double> &x, Col<double> &fvec, int &iflag)
{

	Col<double> xdat = { 2.,4.,6.,8. };
	Col<double> ydat = { 2.,11.,28.,40. };
	
	fvec = x(0) * xdat + x(1) - ydat;

	return 0;
}

int f07_arma(Col<double> &x, Col<double> &fvec, int &iflag)
{
	Col<double> xdat = { 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. };
	Col<double> ydat = { 4.,13.,28.,49.,76.,109.,148.,193.,244., 301. };

	fvec = x(0) + x(1) * pow(xdat, x(2)) - ydat;
	return 0;
}

int f08_arma(Col<double> &x1, Col<double> &fvec1, Col<double> &fjrow1, int &iflag)
{
	double * fjrow = fjrow1.memptr();

	Col<double> xdat = { 2.,4.,6.,8. };
	Col<double> ydat = { 2.,11.,28.,40. };

	if (iflag == 1) {fvec1= x1(0) * xdat + x1(1) - ydat;}

	else {
		fjrow[0] = xdat[iflag - 2];
		fjrow[1] = 1.;
	}
	return 0;
} 

int f09_arma(Col<double> &x1, Col<double> &fvec1, Col<double> &fjrow1, int &iflag)
{
	double * fjrow = fjrow1.memptr();
	Col<double> xdat = { 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. };
	Col<double> ydat = { 4.,13.,28.,49.,76.,109.,148.,193.,244., 301. };

	if (iflag == 1) {fvec1 = x1(0) + x1(1) * pow(xdat, x1(2)) - ydat;}

	else if (2 <= iflag) {
		fjrow[0] = 1.;
		fjrow[1] = pow(xdat(iflag - 2), x1(2));
		fjrow[2] = x1(1) * log(xdat(iflag - 2)) * pow(xdat(iflag - 2), x1(2));
	}
	return 0;
}


// TESTS

int test01_arma(void)
{

	static int n = 5;
	static int m = n;
	
	static Col<double> x(n);
	static Col<double> xp(n);

	static int ldfjac = n;

	Col<double>  fvec(m);
	Col<double>  fvecp(m);


	Col<double>  err(m);
	Mat<double>  fjac(ldfjac, n);
	
	int ido;
	int iflag;
	int mode;
		
	

	cout << "------------------------------------------------------" << endl;
	cout << "-----------------MINPACK TEST 01----------------------" << endl;
	cout << "  CHKDER compares a user supplied jacobian and a finite" << endl;
	cout << "  difference approximation to it and judges whether the" << endl;
	cout << "  jacobian is correct.?" << endl;
	cout << "-------------------------------------------------------" << endl;

	for (ido = 1; ido <= 2; ++ido) {
		if (ido == 1) {
			cout << " " << endl;
			cout << "  On test 1, use a correct jacobian." << endl;


		}
		else if (ido == 2) {
			cout << " " << endl;
			cout << "  On test 2, use a \"bad\" jacobian" << endl;
			cout << "  and see if the routine notices!" << endl;
		}

		//  Set the point at which the test is to be made:

		for (int id = 1; id <= 5; ++id) {
			x(id - 1) = .5;
		}
		x.print("Evaluation point X :");

		mode = 1;
		chkder(&m, &n, x.memptr(), fvec.memptr(), fjac.memptr(), &ldfjac, xp.memptr(), fvecp.memptr(), &mode, err.memptr());
		iflag = 1;


		f01_arma(x, fvec, fjac.memptr(), &ldfjac, iflag);
		f01_arma(xp, fvecp, fjac.memptr(), &ldfjac, iflag);


		cout << " " << endl;
		cout << "Sampled function values F(X) and F(XP)" << endl;
		cout << " " << endl;

		for (int id = 0; id < m; id++) {
			cout << "     " << fvec(id) << "    |   " << fvecp(id) << endl;
		}

		iflag = 2;
		f01_arma(x, fvec, fjac.memptr(), &ldfjac, iflag);

		/*  Here's where we put a mistake into the jacobian, on purpose. */

		if (ido == 2) {
			fjac(0, 0) *= 1.01;
			fjac(1, 2) = -fjac(1, 2);
		}

		cout << " " << endl;
		cout << "Computed jacobian:" << endl;
		cout << " " << endl;
		fjac.print();

		mode = 2;
		chkder(&m, &n, x.memptr(), fvec.memptr(), fjac.memptr(), &ldfjac, xp.memptr(), fvecp.memptr(), &mode, err.memptr());

		cout << " " << endl;
		cout << "CHKDER gradient error estimates: " << endl;
		cout << " > 0.5, probably correct. " << endl;
		cout << " < 0.5, probably incorrect. " << endl;
		cout << " " << endl;

		cout << "     " << "ID" << "    |   " << "ERR" << endl;

		for (int id = 0; id < 5; id++) {
			cout << "     " << id + 1 << "    |   " << err(id) << endl;
		}
	}
	return 0;
}

int test02_arma(void)
{
	static double tol;
	static int info, iflag;

	int n = 2;
	//int lwa = n * (3 * n + 13) / 2;

	static Col<double> x = Col<double>(n);
	//static Col<double>wa = Col<double>(lwa);
	static Col<double> fvec = Col<double>(n);


	cout << "-------------------------------------------------------" << endl;
	cout << "---------------- MINPACK TEST 02 ----------------------" << endl;
	cout << "  HYBRD1 solves a nonlinear system of equations." << endl;
	cout << "-------------------------------------------------------" << endl;

	x(0) = 3.;
	x(1) = 0.;

	cout << "  " << endl;
	x.print("Initial X :");
	cout << "  " << endl;


	iflag = 1;
	f02_arma(x, fvec, iflag);

	cout << "  " << endl;
	fvec.print("F(X) :");
	cout << "  " << endl;

	tol = 1e-5;

	hybrd1((func_nn)f02_arma, n, x, fvec, tol, info);
	cout << "  " << endl;
	cout << " Returned value of INFO = " << info << endl;
	x.print("Final X :");
	fvec.print("F(X) :");

	return 0;
}

int test03_arma(void)
{
	int n = 2;
	int ldfjac = n;
	
	Col<double> x(n);
	Col<double> fvec(n);

	Mat<double> fjac(ldfjac, n);
	
	int iflag;
	int info;
	double tol;
	
	

	cout << "-------------------------------------------------------" << endl;
	cout << "----------------- MINPACK TEST 03 ---------------------" << endl;
	cout << "  HYBRJ1 solves a nonlinear system of equations." << endl;
	cout << "-------------------------------------------------------" << endl;

	x(0) = 3.;
	x(1) = 0.;

	cout << "  " << endl;
	x.print("Initial X :");
	cout << "  " << endl;

	iflag = 1;
	f02_arma(x, fvec, iflag);

	cout << "  " << endl;
	fvec.print("F(X) :");
	cout << "  " << endl;

	tol = 1e-5;

	hybrj1((funcder_nn)f03_arma, n, x, fvec, fjac, ldfjac, tol, info);

	cout << "  " << endl;
	cout << " Returned value of INFO = " << info << endl;
	x.print("Final X :");
	fvec.print("F(X) :");

	return 0;
}

int test04_arma(void)
{
	int m = 4, n = 2, ldfjac = m, lwa = 5 * n + m;
	int iflag, info;
	double tol;

	Mat<double> fjac(ldfjac, n);
	Col<double> fvec(m);
	Col<double> wa(lwa);
	Col<double> x(n);

	Col<int> ipvt(n);

	cout << "-------------------------------------------------------" << endl;
	cout << "---------------- MINPACK TEST 04 ----------------------" << endl;
	cout << "  LMDER1 minimizes M functions in N variables." << endl;
	cout << "-------------------------------------------------------" << endl;

	x(0) = 0.;
	x(1) = 5.;

	cout << "  " << endl;
	x.print("Initial X :");
	cout << "  " << endl;

	iflag = 1;
	tol = 1e-5;

	f04_arma(x, fvec, fjac,iflag);

	cout << "  " << endl;
	fvec.print("F(X) :");
	cout << "  " << endl;


	lmder1((funcder_mn)f04_arma, m, n, x, fvec, fjac, ldfjac, tol, info, ipvt);

	cout << "  " << endl;
	cout << " Returned value of INFO = " << info << endl;
	x.print("Final X :");
	fvec.print("F(X) :");

	return 0;
}

int test05_arma(void)
{
	int m = 10, n = 3, ldfjac = m, lwa = 5 * n + m;

	Mat<double> fjac(ldfjac, n);
	Col<double> fvec(m);
	int iflag, info;
	Col<int> ipvt(n);
	double tol;
	Col<double> wa(lwa);
	Col<double> x(n);


	cout << "-------------------------------------------------------" << endl;
	cout << "---------------- MINPACK TEST 05 ----------------------" << endl;
	cout << "  LMDER1 minimizes M functions in N variables" << endl;
	cout << "-------------------------------------------------------" << endl;


	x(0) = 0.;
	x(1) = 5.;
	x(2) = 1.3;

	cout << "  " << endl;
	x.print("Initial X :");
	cout << "  " << endl;

	iflag = 1;
	tol = 1e-5;

	f05_arma(x, fvec, fjac,iflag);

	cout << "  " << endl;
	fvec.print("F(X) :");
	cout << "  " << endl;


	lmder1((funcder_mn)f05_arma, m, n, x, fvec, fjac, ldfjac, tol, info, ipvt);

	cout << "  " << endl;
	cout << " Returned value of INFO = " << info << endl;
	x.print("Final X :");
	fvec.print("F(X) :");

	return 0;
}

int test06_arma(void)
{
	int m = 4, n = 2, ldfjac = m, lwa = m * n + 5 * n + m;
	Col<double> fvec(m), wa(lwa), x(n);
	int iflag, info;

	Col<int> iwa(n);


	double tol;

	cout << "-------------------------------------------------------" << endl;
	cout << "------------------ INPACK TEST 06 ---------------------" << endl;
	cout << "  LMDIF1 minimizes M functions in N variables." << endl;
	cout << "-------------------------------------------------------" << endl;

	x(0) = 0.0;
	x(1) = 5.0;

	cout << "  " << endl;
	x.print("Initial X :");
	cout << "  " << endl;

	iflag = 1;
	f06_arma(x, fvec, iflag);

	cout << "  " << endl;
	fvec.print("F(X) :");
	cout << "  " << endl;

	tol = 1e-5;
	lmdif1((func_mn)f06_arma, m, n, x, fvec, tol, info, iwa, wa, lwa);

	cout << "  " << endl;
	cout << " Returned value of INFO = " << info << endl;
	x.print("Final X :");
	fvec.print("F(X) :");

	return 0;
}

int test07_arma(void)
{


	int m = 10, n = 3, lwa = m * n + 5 * n + m;
	Col<double> fvec(m), wa(lwa), x(n);
	int iflag, info;
	Col<int> iwa(n);
	double tol;

	cout << "-------------------------------------------------------" << endl;
	cout << "---------------- MINPACK TEST 07 ----------------------" << endl;
	cout << "  LMDIF1 minimizes M functions in N variables." << endl;
	cout << "-------------------------------------------------------" << endl;

	x(0) = 0.0;
	x(1) = 5.0;
	x(2) = 1.3;

	cout << "  " << endl;
	x.print("Initial X :");
	cout << "  " << endl;

	iflag = 1;
	f07_arma(x, fvec, iflag);

	cout << "  " << endl;
	fvec.print("F(X) :");
	cout << "  " << endl;


	tol = 1e-5;
	lmdif1((func_mn)f07_arma, m, n, x, fvec, tol, info, iwa, wa, lwa);

	cout << "  " << endl;
	cout << " Returned value of INFO = " << info << endl;
	x.print("Final X :");
	fvec.print("F(X) :");

	return 0;
}

int test08_arma(void)
{
	int m = 4, n = 2, ldfjac = m, lwa = 5 * n + m;
	Col<double> fvec(m), wa(lwa), x(n);
	int iflag, info;

	Mat<double> fjac(ldfjac, n);
	Col<double> fjrow(n);


	Col<int> ipvt(n);
	double tol;

	cout << "------------------------------------------------------" << endl;
	cout << "---------------- MINPACK TEST 08 ---------------------" << endl;
	cout << " LMSTR1 minimizes M functions in N variables.'" << endl;
	cout << "------------------------------------------------------" << endl;

	x(0) = 0.0;
	x(1) = 5.0;


	cout << "  " << endl;
	x.print("Initial X :");
	cout << "  " << endl;

	iflag = 1;
	f08_arma(x, fvec, fjrow, iflag);

	cout << "  " << endl;
	fvec.print("F(X) :");
	cout << "  " << endl;

	tol = 1e-5;
	lmstr1((funcderstr_mn)f08_arma, m, n, x, fvec, fjac, ldfjac, tol, info, ipvt, wa, lwa);

	cout << "  " << endl;
	cout << " Returned value of INFO = " << info << endl;
	x.print("Final X :");
	fvec.print("F(X) :");

	return 0;
}

int test09_arma(void)
{
	int m = 10, n = 3, ldfjac = m, lwa = 5 * n + m;
	Col<double> fvec(m), wa(lwa), x(n);
	int iflag, info;

	Mat<double> fjac(ldfjac, n);
	Col<double> fjrow(n);


	Col<int> ipvt(n);
	double tol;

	cout << "-------------------------------------------------------" << endl;
	cout << "----------------- MINPACK TEST 09 ---------------------" << endl;
	cout << " LMSTR1 minimizes M functions in N variables.'" << endl;
	cout << "-------------------------------------------------------" << endl;

	x(0) = 0.0;
	x(1) = 5.0;
	x(2) = 1.3;

	cout << "  " << endl;
	x.print("Initial X :");
	cout << "  " << endl;

	iflag = 1;
	f09_arma(x, fvec, fjrow, iflag);


	cout << "  " << endl;
	fvec.print("F(X) :");
	cout << "  " << endl;



	tol = 1e-5;
	lmstr1((funcderstr_mn)f09_arma, m, n, x, fvec, fjac, ldfjac, tol, info, ipvt, wa, lwa);

	cout << "  " << endl;
	cout << " Returned value of INFO = " << info << endl;
	x.print("Final X :");
	fvec.print("F(X) :");

	return 0;
}


