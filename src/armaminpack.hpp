#ifndef ARMAMINPACK_HPP
#define ARMAMINPACK_HPP
#include "config.hpp"


typedef int(*U_fp)(...);
typedef int(*S_fp)(...);

// FUNCTION TYPE DECLERATION
#define type_fcn_nn        typedef int (*func_nn)
#define type_fcnder_nn     typedef int (*funcder_nn)
#define type_fcn_mn        typedef int (*func_mn)
#define type_fcnder_mn     typedef int (*funcder_mn)
#define type_fcnderstr_mn  typedef int (*funcderstr_mn)


#define decl_fcn_nn			func_nn			fcn
#define decl_fcnder_nn		funcder_nn		fcn

#define decl_fcnder_mn      funcder_mn		fcn
#define decl_fcn_mn			func_mn			fcn
#define decl_fcnderstr_mn   funcderstr_mn	fcn


/* for hybrd1 and hybrd: */
/*         calculate the functions at x and */
/*         return this vector in fvec. */
/* return a negative value to terminate hybrd1/hybrd */
type_fcn_nn(Col<double> &x, Col<double> &fvec, int &iflag);

/* for hybrj1 and hybrj */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. do not alter fjac. */
/*         if iflag = 2 calculate the jacobian at x and */
/*         return this matrix in fjac. do not alter fvec. */
/* return a negative value to terminate hybrj1/hybrj */
type_fcnder_nn(Col<double> &x1, Col<double> &fvec1, Mat<double> &fjac,int &iflag);

/* for lmder1 and lmder */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. do not alter fjac. */
/*         if iflag = 2 calculate the jacobian at x and */
/*         return this matrix in fjac. do not alter fvec. */
/* return a negative value to terminate lmder1/lmder */
type_fcnder_mn(Col<double> &x1, Col<double> &fvec1, Mat<double> &fjac, int &iflag);

/* for lmdif1 and lmdif */
/*         calculate the functions at x and */
/*         return this vector in fvec. */
/*         if iflag = 1 the result is used to compute the residuals. */
/*         if iflag = 2 the result is used to compute the Jacobian by finite differences. */
/*         Jacobian computation requires exactly n function calls with iflag = 2. */
/* return a negative value to terminate lmdif1/lmdif */
type_fcn_mn(Col<double> &x, Col<double> &fvec, int &iflag);



/* for lmstr1 and lmstr */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. */
/*         if iflag = i calculate the (i-1)-st row of the */
/*         jacobian at x and return this vector in fjrow. */
/* return a negative value to terminate lmstr1/lmstr */
type_fcnderstr_mn(Col<double> &x1, Col<double> &fvec1, Col<double> &fjrow, int &iflag);

// Interface Functions
int hybrd1(decl_fcn_nn, int &n, Col<double> &x1, Col<double> &fvec1, double &tol, int &info);
int hybrj1(decl_fcnder_nn, int &n, Col<double> &x1, Col<double> &fvec1, mat &fjac1, int &ldfjac1,
	double &tol1, int &info1);
int lmder1(decl_fcnder_mn, int &m, int &n, Col<double> &x1, Col<double> &fvec1, mat &fjac1, int &ldfjac1,
	double &tol1, int &info1, Col<int> &ipvt1);
int lmdif1(decl_fcn_mn, int &m, int &n, Col<double> &x1, Col<double> &fvec1, double &tol1, int &info1,
	Col<int> &iwa1, Col<double> &wa1, int  &lwa1);
int lmstr1(funcderstr_mn, int &m, int &n, Col<double> &x1, Col<double> &fvec1, mat &fjac1, int &ldfjac1,
	double &tol1, int &info1, Col<int> &ipvt1, Col<double> &wa1, int &lwa1);


// Background Functions
int fdjac1(S_fp fcn, int &n1, Col<double> &x1, Col<double> &fvec1, double *fjac, int *ldfjac, int *iflag,
	int *ml, int *mu, double *epsfcn, double *wa1, double *wa2);

int fdjac2(S_fp fcn, int &m1, int &n1, Col<double> &x1, Col<double> &fvec1, double *fjac, int *ldfjac, 
	int *iflag, double *epsfcn, double *wa);


int hybrd(decl_fcn_nn, int &n, Col<double> &x1, Col<double> &fvec, double &tol, int &maxfev, int &ml, int &mu,
	double &epsfcn, Col<double> &diag, int &mode, double &factor,
	int &nprint, int &info, int &nfev, Mat<double> &fjac1, int &ldfjac, Col<double> &r, int &lr, Col<double> &qtf);

int hybrj(decl_fcnder_nn, int &n, Col<double> &x1, Col<double> &fvec, Mat<double> &fjac, int &ldfjac, double &xtol,
	int &maxfev, Col<double> &diag1, int &mode, double &factor,	int &nprint, int &info, int &nfev, int &njev, 
	Col<double> &r, int &lr, Col<double> &qtf1);


int lmder(decl_fcnder_mn, int &m, int &n, Col<double> &x1, Col<double> &fvec, Mat<double> &fjac, int &ldfjac, double &ftol,
	double &xtol, double &gtol, int &maxfev, Col<double> &diag1, int &mode, double &factor, int &nprint, int &info,
	int nfev, int njev, Col<int> &ipvt1, Col<double> &qtf1);



int lmdif(decl_fcn_mn, int &m, int &n, Col<double> &x1, Col<double> &fvec, double &ftol, double &xtol, double &gtol,
	int &maxfev, double &epsfcn, Col<double> &diag, int &mode, double &factor, int &nprint, int &info, int &nfev, 
	double *fjac, int &ldfjac, int *ipvt, Col<double> &qtf);


int lmpar(int *n, double *r__, int *ldr, int *ipvt, double *diag, double *qtb, double *delta, double *par, double *x,
	double *sdiag, double *wa1,
	double *wa2);


int lmstr(funcderstr_mn, int &m, int &n, Col<double> &x1, Col<double> &fvec, double *fjac, int &ldfjac, double &ftol, double &xtol,
	double &gtol, int &maxfev, Col<double> diag1,int &mode, double &factor, int &nprint, int &info, int &nfev, int &njev, 
	Col<int> &ipvt1, Col<double> &qtf);


int qform(int *m, int *n, double *q, int *	ldq, double *wa);


int qrfac(int *m, int *n, double *a, int * lda, bool *pivot, int *ipvt, int *lipvt, double *rdiag, double *acnorm, double *wa);


int qrsolv(int *n, double *r__, int *ldr, int *ipvt, double *diag, double *qtb, double *x, double *sdiag, double *wa);


int r1mpyq(int *m, int *n, double *a, int * lda, double *v, double *w);



int r1updt(int *m, int *n, double *s, int * ls, double *u, double *v, double *w, bool *sing);


int rwupdt(int *n, double *r__, int *ldr, double *w, double *b, double *alpha, double *cos__, double *sin__);


int chkder(int *m, int *n, double *x, double *fvec, double *fjac, int *ldfjac, double *xp, double *fvecp, int *mode, double *err);


int dogleg(int *n, double *r__, int *lr, double *diag, double *qtb, double *delta, double *x, double *wa1, double *wa2);

double dpmpar(int id);

double enorm(int *n, double *x);
double enorm(Col<double> &x);




// ARMAMINPACK CLASS DEFINITION


class minpack
{
protected:
	char _type[100];


public:
	func_nn fcn1;
	funcder_nn		fcn2;
	func_mn			fcn3;
	funcder_mn		fcn4;
	funcderstr_mn	fcn5;

	int m, n;

public:
	minpack() {}
	~minpack() {}
	
	char *getType()
	{
		return _type;
	}
};




#endif