#include "armaminpack.hpp"

// Interface Functions

// COMMON INPUTS
// fcn		is the name of the user - supplied subroutine which calculates the functions.fcn must be declared
//			in an external statement in the user calling program, and should be written as follows.
//			fcn(n, x, fvec, iflag)
//				n, iflag
//				x(n), fvec(n)
//			calculate the functions at x and return this vector in fvec.
//			The value of iflag should not be changed by fcn unless
//			the user wants to terminate execution of hybrd1.
//			In this case set iflag to a negative integer.

// n		is a positive integer input variable set to the number of functions and variables.

// x		is an array of length n.on input x must contain an initial estimate of the solution vector. 
//			On output x contains the final estimate of the solution vector.

// fvec		is an output array of length n which contains the functions evaluated at the output x.
//

// tol		is a nonnegative input variable.termination occurs when the algorithm estimates that the relative error
//			between x and the solution is at most tol.
//

//
// wa		is a work array of length lwa.


int hybrd1(decl_fcn_nn, int &n, Col<double> &x, Col<double> &fvec, double &tol, int &info)
{
	// DESCRIPTION
	// The purpose of hybrd1 is to find a zero of a system of n nonlinear functions in n variables by a modification
	// of the powell hybrid method. this is done by using the more general nonlinear equation solver hybrd.the user
	// must provide a subroutine which calculates the functions. The jacobian is then calculated by a forward - difference
	// approximation.


	// INPUT - 


	// OUTPUT - 
	// info		is an integer output variable. if the user has terminated execution, info is set to the(negative) value 
	//			of iflag.see description of fcn.otherwise, info is set as follows.
	//          info = 0   improper input parameters.
	//			info = 1   algorithm estimates that the relative error
	//                    between x and the solution is at most tol.

	//         info = 2   number of calls to fcn has reached or exceeded
	//	                    200 * (n + 1).
	//
	//         info = 3   tol is too small.no further improvement in
	//                    the approximate solution x is possible.	
	//
	//         info = 4   iteration is not making good progress.


	//
	// lwa		is a positive integer input variable not less than
	//						(n*(3 * n + 13)) / 2.


	// Set Initial Variables

	Col<double> diag(n, fill::ones);

	Col<double> qtf(n);				// An output array of length n which contains the vector(q transpose)*fvec.


	static int lr = n * (n + 1) / 2;
	Col<double> r(lr);				// An output array of length lr which contains the upper triangular matrix produced by the qr factorization
									// of the final approximate jacobian, stored rowwise.

	int ldfjac = n;					//a positive integer input variable not less than n which specifies the leading dimension of the array fjac.

	Mat<double> fjac(n, n);			//an output n by n array which contains the orthogonal matrix q produced by the qr factorization of the final approximate jacobian.

	static int nfev = 0;			// an integer output variable set to the number of calls to fcn



	static double factor = 100.;	// a positive input variable used in determining the initial step bound. this bound is set to the product of factor and the 
									// euclidean norm of diag*x if nonzero, or else to factor itself.in most cases factor should lie in the interval(.1, 100.). 100. 
									// is a generally recommended value.

	static int ml = n - 1;			// a nonnegative integer input variable which specifies the number of subdiagonals within the band of the Jacobian matrix. 
									// If the jacobian is not banded, set ml to at least n - 1.

	static int mu = n - 1;			// is a nonnegative integer input variable which specifies the number of superdiagonals within the band of the jacobian matrix. 
									// If the jacobian is not banded, set mu to at least n - 1.

	static int mode = 2;			// An integer input variable. if mode = 1, the variables will be scaled internally. If mode = 2, the scaling is specified by 
									// the input diag. Other values of mode are equivalent to mode = 1.


	info = 0;



	static int maxfev = (n + 1) * 200; // a positive integer input variable. termination occurs when the number of calls to fcn is at least maxfev by the end of an iteration.


	static double epsfcn = 0.0;		// an input variable used in determining a suitable step length for the forward - difference approximation. this  approximation assumes 
									// that the relative errors in the functions are of the order of epsfcn. if epsfcn is less than the machine precision, it is assumed that 
									// the relative errors in the functions are of the order of the machine precision.


	static int nprint = 0;			// n integer input variable that enables controlled printing of iterates if it is positive.in this case, fcn is called with iflag = 0 
									// at the beginning of the first iteration and every nprint iterations thereafter and immediately prior to return, with x and fvec available
									// for printing. if nprint is not positive, no special calls of fcn with iflag = 0 are made.


	// Check error
	if (n <= 0 || tol < 0.0) { return 0; }

	hybrd(fcn, n, x, fvec, tol, maxfev, ml, mu, epsfcn, diag, mode, factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf);

	if (info == 5) { info = 4; }

	return 0;
}


int hybrj1(decl_fcnder_nn, int &n, Col<double> &x1, Col<double> &fvec1, mat &fjac1, int &ldfjac, double &tol, int &info)
{
	// DESCRIPTION
	// The purpose of hybrj1 is to find a zero of a system of n nonlinear functions in n variables by a modification
	// of the powell hybrid method. this is done by using the more general nonlinear equation solver hybrj.the user
	// must provide a subroutine which calculates the functions and the jacobian.


	// INPUT - 

	// fjac		is an output n by n array which contains the orthogonal matrix q produced by the qr factorization
	//			of the final approximate jacobian.

	// ldfjac	is a positive integer input variable not less than n which specifies the leading dimension of the array fjac.


	// info		is an integer output variable. if the user has
	//	        terminated execution, info is set to the(negative)
	//	        value of iflag.see description of fcn.otherwise,
	//	        info is set as follows.
	//	
	//	        info = 0   improper input parameters.
	//	
	//	        info = 1   algorithm estimates that the relative error
	//	                    between x and the solution is at most tol.
	//	
	//	        info = 2   number of calls to fcn with iflag = 1 has
	//	                    reached 100 * (n + 1).
	//	
	//	        info = 3   tol is too small.no further improvement in
	//	                    the approximate solution x is possible.
	//	
	//	        info = 4   iteration is not making good progress.
	//	
	//	
	//	       lwa is a positive integer input variable not less than
	//		c(n*(n + 13)) / 2.


	//int lwa = (n * (n + 13)) / 2;
	//Col<double> wa1(lwa);


	// Set Initial Variables
	//wa1.ones();

	Col<double> diag(n, fill::ones);
	double *fjac = fjac1.memptr();
	//double *wa = wa1.memptr();

	static double factor = 100.;
	int fjac_dim1 = ldfjac;
	int fjac_offset = 1 + fjac_dim1;
	static int lr = n * (n + 1) / 2;

	Col<double> r(lr);				// An output array of length lr which contains the upper triangular matrix produced by the qr factorization
									// of the final approximate jacobian, stored rowwise.

	Col<double> qtf(n);
	static int mode = 2; 
	static int nfev; 
	static int njev;

	static double xtol;
	static int maxfev = (n + 1) * 100;
	static int nprint = 0;

	fjac -= fjac_offset;
	info = 0;

	if (n <= 0 || ldfjac < n || tol < 0.0) { return 0; }

	xtol = tol;

	hybrj(fcn, n, x1, fvec1, fjac1, ldfjac, xtol, maxfev, diag, mode, factor, nprint, info, nfev, njev, r, lr, qtf);

	if (info == 5) {
		info = 4;
	}
	return 0;
}


int lmder1(decl_fcnder_mn, int &m, int &n, Col<double> &x1, Col<double> &fvec1, mat &fjac1, int &ldfjac, double &tol, int &info, Col<int> &ipvt1)
{
	// DESCRIPTION

	// The purpose of lmder1 is to minimize the sum of the squares of m nonlinear functions in n variables by a modification of the
	// levenberg - marquardt algorithm. this is done by using the more general least - squares solver lmder. The user must provide a
	// subroutine which calculates the functions and the jacobian.


	// fjac		is an output m by n array.the upper n by n submatrix
	//         of fjac contains an upper triangular matrix r with
	//	c         diagonal elements of nonincreasing magnitude such that
	//	c
	//	c                t     t           t
	//	c               p *(jac *jac)*p = r * r,
	//	c
	//	c         where p is a permutation matrix and jac is the final
	//	c         calculated jacobian.column j of p is column ipvt(j)
	//	c(see below) of the identity matrix.the lower trapezoidal
	//	c         part of fjac contains information generated during
	//	c         the computation of r.
	//	c
	//	c       ldfjac is a positive integer input variable not less than m
	//	c         which specifies the leading dimension of the array fjac.


	// ipvt		is an integer output array of length n.ipvt defines a permutation matrix p such that jac*p = q * r,
	//			where jac is the final calculated jacobian, q is orthogonal(not stored), and r is upper triangular
	//			with diagonal elements of nonincreasing magnitude. Column j of p is column ipvt(j) of the identity matrix.

	// wa is a positive integer input variable not less than 5 * n + m.


	// info is an integer output variable. if the user has
	//	         terminated execution, info is set to the(negative)
	//	         value of iflag.see description of fcn.otherwise,
	//	         info is set as follows.
	//	
	//	         info = 0  improper input parameters.
	//	
	//	         info = 1  algorithm estimates that the relative error
	//	                   in the sum of squares is at most tol.
	//	
	//	         info = 2  algorithm estimates that the relative error
	//	                   between x and the solution is at most tol.
	//	
	//	         info = 3  conditions for info = 1 and info = 2 both hold.
	//	
	//	         info = 4  fvec is orthogonal to the columns of the
	//	                   jacobian to machine precision.
	//	
	//	         info = 5  number of calls to fcn with iflag = 1 has
	//	                   reached 100 * (n + 1).
	//	
	//	         info = 6  tol is too small.no further reduction in
	//	                   the sum of squares is possible.
	//	
	//	         info = 7  tol is too small.no further improvement in
	//	                   the approximate solution x is possible.




	//double *x = x1.memptr();
	//double *fvec = fvec1.memptr();
	double *fjac = fjac1.memptr();
	//double * wa = wa1.memptr();
	int * ipvt = ipvt1.memptr();

	static double factor = 100.;
	int fjac_dim1 = ldfjac;
	int	fjac_offset = 1 + fjac_dim1;
	static int mode = 1;
	static int nfev;
	static int njev;

	static double ftol = tol;
	static double gtol = 0.0;
	static double xtol = tol;
	static int maxfev = (n + 1) * 100;
	static int nprint = 0;

	Col<double> diag(n, fill::ones);
	Col<double> qtf(n);

	fjac -= fjac_offset;

	info = 0;

	if (n <= 0 || m < n || ldfjac < m || tol < 0.0) { return 0; }

	lmder(fcn, m, n, x1, fvec1, fjac1, ldfjac, ftol, xtol, gtol, maxfev, diag, mode, factor, nprint,	info, nfev, njev, ipvt1, qtf);

	


	if (info == 8) { info = 4; }
	return 0;
}


int lmdif1(decl_fcn_mn, int &m, int &n, Col<double> &x1, Col<double> &fvec1, double &tol, int &info, Col<int> &iwa1, Col<double> &wa1, int  &lwa)
{
	//DESCRIPTION
	// The purpose of lmdif1 is to minimize the sum of the squares of m nonlinear functions in n variables by a modification of the
	// levenberg - marquardt algorithm. this is done by using the more general least - squares solver lmdif.the user must provide a
	// subroutine which calculates the functions.the jacobian is then calculated by a forward - difference approximation.

	//  info is an integer output variable. if the user has
	//	c         terminated execution, info is set to the(negative)
	//	c         value of iflag.see description of fcn.otherwise,
	//	c         info is set as follows.
	//	c
	//	c         info = 0  improper input parameters.
	//	c
	//	c         info = 1  algorithm estimates that the relative error
	//	c                   in the sum of squares is at most tol.
	//	c
	//	c         info = 2  algorithm estimates that the relative error
	//	c                   between x and the solution is at most tol.
	//	c
	//	c         info = 3  conditions for info = 1 and info = 2 both hold.
	//	c
	//	c         info = 4  fvec is orthogonal to the columns of the
	//	c                   jacobian to machine precision.
	//	c
	//	c         info = 5  number of calls to fcn has reached or
	//	c                   exceeded 200 * (n + 1).
	//	c
	//	c         info = 6  tol is too small.no further reduction in
	//	c                   the sum of squares is possible.
	//	c
	//	c         info = 7  tol is too small.no further improvement in
	//	c                   the approximate solution x is possible.
	//	c
	//	c       iwa is an integer work array of length n.
	//	c
	//	c       wa is a work array of length lwa.
	//	c
	//	c       lwa is a positive integer input variable not less than
	//	c         m*n + 5 * n + m.

	double *x = x1.memptr();
	double *fvec = fvec1.memptr();
	double * wa = wa1.memptr();
	int * iwa = iwa1.memptr();

	static double factor = 100.;
	static int mp5n = m + n * 5;
	static int	mode = 1;
	static int nfev;
	static double ftol = tol;
	static double gtol = 0.0;
	static double xtol = tol;
	static double epsfcn = 0.0;
	static int maxfev = (n + 1) * 200;
	static int nprint = 0;

	Col<double> diag(n, fill::ones);
	Col<double> qtf(n);

	int ldfjac = m;

	info = 0;

	if (n <= 0 || m < n || tol < 0.0 || lwa < m * n + n * 5 + m) { return 0; }

	lmdif(fcn, m, n, x1, fvec1, ftol, xtol, gtol, maxfev, epsfcn, diag, mode,
		factor, nprint, info, nfev, &wa[mp5n + 0], ldfjac, &iwa[0], qtf);
	
	if (info == 8) { info = 4; }
	return 0;
}


int lmstr1(funcderstr_mn fcn, int &m, int &n, Col<double> &x1, Col<double> &fvec1, mat &fjac1, int &ldfjac, double &tol, int &info, Col<int> &ipvt1, Col<double> &wa1, int &lwa)
{
	// DESCRIPTION - 
	// The purpose of lmstr1 is to minimize the sum of the squares of m nonlinear functions in n variables by a modification of
	// the levenberg - marquardt algorithm which uses minimal storage. This is done by using the more general least - squares solver
	// lmstr.the user must provide a subroutine which calculates the functions and the rows of the jacobian.

	//       m is a positive integer input variable set to the number
	//	         of functions.
	//	
	//	       n is a positive integer input variable set to the number
	//	         of variables.n must not exceed m.
	//	
	//       fjac is an output n by n array.the upper triangle of fjac
	//	         contains an upper triangular matrix r such that
	//	
	//	                t     t           t
	//	               p *(jac *jac)*p = r * r,
	//	
	//	         where p is a permutation matrix and jac is the final
	//	         calculated jacobian.column j of p is column ipvt(j)
	//	(see below) of the identity matrix.the lower triangular
	//	         part of fjac contains information generated during
	//	         the computation of r.

	//       info is an integer output variable. if the user has
	//	         terminated execution, info is set to the(negative)
	//	         value of iflag.see description of fcn.otherwise,
	//	         info is set as follows.
	//	
	//	         info = 0  improper input parameters.
	//	
	//	         info = 1  algorithm estimates that the relative error
	//	                   in the sum of squares is at most tol.
	//	
	//	         info = 2  algorithm estimates that the relative error
	//	                   between x and the solution is at most tol.
	//	
	//	         info = 3  conditions for info = 1 and info = 2 both hold.
	//	
	//	         info = 4  fvec is orthogonal to the columns of the
	//	                   jacobian to machine precision.
	//	
	//	         info = 5  number of calls to fcn with iflag = 1 has
	//	                   reached 100 * (n + 1).
	//	
	//	         info = 6  tol is too small.no further reduction in
	//	                   the sum of squares is possible.
	//	
	//	         info = 7  tol is too small.no further improvement in
	//	                   the approximate solution x is possible.
	//	
	//	       ipvt is an integer output array of length n.ipvt
	//	         defines a permutation matrix p such that jac*p = q * r,
	//	         where jac is the final calculated jacobian, q is
	//	         orthogonal(not stored), and r is upper triangular.
	//	         column j of p is column ipvt(j) of the identity matrix.
	//	
	//	       wa is a work array of length lwa.
	//	
	//	       lwa is a positive integer input variable not less than 5 * n + m.

	double *x = x1.memptr();
	double *fvec = fvec1.memptr();
	double *fjac = fjac1.memptr();
	int * ipvt = ipvt1.memptr();
	double *wa = wa1.memptr();

	static double factor = 100.;
	int fjac_dim1 = ldfjac;
	int fjac_offset = 1 + fjac_dim1;
	static int mode = 1;
	static int  nfev;
	static int  njev;
	static double ftol = tol;
	static double gtol = 0.0;
	static double xtol = tol;
	static int maxfev = (n + 1) * 100;
	static int nprint = 0;

	Col<double> diag(n, fill::ones);
	Col<double> qtf(n);

	fjac -= fjac_offset;
	info = 0;

	if (n <= 0 || m < n || ldfjac < n || tol < 0.0 || lwa < n * 5 + m) { return 0; }


	lmstr(fcn, m, n, x1, fvec1, &fjac[fjac_offset], ldfjac, ftol, xtol, gtol, maxfev, diag, mode, factor, nprint,	info, nfev, njev, ipvt1, qtf);

	if (info == 8) {info = 4;}

	return 0;
}


// Background functions

int hybrd(decl_fcn_nn, int &n, Col<double> &x, Col<double> &fvec, double &tol, int &maxfev, int &ml, int &mu, double &epsfcn, Col<double> &diag,
	int &mode, double &factor,int &nprint, int &info, int &nfev, Mat<double> &fjac1, int &ldfjac, Col<double> &r, int &lr, Col<double> &qtf)
{
	int c1 = 1;
	bool ff = false;

	double * fjac = fjac1.memptr();
	int fjac_dim1 = ldfjac;
	int fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;

	static int l;
	static int jm1;
	static int iwa[1];
	static double sum;
	static bool sing;
	static int iter;
	static double temp;
	static int msum;
	static int iflag;
	static double delta;
	static bool jeval;
	static int ncsuc;
	static double ratio;
	static double fnorm;
	static double pnorm;
	static double xnorm;
	static double fnorm1;
	static int nslow1;
	static int nslow2;
	static int ncfail;
	static double actred;
	static double prered;

	static double epsmch = dpmpar(1);

	Col<double> w1(n), w2(n), w3(n), w4(n);  // Working vectors

	info = 0;
	iflag = 0;
	
	// check the input parameters for errors.
	if (n <= 0 || tol < 0.0 || maxfev <= 0 || ml < 0 || mu < 0 || 
		factor <= 0.0 || ldfjac < n || lr < n * (n + 1) / 2) {
		goto L300;
	}
	
	if (mode == 1 && any(diag <= 0.0)) {
		goto L300;
	}

	// Evaluate the function at the starting point and calculate its norm.

	iflag = 1;
	(*fcn)(x, fvec, iflag);++nfev;
	if (iflag < 0) {goto L300;}

	fnorm = enorm(fvec);

	// Determine the number of calls to fcn needed to compute the jacobian matrix.
	msum = min(ml + mu + 1, n);

	// initialize iteration counter and monitors. 
	iter = 1;
	ncsuc = 0;
	ncfail = 0;
	nslow1 = 0;
	nslow2 = 0;

	// Beginning of the outer loop
L30:
	jeval = true;

	// Calculate the jacobian matrix

	iflag = 2;
	fdjac1((S_fp)fcn, n, x, fvec, &fjac[fjac_offset], &ldfjac, &iflag, &ml, &mu, &epsfcn, w1.memptr(), w2.memptr());
	nfev += msum;
	if (iflag < 0) {goto L300;}
	
	// compute the qr factorization of the jacobian.
	qrfac(&n, &n, &fjac[fjac_offset], &ldfjac, &ff, iwa, &c1, w1.memptr(), w2.memptr(), w3.memptr());

	// On the first iteration and if mode is 1, scale according
	// to the norms of the columns of the initial jacobian.

	if (iter == 1) {
		if (mode == 1) {
			diag = w2;
			diag.elem(find(w2 == 0)).ones();
		}

		// on the first iteration, calculate the norm of the scaled x
		// and initialize the step bound delta.
		w3 = diag % x;
		xnorm = enorm(w3);
		delta = factor * xnorm;
		if (delta == 0.0) {delta = factor;}
	}

	// form (q transpose)*fvec and store in qtf.

	qtf = fvec;

	for (int jd = 1; jd <= n; ++jd) {

		if (fjac[jd + jd * fjac_dim1] != 0.0) {
			sum = 0.0;
			for (int id = jd; id <= n; ++id) {
				sum += fjac[id + jd * fjac_dim1] * qtf(id - 1);
			}
			temp = -sum / fjac[jd + jd * fjac_dim1];
			for (int id = jd; id <= n; ++id) {
				qtf(id - 1) += fjac[id + jd * fjac_dim1] * temp;
			}
		}
	}

	// copy the triangular factor of the qr factorization into r.
	sing = false;
	for (int j = 1; j <= n; ++j) {
		l = j;
		jm1 = j - 1;
		if (jm1 >= 1) {
			for (int id = 1; id <= jm1; ++id) {
				r(l - 1) = fjac[id + j * fjac_dim1];
				l = l + n - id;
			}
		}
		r(l - 1) = w1(j - 1);
		if (w1(j - 1) == 0.0) {
			sing = true;
		}
	}

	// accumulate the orthogonal factor in fjac.

	qform(&n, &n, &fjac[fjac_offset], &ldfjac, w1.memptr());

	// rescale if necessary.
	if (mode == 1) {diag = arma::max(diag, w2);}
	
	// beginning of the inner loop.
L180:

	// if requested, call fcn to enable printing of iterates.
	if (nprint > 0) {
		iflag = 0;
		if ((iter - 1) % nprint == 0) {
			(*fcn)(x, fvec, iflag);
		}
		if (iflag < 0) {goto L300;}
	}

	// determine the direction p.
	dogleg(&n, r.memptr(), &lr, diag.memptr(), qtf.memptr(), &delta, w1.memptr(), w2.memptr(), w3.memptr());

	// store the direction p and x + p. calculate the norm of p.
	w1 = -w1;
	w2 = x + w1;
	w3 = diag % w1;
	pnorm = enorm(w3);

	// On the first iteration, adjust the initial step bound. */
	if (iter == 1) {delta = min(delta, pnorm);}

	// evaluate the function at x + p and calculate its norm. */
	iflag = 1;
	(*fcn)(w2, w4, iflag);++(nfev);

	if (iflag < 0) {goto L300;}

	fnorm1 = enorm(w4);

	// compute the scaled actual reduction. */
	actred = -1.0;
	if (fnorm1 < fnorm) {
		actred = 1.0 - pow(fnorm1 / fnorm, 2);
	}

	// compute the scaled predicted reduction.
	l = 1;
	for (int id = 1; id <= n; ++id) {
		sum = 0.0;
		for (int j = id; j <= n; ++j) {
			sum += r(l - 1) * w1(j - 1);
			++l;
		}
		w3(id - 1) = qtf(id - 1) + sum;
	}

	temp = enorm(w3);

	prered = 0.0;
	if (temp < fnorm) {	prered = 1.0 - pow(temp / fnorm, 2);}

	// compute the ratio of the actual to the predicted reduction.
	ratio = 0.0;
	if (prered > 0.0) {
		ratio = actred / prered;
	}

	// update the step bound.
	if (ratio >= 0.1) {
		ncfail = 0;
		++ncsuc;
		if (ratio >= 0.5 || ncsuc > 1) {
			delta = max(delta, 2 * pnorm);
		}
		if (abs(ratio - 1.0) <= 0.1) {
			delta = 2 * pnorm;
		}
	}

	ncsuc = 0;
	++ncfail;
	delta = 0.5* delta;

	// Test for successful iteration.

	if (ratio >= 0.0001) {	// Successful iteration. update x, fvec, and their norms. */
		x = w2;
		w2 = diag % x;
		fvec = w4;
		xnorm = enorm(w2);
		fnorm = fnorm1;
		++iter;
	}

	// determine the progress of the iteration.
	++nslow1;
	if (actred >= 0.001) {nslow1 = 0;}

	if (jeval) {++nslow2;}

	if (actred >= 0.1) {nslow2 = 0;}
	   
	// test for convergence
	if (delta <= tol * xnorm || fnorm == 0.0) {	info = 1;}

	if (info != 0) {goto L300;}

	// Tests for termination and stringent tolerances. */
	if (nfev >= maxfev) {info = 2;}

	/* Computing MAX */
	if (0.1 * max(0.1 * delta, pnorm) <= epsmch * xnorm) {info = 3;}

	if (nslow2 == 5) {info = 4;}

	if (nslow1 == 10) {info = 5;}

	if (info != 0) {goto L300;}


	// criterion for recalculating jacobian approximation by forward differences.
	if (ncfail == 2) {
		goto L290;
	}

	// calculate the rank one modification to the jacobian and update qtf if necessary

	for (int j = 1; j <= n; ++j) {
		sum = 0.0;
		for (int id = 1; id <= n; ++id) {
			sum += fjac[id + j * fjac_dim1] * w4(id - 1);
		}
		w2(j - 1) = (sum - w3(j - 1)) / pnorm;
		w1(j - 1) = diag(j - 1) * (diag(j - 1) * w1(j - 1) / pnorm);
		if (ratio >= 0.0001) {
			qtf(j - 1) = sum;
		}
	}

	// compute the qr factorization of the updated jacobian.
	r1updt(&n, &n, r.memptr(), &lr, w1.memptr(), w2.memptr(), w3.memptr(), &sing);
	r1mpyq(&n, &n, &fjac[fjac_offset], &ldfjac, w2.memptr(), w3.memptr());
	r1mpyq(&c1, &n, qtf.memptr(), &c1, w2.memptr(), w3.memptr());

	// End of the inner loop. */
	jeval = false;
	goto L180;

L290:

	// End of the outer loop. */
	goto L30;
L300:

	// Termination, either normal or user imposed. */
	if (iflag < 0) {
		info = iflag;
	}
	iflag = 0;
	if (nprint > 0) { (*fcn)(x, fvec, iflag); }

	return 0;
}


int hybrj(decl_fcnder_nn, int &n, Col<double> &x, Col<double> &fvec, Mat<double> &fjac1, int &ldfjac, double &xtol, int &maxfev, Col<double> &diag, int &mode, double &factor,
	int &nprint, int &info, int &nfev, int &njev, Col<double> &r, int &lr, Col<double> &qtf)
{
	int c1 = 1;
	bool ff = false;
	double *fjac = fjac1.memptr();

	int fjac_dim1, fjac_offset;
	double d__1;

	static int i__, j, l, jm1, iwa[1];
	static double sum;
	static bool sing;
	static int iter;
	static double temp;
	static int iflag;
	static double delta;
	static bool jeval;
	static int ncsuc;
	static double ratio;
	static double fnorm;
	static double pnorm;
	static double xnorm;
	static double fnorm1;
	static int nslow1;
	static int nslow2;
	static int ncfail;
	static double actred;
	static double  epsmch;
	static double  prered;

	Col<double> w1(n), w2(n), w3(n), w4(n);  // Working vectors

	fjac_dim1 = ldfjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;

	// Epsmch is the machine precision.
	epsmch = dpmpar(1);

	info = 0;
	iflag = 0;
	nfev = 0;
	njev = 0;

	// Check the input parameters for errors.

	if (n <= 0 || ldfjac < n || xtol < 0.0 || maxfev <= 0 || factor <= 0.0 || lr < n * (n + 1) / 2) {
		goto L300;
	}

	if (mode == 1 && any(diag <= 0.0)) {
		goto L300;
	}

	//  evaluate the function at the starting point and calculate its norm.

	iflag = 1;
	(*fcn)(x, fvec, fjac1, iflag);++nfev;
	if (iflag < 0) {goto L300;}
	fnorm = enorm(fvec);

	// Initialize iteration counter and monitors.
	iter = 1; 
	ncsuc = 0;
	ncfail = 0;
	nslow1 = 0;
	nslow2 = 0;

	// Beginning of the outer loop.
L30:
	jeval = true;

	// Calculate the jacobian matrix.

	iflag = 2;
	(*fcn)(x, fvec, fjac1, iflag); ++(njev);
	if (iflag < 0) {
		goto L300;
	}

	// Compute the qr factorization of the jacobian.
	qrfac(&n, &n, &fjac[fjac_offset], &ldfjac, &ff, iwa, &c1, w1.memptr(), w2.memptr(), w3.memptr());

	// On the first iteration and if mode is 1, scale according
	// to the norms of the columns of the initial jacobian.
	if (iter == 1) {
		if (mode == 1) {
			diag = w2;
			diag.elem(find(w2 == 0)).ones();
		}

		// on the first iteration, calculate the norm of the scaled x
		// and initialize the step bound delta.
		w3 = diag % x;
		xnorm = enorm(w3);
		delta = factor * xnorm;
		if (delta == 0.0) { delta = factor; }
	}

	// form (q transpose)*fvec and store in qtf.
	qtf = fvec;

	for (int jd = 1; jd <= n; ++jd) {

		if (fjac[jd + jd * fjac_dim1] != 0.0) {
			sum = 0.0;
			for (int id = jd; id <= n; ++id) {
				sum += fjac[id + jd * fjac_dim1] * qtf(id - 1);
			}
			temp = -sum / fjac[jd + jd * fjac_dim1];
			for (int id = jd; id <= n; ++id) {
				qtf(id - 1) += fjac[id + jd * fjac_dim1] * temp;
			}
		}
	}

	// copy the triangular factor of the qr factorization into r.
	sing = false;
	for (int j = 1; j <= n; ++j) {
		l = j;
		jm1 = j - 1;
		if (jm1 >= 1) {
			for (int id = 1; id <= jm1; ++id) {
				r(l - 1) = fjac[id + j * fjac_dim1];
				l = l + n - id;
			}
		}
		r(l - 1) = w1(j - 1);
		if (w1(j - 1) == 0.0) {
			sing = true;
		}
	}

	// Accumulate the orthogonal factor in fjac.
	qform(&n, &n, &fjac[fjac_offset], &ldfjac, w1.memptr());

	// rescale if necessary.
	if (mode == 1) {diag = arma::max(diag, w2);}

	// beginning of the inner loop.

L180:

	// If requested, call fcn to enable printing of iterates.

	if (nprint > 0) {
		iflag = 0;
		if ((iter - 1) % nprint == 0) {
			(*fcn)(x, fvec, fjac1, iflag);
		}
		if (iflag < 0) {goto L300;}
	}

	// Determine the direction p.

	dogleg(&n, r.memptr(), &lr, diag.memptr(), qtf.memptr(), &delta, w1.memptr(), w2.memptr(), w3.memptr());

	// store the direction p and x + p. calculate the norm of p.
	w1 = -w1;
	w2 = x + w1;
	w3 = diag % w1;
	pnorm = enorm(w3);

	// on the first iteration, adjust the initial step bound. */
	if (iter == 1) {delta = min(delta, pnorm);}

	// evaluate the function at x + p and calculate its norm. */
	iflag = 1;
	(*fcn)(w2, w4, fjac1, iflag); ++(nfev);
	if (iflag < 0) {goto L300;}
	fnorm1 = enorm(w4);

	// compute the scaled actual reduction. */
	actred = -1.0;
	if (fnorm1 < fnorm) {
		/* Computing 2nd power */
		d__1 = fnorm1 / fnorm;
		actred = 1.0 - d__1 * d__1;
	}

	/*           compute the scaled predicted reduction. */

	l = 1;

	for (i__ = 1; i__ <= n; ++i__) {
		sum = 0.0;

		for (j = i__; j <= n; ++j) {
			sum += r(l-1) * w1(j - 1);
			++l;
		}
		w3[i__-1] = qtf[i__-1] + sum;
	}
	temp = enorm(w3);
	prered = 0.0;
	if (temp < fnorm) {
		/* Computing 2nd power */
		d__1 = temp / fnorm;
		prered = 1.0 - d__1 * d__1;
	}

	/*           compute the ratio of the actual to the predicted */
	/*           reduction. */

	ratio = 0.0;
	if (prered > 0.0) {
		ratio = actred / prered;
	}

	// update the step bound.
	if (ratio >= 0.1) {
		ncfail = 0;
		++ncsuc;
		if (ratio >= 0.5 || ncsuc > 1) {
			delta = max(delta, 2 * pnorm);
		}
		if (abs(ratio - 1.0) <= 0.1) {
			delta = 2 * pnorm;
		}
	}

	ncsuc = 0;
	++ncfail;
	delta = 0.5* delta;

	// test for successful iteration. 
	if (ratio >= 0.0001) {	// Successful iteration. update x, fvec, and their norms. */
		x = w2;
		w2 = diag % x;
		fvec = w4;
		xnorm = enorm(w2);
		fnorm = fnorm1;
		++iter;
	}

	
	/*           determine the progress of the iteration. */

	++nslow1;
	if (actred >= 0.001) {
		nslow1 = 0;
	}
	if (jeval) {
		++nslow2;
	}
	if (actred >= 0.1) {
		nslow2 = 0;
	}

	// Test for convergence. */

	if (delta <= xtol * xnorm || fnorm == 0.0) {info = 1;}
	if (info != 0) {
		goto L300;
	}

	/*           tests for termination and stringent tolerances. */

	if (nfev >= maxfev) {info = 2;}

	/* Computing MAX */
	d__1 = 0.1 * delta;
	if (0.1 * max(d__1, pnorm) <= epsmch * xnorm) {
		info = 3;
	}
	if (nslow2 == 5) {
		info = 4;
	}
	if (nslow1 == 10) {
		info = 5;
	}
	if (info != 0) {
		goto L300;
	}

	/*           criterion for recalculating jacobian. */

	if (ncfail == 2) {
		goto L290;
	}

	/*           calculate the rank 1.0 modification to the jacobian */
	/*           and update qtf if necessary. */

	for (j = 1; j <= n; ++j) {
		sum = 0.0;
		for (i__ = 1; i__ <= n; ++i__) {
			sum += fjac[i__ + j * fjac_dim1] * w4(i__-1);
		}
		w2(j - 1) = (sum - w3(j - 1)) / pnorm;
		w1(j - 1) = diag(j-1) * (diag(j-1) * w1(j - 1) / pnorm);
		if (ratio >= 0.0001) {
			qtf[j-1] = sum;
		}
	}

	/*           compute the qr factorization of the updated jacobian. */

	r1updt(&n, &n, r.memptr(), &lr, w1.memptr(), w2.memptr(), w3.memptr(), &sing);
	r1mpyq(&n, &n, &fjac[fjac_offset], &ldfjac, w2.memptr(), w3.memptr());
	r1mpyq(&c1, &n, qtf.memptr(), &c1, w2.memptr(), w3.memptr());

	// End of the inner loop

	jeval = false;
	goto L180;
L290:
	// End of the outer loop
	goto L30;
L300:

	/*     termination, either normal or user imposed. */

	if (iflag < 0) {
		info = iflag;
	}
	iflag = 0;
	if (nprint > 0) {
		(*fcn)(x, fvec, fjac1, iflag);
	}
	return 0;
}

int lmder(decl_fcnder_mn, int &m, int &n, Col<double> &x, Col<double> &fvec, Mat<double> &fjac1, int &ldfjac, double &ftol, double &xtol, double &gtol, int &maxfev,
	Col<double> &diag, int &mode, double &factor, int &nprint, int &info, int nfev, int njev, Col<int> &ipvt, Col<double> &qtf)
{
	double * fjac = fjac1.memptr();
	bool ft = true;
	int fjac_dim1, fjac_offset;
	static int i__, j, l;
	static double par, sum;
	static int iter;
	static double temp, temp1, temp2;
	static int iflag;
	static double delta;
	static double ratio;
	static double fnorm, gnorm, pnorm, xnorm, fnorm1, actred, dirder, epsmch, prered;



	Col<double> w1(n), w2(n), w3(n), w4(m);  // Working vectors

	fjac_dim1 = ldfjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;

	// epsmch is the machine precision.
	epsmch = dpmpar(1);

	info = 0;
	iflag = 0;
	nfev = 0;
	njev = 0;

	// check the input parameters for errors. */

	if (n <= 0 || m < n || ldfjac < m || ftol < 0.0 || xtol < 0.0 || gtol < 0.0 || maxfev <= 0 || factor <= 0.0) {
		goto L300;
	}

	if (mode == 1 && any(diag <= 0.0)) {
		goto L300;
	}


	// evaluate the function at the starting point
	// and calculate its norm.

	iflag = 1;
	(*fcn)(x, fvec, fjac1, iflag);++nfev;
	if (iflag < 0) {goto L300;}
	fnorm = enorm(fvec);

	//  initialize levenberg-marquardt parameter and iteration counter.
	par = 0.0;
	iter = 1;

	// Beginning of the outer loop.
L30:

	// calculate the jacobian matrix.
	iflag = 2;
	(*fcn)(x, fvec, fjac1,iflag); ++(njev);
	if (iflag < 0) { goto L300; }

	// if requested, call fcn to enable printing of iterates.
	if (nprint > 0) {
		iflag = 0;
		if ((iter - 1) % nprint == 0) {
			(*fcn)(x, fvec, fjac1,iflag);
		}
		if (iflag < 0) {goto L300;}
	}


	// compute the qr factorization of the jacobian.
	qrfac(&m, &n, &fjac[fjac_offset], &ldfjac, &ft, ipvt.memptr(), &n, w1.memptr(), w2.memptr(), w3.memptr());

	// On the first iteration and if mode is 1, scale according to the norms of the columns of the initial jacobian. */

	if (iter == 1) {
		if (mode == 1) {
			diag = w2;
			diag.elem(find(w2 == 0)).ones();
		}

		// On the first iteration, calculate the norm of the scaled x and initialize the step bound delta. */
		w3= diag % x;
		xnorm = enorm(w3);
		delta = factor * xnorm;
		if (delta == 0.0) {	delta = factor;}
	}


	//form (q transpose)*fvec and store the first n com0.1.0nts in qtf.
	w4 = fvec;
	for (j = 1; j <= n; ++j) {
		if (fjac[j + j * fjac_dim1] != 0.0) {
			sum = 0.0;
			for (i__ = j; i__ <= m; ++i__) {
				sum += fjac[i__ + j * fjac_dim1] * w4(i__ - 1);
			}
			temp = -sum / fjac[j + j * fjac_dim1];
			for (i__ = j; i__ <= m; ++i__) {
				w4(i__ - 1) += fjac[i__ + j * fjac_dim1] * temp;
			}
		}
		fjac[j + j * fjac_dim1] = w1(j-1);
		qtf(j - 1) = w4(j-1);
	}


	// compute the norm of the scaled gradient.

	gnorm = 0.0;
	if (fnorm != 0.0) {

		for (j = 1; j <= n; ++j) {
			l = ipvt(j - 1);
			if (w2(l - 1) != 0.0) {
				sum = 0.0;
				for (i__ = 1; i__ <= j; ++i__) {
					sum += fjac[i__ + j * fjac_dim1] * (qtf(i__ - 1) / fnorm);
				}
				gnorm = max(gnorm, abs(sum / w2(l - 1)));
			}
		}
	}

	// Test for convergence of the gradient norm.

	if (gnorm <= gtol) { info = 4; }

	if (info != 0) { goto L300; }

	// rescale if necessary.
	if (mode == 1) {diag = arma::max(diag, w2);	}

	// Beginning of the inner loop.

L200:

	// Determine the levenberg-marquardt parameter.
	lmpar(&n, &fjac[fjac_offset], &ldfjac, ipvt.memptr(), diag.memptr(), qtf.memptr(),
		&delta, &par, w1.memptr(), w2.memptr(), w3.memptr(), w4.memptr());

	// Store the direction p and x + p. calculate the norm of p.
	w1 = -w1;
	w2 = x + w1;
	w3 = diag % w1;
	pnorm = enorm(w3);

	// On the first iteration, adjust the initial step bound. 
	if (iter == 1) {delta = min(delta, pnorm);}

	// Evaluate the function at x + p and calculate its norm.
	iflag = 1;
	(*fcn)(w2, w4, fjac1,iflag); ++(nfev);

	if (iflag < 0) {goto L300;}

	fnorm1 = enorm(w4);

	// compute the scaled actual reduction. */
	actred = -1.0;
	if (0.1 * fnorm1 < fnorm) {
		actred = 1.0 - pow(fnorm1 / fnorm,2);
	}

	// Compute the scaled predicted reduction and  the scaled directional derivative.
	w3.zeros();
	for (int jd = 1; jd <= n; ++jd) {
		for (int id = 1; id <= jd; ++id) {
			w3(id-1) += fjac[id + jd * fjac_dim1] * w1(ipvt(jd - 1) - 1);
		}
	}

	temp1 = pow(enorm(w3) / fnorm,2);
	temp2 = pow(sqrt(par) * pnorm / fnorm,2);
	prered = temp1 + temp2 / 0.5;
	dirder = -(temp1 + temp2);

	// Compute the ratio of the actual to the predicted reduction.

	ratio = 0.0;
	if (prered != 0.0) {ratio = actred / prered;}

	// Update the step bound.

	if (ratio <= 0.25) {
		if (actred >= 0.0) { temp = 0.5; }

		if (actred < 0.0) { temp = 0.5* dirder / (dirder + 0.5* actred); }

		if (0.1 * fnorm1 >= fnorm || temp < 0.1) { temp = 0.1; }

		delta = temp * min(delta, pnorm / 0.1);
		par /= temp;

	}
	else {
		if (par != 0.0 && ratio < 0.75) {
			// DO NOTHING
		}
		else {
			delta = pnorm / 0.5;
			par = 0.5* par;
		}
	}


	// Test for successful iteration.

	if (ratio >= 0.0001) {//  successful iteration. update x, fvec, and their norms. */
		x = w2;
		w2 = diag % x;
		fvec = w4;
		xnorm = enorm(w2);
		fnorm = fnorm1;
		++iter;
	}

	// tests for convergence.
	if (abs(actred) <= ftol && prered <= ftol && 0.5* ratio <= 1.0) { info = 1; }

	if (delta <= xtol * xnorm) { info = 2; }

	if (abs(actred) <= ftol && prered <= ftol && 0.5* ratio <= 1.0 && info == 2) { info = 3; }

	if (info != 0) { goto L300; }

	// tests for termination and stringent tolerances. */

	if (nfev >= maxfev) { info = 5; }

	if (abs(actred) <= epsmch && prered <= epsmch && 0.5* ratio <= 1.0) { info = 6; }

	if (delta <= epsmch * xnorm) { info = 7; }

	if (gnorm <= epsmch) { info = 8; }

	if (info != 0) { goto L300; }

	// End of the inner loop. repeat if iteration unsuccessful. */

	if (ratio < 0.0001) { goto L200; }

	// End of the outer loop.

	goto L30;
L300:

	/*     termination, either normal or user imposed. */

	if (iflag < 0) {
		info = iflag;
	}
	iflag = 0;
	if (nprint > 0) {
		(*fcn)(x, fvec, fjac1,iflag);
	}
	return 0;
}


int lmdif(decl_fcn_mn, int &m, int &n, Col<double> &x, Col<double> &fvec, double &ftol, double &xtol, double &gtol, int &maxfev, double &epsfcn, Col<double> &diag,
	int &mode, double &factor, int &nprint, int &info, int &nfev, double *fjac, int &ldfjac, int *ipvt, Col<double> &qtf)
{
	bool ft = true;
	int fjac_dim1, fjac_offset;
	double d__1, d__2, d__3;
	static int i__, j, l;
	static double par, sum;
	static int iter;
	static double temp, temp1, temp2;
	static int iflag;
	static double delta;
	static double ratio;
	static double fnorm, gnorm;
	static double pnorm;
	static double xnorm;
	static double fnorm1;
	static double actred;
	static double dirder;
	static double epsmch;
	static double prered;

	Col<double> w1(n), w2(n), w3(n), w4(m);  // Working vectors

	fjac_dim1 = ldfjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;

	// epsmch is the machine precision.
	epsmch = dpmpar(1);

	info = 0;
	iflag = 0;
	nfev = 0;

	// check the input parameters for errors.

	if (n <= 0 || m < n || ldfjac < m || ftol < 0.0 || xtol < 0.0 ||	gtol < 0.0 || maxfev <= 0 || factor <= 0.0) {
		goto L300;
	}

	if (mode == 1 && any(diag <= 0.0)) {
		goto L300;
	}

	// Evaluate the function at the starting point and calculate its norm.
	iflag = 1;
	(*fcn)(x, fvec, iflag); ++nfev;
	if (iflag < 0) {goto L300;}
	fnorm = enorm(fvec);

	// Initialize levenberg-marquardt parameter and iteration counter.
	par = 0.0;
	iter = 1;

	//Beginning of the outer loop.
L30:

	// Calculate the jacobian matrix.

	iflag = 2;
	fdjac2((S_fp)fcn, m, n, x, fvec, &fjac[fjac_offset], &ldfjac, &iflag, &epsfcn, w4.memptr()); nfev += n;
	if (iflag < 0) {goto L300;}


	// if requested, call fcn to enable printing of iterates.
	if (nprint > 0) {
		iflag = 0;
		if ((iter - 1) % nprint == 0) {
			(*fcn)(x, fvec, iflag);
		}
		if (iflag < 0) {goto L300;}
	}
	


	// compute the qr factorization of the jacobian.
	qrfac(&m, &n, &fjac[fjac_offset], &ldfjac, &ft, ipvt, &n, w1.memptr(), w2.memptr(), w3.memptr());

	// on the first iteration and if mode is 1, scale according
	// to the norms of the columns of the initial jacobian.

	if (iter == 1) {
		if (mode == 1) {
			diag = w2;
			diag.elem(find(w2 == 0)).ones();
		}

		// On the first iteration, calculate the norm of the scaled x and initialize the step bound delta. */
		w3 = diag % x;
		xnorm = enorm(w3);
		delta = factor * xnorm;
		if (delta == 0.0) { delta = factor; }
	}

	// form (q transpose)*fvec and store the first n com0.1.0nts in 
	// qtf.


	for (i__ = 1; i__ <= m; ++i__) {
		w4(i__-1) = fvec(i__-1);
	}
	for (j = 1; j <= n; ++j) {
		if (fjac[j + j * fjac_dim1] == 0.0) {
			goto L120;
		}
		sum = 0.0;
		for (i__ = j; i__ <= m; ++i__) {
			sum += fjac[i__ + j * fjac_dim1] * w4(i__-1);
		}
		temp = -sum / fjac[j + j * fjac_dim1];

		for (i__ = j; i__ <= m; ++i__) {
			w4(i__-1) += fjac[i__ + j * fjac_dim1] * temp;
		}
	L120:
		fjac[j + j * fjac_dim1] = w1(j-1);
		qtf(j-1) = w4(j-1);
	}

	// compute the norm of the scaled gradient.

	gnorm = 0.0;
	if (fnorm == 0.0) {
		goto L170;
	}

	for (j = 1; j <= n; ++j) {
		l = ipvt[j-1];
		if (w2(l-1) == 0.0) {
			goto L150;
		}
		sum = 0.0;
		for (i__ = 1; i__ <= j; ++i__) {
			sum += fjac[i__ + j * fjac_dim1] * (qtf(i__-1) / fnorm);
		}
		/* Computing MAX */
		d__2 = gnorm, d__3 = (d__1 = sum / w2(l-1), abs(d__1));
		gnorm = max(d__2, d__3);
	L150:
		;
	}
L170:

	/*        test for convergence of the gradient norm. */

	if (gnorm <= gtol) {info = 4;}
	if (info != 0) {goto L300;}

	// rescale if necessary. */
	if (mode == 1) { diag = arma::max(diag, w2); }

	/*        beginning of the inner loop. */

L200:

	// determine the levenberg-marquardt parameter. */

	lmpar(&n, &fjac[fjac_offset], &ldfjac, ipvt, diag.memptr(), qtf.memptr(), &delta,
		&par, w1.memptr(), w2.memptr(), w3.memptr(), w4.memptr());

	// store the direction p and x + p. calculate the norm of p.
	w1 = -w1;
	w2 = x + w1;
	w3 = diag % w1;
	pnorm = enorm(w3);
	
	// on the first iteration, adjust the initial step bound.
	if (iter == 1) {delta = min(delta, pnorm);}

	//  evaluate the function at x + p and calculate its norm.

	iflag = 1;
	(*fcn)(w2, w4, iflag); ++(nfev);

	if (iflag < 0) {goto L300;}

	fnorm1 = enorm(w4);

	// compute the scaled actual reduction.

	actred = -1.0;
	if (0.1 * fnorm1 < fnorm) {
		/* Computing 2nd power */
		d__1 = fnorm1 / fnorm;
		actred = 1.0 - d__1 * d__1;
	}

	/*           compute the scaled predicted reduction and */
	/*           the scaled directional derivative. */

	for (j = 1; j <= n; ++j) {
		w3(j-1) = 0.0;
		l = ipvt[j-1];
		temp = w1(l-1);
		for (i__ = 1; i__ <= j; ++i__) {
			w3(i__-1) += fjac[i__ + j * fjac_dim1] * temp;
		}
	}
	temp1 = enorm(w3) / fnorm;
	temp2 = sqrt(par) * pnorm / fnorm;
	/* Computing 2nd power */
	d__1 = temp1;
	/* Computing 2nd power */
	d__2 = temp2;
	prered = d__1 * d__1 + d__2 * d__2 / 0.5;
	/* Computing 2nd power */
	d__1 = temp1;
	/* Computing 2nd power */
	d__2 = temp2;
	dirder = -(d__1 * d__1 + d__2 * d__2);

	// compute the ratio of the actual to the predicted */
	// reduction.

	ratio = 0.0;
	if (prered != 0.0) {
		ratio = actred / prered;
	}

	// update the step bound.
	if (ratio <= 0.25) {
		if (actred >= 0.0) { temp = 0.5; }

		if (actred < 0.0) { temp = 0.5* dirder / (dirder + 0.5* actred); }

		if (0.1 * fnorm1 >= fnorm || temp < 0.1) { temp = 0.1; }

		delta = temp * min(delta, pnorm / 0.1);
		par /= temp;

	}
	else {
		if (par != 0.0 && ratio < 0.75) {
			// DO NOTHING
		}
		else {
			delta = pnorm / 0.5;
			par = 0.5* par;
		}
	}

	// Test for successful iteration.
	if (ratio >= 0.0001) {//  successful iteration. update x, fvec, and their norms. */
		x = w2;
		w2 = diag % x;
		fvec = w4;
		xnorm = enorm(w2);
		fnorm = fnorm1;
		++iter;
	}


	// tests for convergence. */

	if (abs(actred) <= ftol && prered <= ftol && 0.5* ratio <= 1.0) {info = 1;}

	if (delta <= xtol * xnorm) {info = 2;}

	if (abs(actred) <= ftol && prered <= ftol && 0.5* ratio <= 1.0 && info == 2) { info = 3;}

	if (info != 0) {goto L300;}

	// tests for termination and stringent tolerances. */

	if (nfev >= maxfev) {info = 5;}

	if (abs(actred) <= epsmch && prered <= epsmch && 0.5* ratio <= 1.0) {info = 6;}

	if (delta <= epsmch * xnorm) {info = 7;}

	if (gnorm <= epsmch) {info = 8;}

	if (info != 0) {goto L300;}

	// end of the inner loop. repeat if iteration unsuccessful. */

	if (ratio < 0.0001) {goto L200;}

	// end of the outer loop.

	goto L30;
L300:

	// termination, either normal or user imposed. */

	if (iflag < 0) {info = iflag;}

	iflag = 0;
	if (nprint > 0) {
		(*fcn)(x, fvec, iflag);
	}
	return 0;

}


int lmstr(funcderstr_mn fcn, int &m, int &n, Col<double> &x, Col<double> &fvec, double *fjac, int &ldfjac, double &ftol, double &xtol, double &gtol, int &maxfev, Col<double> diag,
	int &mode, double &factor, int &nprint, int &info, int &nfev, int &njev, Col<int> &ipvt, Col<double> &qtf)
{
	bool ft = true;
	int fjac_dim1, fjac_offset;
	double d__1, d__2;
	static int i__, j, l;
	static double par, sum;
	static bool sing;
	static int iter;
	static double temp, temp1, temp2;
	static int iflag;
	static double delta;
	static double ratio;
	static double fnorm;
	static double gnorm;
	static double pnorm;
	static double xnorm; 
	static double fnorm1; 
	static double actred; 
	static double dirder; 
	static double epsmch; 
	static double prered;

	Col<double> w1(n), w2(n), w3(n), w4(m);  // Working vectors

	fjac_dim1 = ldfjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;

	// epsmch is the machine precision.
	epsmch = dpmpar(1);

	info = 0;
	iflag = 0;
	nfev = 0;
	njev = 0;

	// check the input parameters for errors.

	if (n <= 0 || m < n || ldfjac < n || ftol < 0.0 || xtol < 0.0 ||
		gtol < 0.0 || maxfev <= 0 || factor <= 0.0) {
		goto L340;
	}

	if (mode == 1 && any(diag <= 0.0)) {
		goto L300;
	}

	// Evaluate the function at the starting point and calculate its norm.

	iflag = 1;
	(*fcn)(x, fvec, w3, iflag);++nfev;
	if (iflag < 0) {goto L340;}
	fnorm = enorm(fvec);

	// initialize levenberg-marquardt parameter and iteration counter.
	par = 0.0;
	iter = 1;

	// beginning of the outer loop.
L30:

	// if requested, call fcn to enable printing of iterates.

	if (nprint > 0) {
		iflag = 0;
		if ((iter - 1) % nprint == 0) {
			(*fcn)(x, fvec, w3, iflag);
		}
		if (iflag < 0) {goto L340;
		}
	}

	// compute the qr factorization of the jacobian matrix 
	// calculated 1.0 row at a time, while simultaneously
	// forming (q transpose)*fvec and storing the first 
	// n components nts in qtf. */
	qtf.zeros();
	for (j = 1; j <= n; ++j) {
		for (i__ = 1; i__ <= n; ++i__) {
			fjac[i__ + j * fjac_dim1] = 0.0;
		}
	}

	iflag = 2;
	for (uword id = 0; id < m; ++id) {
		(*fcn)(x, fvec, w3, iflag);
		if (iflag < 0) {goto L340;}

		temp = fvec(id);
		rwupdt(&n, &fjac[fjac_offset], &ldfjac, w3.memptr(), qtf.memptr(), &temp, w1.memptr(), w2.memptr());
		++iflag;
	}
	++(njev);

	// if the jacobian is rank deficient, call qrfac to
	// reorder its columns and update the com0.1.0nts of qtf.

	sing = false;

	for (j = 1; j <= n; ++j) {
		if (fjac[j + j * fjac_dim1] == 0.0) {
			sing = true;
		}
		ipvt(j-1) = j;
		w2(j-1) = enorm(&j, &fjac[j * fjac_dim1 + 1]);
	}

	if (!sing) {goto L130;}

	qrfac(&n, &n, &fjac[fjac_offset], &ldfjac, &ft, ipvt.memptr(), &n, w1.memptr(), w2.memptr(), w3.memptr());

	for (j = 1; j <= n; ++j) {
		if (fjac[j + j * fjac_dim1] == 0.0) {
			goto L110;
		}
		sum = 0.0;
		for (i__ = j; i__ <= n; ++i__) {
			sum += fjac[i__ + j * fjac_dim1] * qtf(i__-1);
		}
		temp = -sum / fjac[j + j * fjac_dim1];

		for (i__ = j; i__ <= n; ++i__) {
			qtf(i__-1) += fjac[i__ + j * fjac_dim1] * temp;
		}
	L110:
		fjac[j + j * fjac_dim1] = w1(j-1);
	}
L130:

	/*        on the first iteration and if mode is 1, scale according */
	/*        to the norms of the columns of the initial jacobian. */

	if (iter == 1) {
		if (mode == 1) {
			diag = w2;
			diag.elem(find(w2 == 0)).ones();
		}

		// on the first iteration, calculate the norm of the scaled x */
		// and initialize the step bound delta. */
		w3 = diag % x;
		xnorm = enorm(w3);
		delta = factor * xnorm;
		if (delta == 0.0) { delta = factor; }

	}

	// compute the norm of the scaled gradient. */

	gnorm = 0.0;
	if (fnorm != 0.0) {

		for (j = 1; j <= n; ++j) {
			l = ipvt(j - 1);
			if (w2(l - 1) != 0.0) {
				sum = 0.0;
				for (i__ = 1; i__ <= j; ++i__) {
					sum += fjac[i__ + j * fjac_dim1] * (qtf(i__ - 1) / fnorm);
				}
				gnorm = max(gnorm, abs(sum / w2(l - 1)));
			}
		}
	}

	// test for convergence of the gradient norm. */
	if (gnorm <= gtol) {info = 4;}

	if (info != 0) {goto L340;}


	//  rescale if necessary.
	if (mode == 1) { diag = arma::max(diag, w2);}

	// beginning of the inner loop

L240:

	// determine the levenberg-marquardt parameter. */

	lmpar(&n, &fjac[fjac_offset], &ldfjac, ipvt.memptr(), diag.memptr(), qtf.memptr(), &delta,
		&par, w1.memptr(), w2.memptr(), w3.memptr(), w4.memptr());

	// store the direction p and x + p. calculate the norm of p. */
	w1 = -w1;
	w2 = x + w1;
	w3 = diag % w1;
	pnorm = enorm(w3);

	// on the first iteration, adjust the initial step bound.
	if (iter == 1) {delta = min(delta, pnorm);}

	// evaluate the function at x + p and calculate its norm.
	iflag = 1;
	(*fcn)(w2, w4, w3, iflag);
	++(nfev);
	if (iflag < 0) {
		goto L340;
	}
	fnorm1 = enorm(w4);

	// compute the scaled actual reduction. */

	actred = -1.0;
	if (0.1 * fnorm1 < fnorm) {
		/* Computing 2nd power */
		d__1 = fnorm1 / fnorm;
		actred = 1.0 - d__1 * d__1;
	}

	// compute the scaled predicted reduction and */
	// the scaled directional derivative. */

	for (int jd = 1; jd <= n; ++jd) {
		w3(jd-1) = 0.0;
		l = ipvt[jd-1];
		temp = w1(l-1);
		for (int id = 1; id <= jd; ++id) {
			w3(id-1) += fjac[id + jd * fjac_dim1] * temp;
		}
	}

	temp1 = enorm(w3) / fnorm;
	temp2 = sqrt(par) * pnorm / fnorm;
	/* Computing 2nd power */
	d__1 = temp1;
	/* Computing 2nd power */
	d__2 = temp2;
	prered = d__1 * d__1 + d__2 * d__2 / 0.5;
	/* Computing 2nd power */
	d__1 = temp1;
	/* Computing 2nd power */
	d__2 = temp2;
	dirder = -(d__1 * d__1 + d__2 * d__2);

	// compute the ratio of the actual to the predicted reduction.

	ratio = 0.0;
	if (prered != 0.0) {
		ratio = actred / prered;
	}

	// update the step bound. 

	if (ratio <= 0.25) {
		if (actred >= 0.0) { temp = 0.5; }

		if (actred < 0.0) { temp = 0.5* dirder / (dirder + 0.5* actred); }

		if (0.1 * fnorm1 >= fnorm || temp < 0.1) { temp = 0.1; }

		delta = temp * min(delta, pnorm / 0.1);
		par /= temp;

	}
	else {
		if (par != 0.0 && ratio < 0.75) {
			// DO NOTHING
		}
		else {
			delta = pnorm / 0.5;
			par = 0.5* par;
		}
	}


L300:
	// Test for successful iteration.

	if (ratio >= 0.0001) {//  successful iteration. update x, fvec, and their norms. */
		x = w2;
		w2 = diag % x;
		fvec = w4;
		xnorm = enorm(w2);
		fnorm = fnorm1;
		++iter;
	}

	// tests for convergence.

	if (abs(actred) <= ftol && prered <= ftol && 0.5* ratio <= 1.0) {info = 1;}

	if (delta <= xtol * xnorm) {info = 2;}

	if (abs(actred) <= ftol && prered <= ftol && 0.5* ratio <= 1.0 && info == 2) {info = 3;}

	if (info != 0) {goto L340;}

	// tests for termination and stringent tolerances.
	if (nfev >= maxfev) {info = 5;}

	if (abs(actred) <= epsmch && prered <= epsmch && 0.5* ratio <= 1.0) {info = 6;}

	if (delta <= epsmch * xnorm) {info = 7;}

	if (gnorm <= epsmch) {info = 8;}

	if (info != 0) {goto L340;}

	// End of the inner loop. repeat if iteration unsuccessful.
	if (ratio < 0.0001) {goto L240;}

	// End of the outer loop.
	goto L30;
L340:

	// termination, either normal or user imposed.

	if (iflag < 0) { info = iflag; }
	iflag = 0;
	if (nprint > 0) { (*fcn)(x, fvec, w3, iflag); }

	return 0;
}


// Supporting functions

int fdjac2(S_fp fcn, int &m1, int &n1, Col<double> &x1, Col<double> &fvec1, double *fjac, int *ldfjac, int *iflag, double *epsfcn, double *wa)
{

	int *m = &m1;
	int *n = &n1;
	double *x = x1.memptr();
	double *fvec = fvec1.memptr();

	int fjac_dim1, fjac_offset, i__1, i__2;

	static double h__;
	static int i__, j;
	static double eps, temp, epsmch;

	Col<double> w1(m1);
	wa = w1.memptr();
	--wa;
	--fvec;
	--x;
	fjac_dim1 = *ldfjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;

	// epsmch is the machine precision.

	epsmch = dpmpar(1);

	eps = sqrt((max(*epsfcn, epsmch)));
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		temp = x[j];
		h__ = eps * abs(temp);
		if (h__ == 0.0) {
			h__ = eps;
		}
		x[j] = temp + h__;
		//int m1 = *m;
		//int n1 = *n;

		(*fcn)(&x1, &w1, iflag);
		if (*iflag < 0) {
			return 0;
		}
		x[j] = temp;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
			fjac[i__ + j * fjac_dim1] = (wa[i__] - fvec[i__]) / h__;
		}
	}
	return 0;
}

int fdjac1(S_fp fcn, int &n1, Col<double> &x1, Col<double> &fvec1, double *fjac, int *ldfjac, int *iflag, int *ml, int *mu, double *epsfcn, double *wa1, double *wa2)
{
	int *n = &n1;
	double *x = x1.memptr();
	double *fvec = fvec1.memptr();

	int fjac_dim1, fjac_offset, i__1, i__2, i__3, i__4;
	double d__1;
	static double h__;
	static int i__, j, k;
	static double eps, temp;
	static int msum;
	static double epsmch;
	Col<double> w1(n1);
	wa1 = w1.memptr();

	Col<double> w2(n1);
	wa2 = w2.memptr();

	--wa2;
	--wa1;
	--fvec;
	--x;
	fjac_dim1 = *ldfjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;

	//     epsmch is the machine precision.

	epsmch = dpmpar(1);

	eps = sqrt((max(*epsfcn, epsmch)));
	msum = *ml + *mu + 1;
	if (msum < *n) {
		goto L40;
	}

	//  computation of dense approximate jacobian. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		temp = x[j];
		h__ = eps * abs(temp);
		if (h__ == 0.0) {
			h__ = eps;
		}
		x[j] = temp + h__;
		//int n1= *n;
		(*fcn)(&x1, &w1, iflag);
		if (*iflag < 0) {
			return 0;
		}
		x[j] = temp;
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			fjac[i__ + j * fjac_dim1] = (wa1[i__] - fvec[i__]) / h__;
		}
	}
	return 0;
L40:

	//computation of banded approximate jacobian. 

	i__1 = msum;
	for (k = 1; k <= i__1; ++k) {
		i__2 = *n;
		i__3 = msum;
		for (j = k; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3) {
			wa2[j] = x[j];
			h__ = eps * (d__1 = wa2[j], abs(d__1));
			if (h__ == 0.0) {
				h__ = eps;
			}
			x[j] = wa2[j] + h__;
			/* L60: */
		}
		//int n1 = *n;
		(*fcn)(&x1, &w1, iflag);
		if (*iflag < 0) {
			return 0;
		}
		i__3 = *n;
		i__2 = msum;
		for (j = k; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) {
			x[j] = wa2[j];
			h__ = eps * (d__1 = wa2[j], abs(d__1));
			if (h__ == 0.0) {
				h__ = eps;
			}
			i__4 = *n;
			for (i__ = 1; i__ <= i__4; ++i__) {
				fjac[i__ + j * fjac_dim1] = 0.0;
				if (i__ >= j - *mu && i__ <= j + *ml) {
					fjac[i__ + j * fjac_dim1] = (wa1[i__] - fvec[i__]) / h__;
				}
			}
		}
	}
	return 0;
}



int lmpar(int *n, double *r__, int *ldr, int *ipvt, double *diag, double *qtb, double *delta, double *par, double *x, double *sdiag, double *wa1,
	double *wa2)
{
	int r_dim1, r_offset, i__1, i__2;
	double d__1, d__2;

	static int i__, j, k, l;
	static double fp;
	static int jm1, jp1;
	static double sum, parc, parl;
	static int iter;
	static double temp, paru, dwarf;
	static int nsing;
	static double gnorm;
	static double dxnorm;

	--wa2;
	--wa1;
	--sdiag;
	--x;
	--qtb;
	--diag;
	--ipvt;
	r_dim1 = *ldr;
	r_offset = 1 + r_dim1;
	r__ -= r_offset;


	// dwarf is the smallest positive magnitude. */
	dwarf = dpmpar(2);

	// compute and store in x the gauss-newton direction. if the */
	// jacobian is rank-deficient, obtain a least squares solution. */

	nsing = *n;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		wa1[j] = qtb[j];
		if (r__[j + j * r_dim1] == 0.0 && nsing == *n) {
			nsing = j - 1;
		}
		if (nsing < *n) {
			wa1[j] = 0.0;
		}
	}
	if (nsing < 1) {
		goto L50;
	}
	i__1 = nsing;
	for (k = 1; k <= i__1; ++k) {
		j = nsing - k + 1;
		wa1[j] /= r__[j + j * r_dim1];
		temp = wa1[j];
		jm1 = j - 1;
		if (jm1 < 1) {
			goto L30;
		}
		i__2 = jm1;
		for (i__ = 1; i__ <= i__2; ++i__) {
			wa1[i__] -= r__[i__ + j * r_dim1] * temp;
		}
	L30:
		;
	}
L50:
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		l = ipvt[j];
		x[l] = wa1[j];
	}

	/*     initialize the iteration counter. */
	/*     evaluate the function at the origin, and test */
	/*     for acceptance of the gauss-newton direction. */

	iter = 0;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		wa2[j] = diag[j] * x[j];
	}
	dxnorm = enorm(n, &wa2[1]);
	fp = dxnorm - *delta;
	if (fp <= 0.1 * *delta) {
		goto L220;
	}

	/*     if the jacobian is not rank deficient, the newton */
	/*     step provides a lower bound, parl, for the 0.0 of */
	/*     the function. otherwise set this bound to 0.0. */

	parl = 0.0;
	if (nsing < *n) {
		goto L120;
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		l = ipvt[j];
		wa1[j] = diag[l] * (wa2[l] / dxnorm);
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		sum = 0.0;
		jm1 = j - 1;
		if (jm1 < 1) {
			goto L100;
		}
		i__2 = jm1;
		for (i__ = 1; i__ <= i__2; ++i__) {
			sum += r__[i__ + j * r_dim1] * wa1[i__];
		}
	L100:
		wa1[j] = (wa1[j] - sum) / r__[j + j * r_dim1];
	}
	temp = enorm(n, &wa1[1]);
	parl = fp / *delta / temp / temp;
L120:

	/*     calculate an upper bound, paru, for the 0.0 of the function. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		sum = 0.0;
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
			sum += r__[i__ + j * r_dim1] * qtb[i__];
		}
		l = ipvt[j];
		wa1[j] = sum / diag[l];
	}
	gnorm = enorm(n, &wa1[1]);
	paru = gnorm / *delta;
	if (paru == 0.0) {
		paru = dwarf / min(*delta, 0.1);
	}

	/*     if the input par lies outside of the interval (parl,paru), */
	/*     set par to the closer endpoint. */

	*par = max(*par, parl);
	*par = min(*par, paru);
	if (*par == 0.0) {
		*par = gnorm / dxnorm;
	}

	/*     beginning of an iteration. */

L150:
	++iter;

	// evaluate the function at the current value of par.

	if (*par == 0.0) {
		/* Computing MAX */
		d__1 = dwarf, d__2 = 0.001 * paru;
		*par = max(d__1, d__2);
	}
	temp = sqrt(*par);
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		wa1[j] = temp * diag[j];
	}
	qrsolv(n, &r__[r_offset], ldr, &ipvt[1], &wa1[1], &qtb[1], &x[1], &sdiag[1], &wa2[1]);
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		wa2[j] = diag[j] * x[j];
	}
	dxnorm = enorm(n, &wa2[1]);
	temp = fp;
	fp = dxnorm - *delta;

	/*        if the function is small enough, accept the current value */
	/*        of par. also test for the exceptional cases where parl */
	/*        is 0.0 or the number of iterations has reached 10. */

	if (abs(fp) <= 0.1 * *delta || parl == 0.0 && fp <= temp && temp < 0.0 ||
		iter == 10) {
		goto L220;
	}

	/*        compute the newton correction. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		l = ipvt[j];
		wa1[j] = diag[l] * (wa2[l] / dxnorm);
		/* L180: */
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		wa1[j] /= sdiag[j];
		temp = wa1[j];
		jp1 = j + 1;
		if (*n < jp1) {
			goto L200;
		}
		i__2 = *n;
		for (i__ = jp1; i__ <= i__2; ++i__) {
			wa1[i__] -= r__[i__ + j * r_dim1] * temp;
		}
	L200:
		;
	}
	temp = enorm(n, &wa1[1]);
	parc = fp / *delta / temp / temp;

	/*        depending on the sign of the function, update parl or paru. */

	if (fp > 0.0) {
		parl = max(parl, *par);
	}
	if (fp < 0.0) {
		paru = min(paru, *par);
	}

	/*        compute an improved estimate for par. */

	/* Computing MAX */
	d__1 = parl, d__2 = *par + parc;
	*par = max(d__1, d__2);

	/*        end of an iteration. */

	goto L150;
L220:

	if (iter == 0) {
		*par = 0.0;
	}
	return 0;
}


int qform(int *m, int *n, double *q, int *	ldq, double *wa)
{
	int q_dim1, q_offset, i__1, i__2, i__3;
	static int i__, j, k, l, jm1, np1;
	static double sum, temp;
	static int minmn;

	--wa;
	q_dim1 = *ldq;
	q_offset = 1 + q_dim1;
	q -= q_offset;

	// 0.0 out upper triangle of q in the first min(m,n) columns. */

	minmn = min(*m, *n);
	if (minmn < 2) {
		goto L30;
	}
	i__1 = minmn;
	for (j = 2; j <= i__1; ++j) {
		jm1 = j - 1;
		i__2 = jm1;
		for (i__ = 1; i__ <= i__2; ++i__) {
			q[i__ + j * q_dim1] = 0.0;
		}
	}
L30:
	// initialize remaining columns to those of the identity matrix.

	np1 = *n + 1;
	if (*m < np1) {
		goto L60;
	}
	i__1 = *m;
	for (j = np1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
			q[i__ + j * q_dim1] = 0.0;
		}
		q[j + j * q_dim1] = 1.0;
	}
L60:

	//     accumulate q from its factored form.

	i__1 = minmn;
	for (l = 1; l <= i__1; ++l) {
		k = minmn - l + 1;
		i__2 = *m;
		for (i__ = k; i__ <= i__2; ++i__) {
			wa[i__] = q[i__ + k * q_dim1];
			q[i__ + k * q_dim1] = 0.0;
			/* L70: */
		}
		q[k + k * q_dim1] = 1.0;
		if (wa[k] == 0.0) {
			goto L110;
		}
		i__2 = *m;
		for (j = k; j <= i__2; ++j) {
			sum = 0.0;
			i__3 = *m;
			for (i__ = k; i__ <= i__3; ++i__) {
				sum += q[i__ + j * q_dim1] * wa[i__];
				/* L80: */
			}
			temp = sum / wa[k];
			i__3 = *m;
			for (i__ = k; i__ <= i__3; ++i__) {
				q[i__ + j * q_dim1] -= temp * wa[i__];
			}
		}
	L110:
		;
	}
	return 0;
}


int qrfac(int *m, int *n, double *a, int * lda, bool *pivot, int *ipvt, int *lipvt, double *rdiag, double *acnorm, double *wa)
{
	int a_dim1, a_offset, i__1, i__2, i__3;
	double d__1, d__2, d__3;
	static int i__, j, k, jp1;
	static double sum;
	static int kmax;
	static double temp;
	static int minmn;
	static double epsmch;
	static double ajnorm;

	--wa;
	--acnorm;
	--rdiag;
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--ipvt;


	// epsmch is the machine precision. 

	epsmch = dpmpar(1);

	// compute the initial column norms and initialize several arrays.

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		acnorm[j] = enorm(m, &a[j * a_dim1 + 1]);
		rdiag[j] = acnorm[j];
		wa[j] = rdiag[j];
		if (*pivot) {
			ipvt[j] = j;
		}
	}

	// reduce a to r with householder transformations.

	minmn = min(*m, *n);
	i__1 = minmn;
	for (j = 1; j <= i__1; ++j) {
		if (!(*pivot)) {
			goto L40;
		}

		// bring the column of largest norm into the pivot position.

		kmax = j;
		i__2 = *n;
		for (k = j; k <= i__2; ++k) {
			if (rdiag[k] > rdiag[kmax]) {
				kmax = k;
			}
			/* L20: */
		}
		if (kmax == j) {
			goto L40;
		}
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
			temp = a[i__ + j * a_dim1];
			a[i__ + j * a_dim1] = a[i__ + kmax * a_dim1];
			a[i__ + kmax * a_dim1] = temp;
		}
		rdiag[kmax] = rdiag[j];
		wa[kmax] = wa[j];
		k = ipvt[j];
		ipvt[j] = ipvt[kmax];
		ipvt[kmax] = k;
	L40:

		// compute the householder transformation to reduce the
		// j-th column of a to a multiple of the j-th unit vector.

		i__2 = *m - j + 1;
		ajnorm = enorm(&i__2, &a[j + j * a_dim1]);
		if (ajnorm == 0.0) {
			goto L100;
		}
		if (a[j + j * a_dim1] < 0.0) {
			ajnorm = -ajnorm;
		}
		i__2 = *m;
		for (i__ = j; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] /= ajnorm;
		}
		a[j + j * a_dim1] += 1.0;

		// apply the transformation to the remaining columns
		// and update the norms.

		jp1 = j + 1;
		if (*n < jp1) {
			goto L100;
		}
		i__2 = *n;
		for (k = jp1; k <= i__2; ++k) {
			sum = 0.0;
			i__3 = *m;
			for (i__ = j; i__ <= i__3; ++i__) {
				sum += a[i__ + j * a_dim1] * a[i__ + k * a_dim1];
				/* L60: */
			}
			temp = sum / a[j + j * a_dim1];
			i__3 = *m;
			for (i__ = j; i__ <= i__3; ++i__) {
				a[i__ + k * a_dim1] -= temp * a[i__ + j * a_dim1];
				/* L70: */
			}
			if (!(*pivot) || rdiag[k] == 0.0) {
				goto L80;
			}
			temp = a[j + k * a_dim1] / rdiag[k];
			// Computing MAX
			// Computing 2nd power
			d__3 = temp;
			d__1 = 0.0, d__2 = 1.0 - d__3 * d__3;
			rdiag[k] *= sqrt((max(d__1, d__2)));
			// Computing 2nd power
			d__1 = rdiag[k] / wa[k];
			if (0.05 * (d__1 * d__1) > epsmch) {
				goto L80;
			}
			i__3 = *m - j;
			rdiag[k] = enorm(&i__3, &a[jp1 + k * a_dim1]);
			wa[k] = rdiag[k];
		L80:
			;
		}
	L100:
		rdiag[j] = -ajnorm;
	}
	return 0;
}


int qrsolv(int *n, double *r__, int *ldr, int *ipvt, double *diag, double *qtb, double *x, double *sdiag, double *wa)
{
	int r_dim1, r_offset, i__1, i__2, i__3;
	double d__1, d__2;

	static int i__, j, k, l, jp1, kp1;
	static double tan__, cos__, sin__, sum, temp, cotan;
	static int nsing;
	static double qtbpj;


	--wa;
	--sdiag;
	--x;
	--qtb;
	--diag;
	--ipvt;
	r_dim1 = *ldr;
	r_offset = 1 + r_dim1;
	r__ -= r_offset;

	/*     copy r and (q transpose)*b to preserve input and initialize s. */
	/*     in particular, save the diagonal elements of r in x. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
			r__[i__ + j * r_dim1] = r__[j + i__ * r_dim1];
			/* L10: */
		}
		x[j] = r__[j + j * r_dim1];
		wa[j] = qtb[j];
	}

	/*     eliminate the diagonal matrix d using a givens rotation. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

		/*        prepare the row of d to be eliminated, locating the */
		/*        diagonal element using p from the qr factorization. */

		l = ipvt[j];
		if (diag[l] == 0.0) {
			goto L90;
		}
		i__2 = *n;
		for (k = j; k <= i__2; ++k) {
			sdiag[k] = 0.0;
		}
		sdiag[j] = diag[l];

		/*        the transformations to eliminate the row of d */
		/*        modify only a single element of (q transpose)*b */
		/*        beyond the first n, which is initially 0.0. */

		qtbpj = 0.0;
		i__2 = *n;
		for (k = j; k <= i__2; ++k) {

			/*           determine a givens rotation which eliminates the */
			/*           appropriate element in the current row of d. */

			if (sdiag[k] == 0.0) {
				goto L70;
			}
			if ((d__1 = r__[k + k * r_dim1], abs(d__1)) >= (d__2 = sdiag[k],
				abs(d__2))) {
				goto L40;
			}
			cotan = r__[k + k * r_dim1] / sdiag[k];
			/* Computing 2nd power */
			d__1 = cotan;
			sin__ = 0.5 / sqrt(0.25 + 0.25 * (d__1 * d__1));
			cos__ = sin__ * cotan;
			goto L50;
		L40:
			tan__ = sdiag[k] / r__[k + k * r_dim1];
			/* Computing 2nd power */
			d__1 = tan__;
			cos__ = 0.5 / sqrt(0.25 + 0.25 * (d__1 * d__1));
			sin__ = cos__ * tan__;
		L50:

			/*           compute the modified diagonal element of r and */
			/*           the modified element of ((q transpose)*b,0). */

			r__[k + k * r_dim1] = cos__ * r__[k + k * r_dim1] + sin__ * sdiag[
				k];
			temp = cos__ * wa[k] + sin__ * qtbpj;
			qtbpj = -sin__ * wa[k] + cos__ * qtbpj;
			wa[k] = temp;

			/*           accumulate the tranformation in the row of s. */

			kp1 = k + 1;
			if (*n < kp1) {
				goto L70;
			}
			i__3 = *n;
			for (i__ = kp1; i__ <= i__3; ++i__) {
				temp = cos__ * r__[i__ + k * r_dim1] + sin__ * sdiag[i__];
				sdiag[i__] = -sin__ * r__[i__ + k * r_dim1] + cos__ * sdiag[
					i__];
				r__[i__ + k * r_dim1] = temp;
			}
		L70:
			;
		}
	L90:

		/*        store the diagonal element of s and restore */
		/*        the corresponding diagonal element of r. */

		sdiag[j] = r__[j + j * r_dim1];
		r__[j + j * r_dim1] = x[j];
	}

	/*     solve the triangular system for z. if the system is */
	/*     singular, then obtain a least squares solution. */

	nsing = *n;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		if (sdiag[j] == 0.0 && nsing == *n) {
			nsing = j - 1;
		}
		if (nsing < *n) {
			wa[j] = 0.0;
		}
		/* L110: */
	}
	if (nsing < 1) {
		goto L150;
	}
	i__1 = nsing;
	for (k = 1; k <= i__1; ++k) {
		j = nsing - k + 1;
		sum = 0.0;
		jp1 = j + 1;
		if (nsing < jp1) {
			goto L130;
		}
		i__2 = nsing;
		for (i__ = jp1; i__ <= i__2; ++i__) {
			sum += r__[i__ + j * r_dim1] * wa[i__];
		}
	L130:
		wa[j] = (wa[j] - sum) / sdiag[j];
	}
L150:

	/*     permute the com0.1.0nts of z back to com0.1.0nts of x. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		l = ipvt[j];
		x[l] = wa[j];
		/* L160: */
	}
	return 0;
}


int r1mpyq(int *m, int *n, double *a, int * lda, double *v, double *w)
{
	int a_dim1, a_offset, i__1, i__2;
	double d__1, d__2;
	static int i__, j, nm1, nmj;
	static double cos__, sin__, temp;

	--w;
	--v;
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	// apply the first set of givens rotations to a.

	nm1 = *n - 1;
	if (nm1 < 1) {
		goto L50;
	}
	i__1 = nm1;
	for (nmj = 1; nmj <= i__1; ++nmj) {
		j = *n - nmj;
		if ((d__1 = v[j], abs(d__1)) > 1.0) {
			cos__ = 1.0 / v[j];
		}
		if ((d__1 = v[j], abs(d__1)) > 1.0) {
			/* Computing 2nd power */
			d__2 = cos__;
			sin__ = sqrt(1.0 - d__2 * d__2);
		}
		if ((d__1 = v[j], abs(d__1)) <= 1.0) {
			sin__ = v[j];
		}
		if ((d__1 = v[j], abs(d__1)) <= 1.0) {
			/* Computing 2nd power */
			d__2 = sin__;
			cos__ = sqrt(1.0 - d__2 * d__2);
		}
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
			temp = cos__ * a[i__ + j * a_dim1] - sin__ * a[i__ + *n * a_dim1];
			a[i__ + *n * a_dim1] = sin__ * a[i__ + j * a_dim1] + cos__ * a[
				i__ + *n * a_dim1];
			a[i__ + j * a_dim1] = temp;
		}
	}

	/*     apply the second set of givens rotations to a. */

	i__1 = nm1;
	for (j = 1; j <= i__1; ++j) {
		if ((d__1 = w[j], abs(d__1)) > 1.0) {
			cos__ = 1.0 / w[j];
		}
		if ((d__1 = w[j], abs(d__1)) > 1.0) {
			/* Computing 2nd power */
			d__2 = cos__;
			sin__ = sqrt(1.0 - d__2 * d__2);
		}
		if ((d__1 = w[j], abs(d__1)) <= 1.0) {
			sin__ = w[j];
		}
		if ((d__1 = w[j], abs(d__1)) <= 1.0) {
			/* Computing 2nd power */
			d__2 = sin__;
			cos__ = sqrt(1.0 - d__2 * d__2);
		}
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
			temp = cos__ * a[i__ + j * a_dim1] + sin__ * a[i__ + *n * a_dim1];
			a[i__ + *n * a_dim1] = -sin__ * a[i__ + j * a_dim1] + cos__ * a[
				i__ + *n * a_dim1];
			a[i__ + j * a_dim1] = temp;
		}
	}
L50:
	return 0;
}


int r1updt(int *m, int *n, double *s, int * ls, double *u, double *v, double *w, bool *sing)
{
	int i__1, i__2;
	double d__1, d__2;
	static int i__, j, l, jj, nm1;
	static double tan__;
	static int nmj;
	static double cos__, sin__, tau, temp, giant, cotan;
	--w;
	--u;
	--v;
	--s;
	// Giant is the largest magnitude. */

	giant = dpmpar(3);

	/*     initialize the diagonal element pointer. */

	jj = *n * ((*m << 1) - *n + 1) / 2 - (*m - *n);

	/*     move the nontrivial part of the last column of s into w. */

	l = jj;
	i__1 = *m;
	for (i__ = *n; i__ <= i__1; ++i__) {
		w[i__] = s[l];
		++l;
	}

	// rotate the vector v into a multiple of the n-th unit vector
	// in such a way that a spike is introduced into w.

	nm1 = *n - 1;
	if (nm1 < 1) {
		goto L70;
	}
	i__1 = nm1;
	for (nmj = 1; nmj <= i__1; ++nmj) {
		j = *n - nmj;
		jj -= *m - j + 1;
		w[j] = 0.0;
		if (v[j] == 0.0) {
			goto L50;
		}

		/*        determine a givens rotation which eliminates the */
		/*        j-th element of v. */

		if ((d__1 = v[*n], abs(d__1)) >= (d__2 = v[j], abs(d__2))) {
			goto L20;
		}
		cotan = v[*n] / v[j];
		/* Computing 2nd power */
		d__1 = cotan;
		sin__ = 0.5 / sqrt(0.25 + 0.25 * (d__1 * d__1));
		cos__ = sin__ * cotan;
		tau = 1.0;
		if (abs(cos__) * giant > 1.0) {
			tau = 1.0 / cos__;
		}
		goto L30;
	L20:
		tan__ = v[j] / v[*n];
		/* Computing 2nd power */
		d__1 = tan__;
		cos__ = 0.5 / sqrt(0.25 + 0.25 * (d__1 * d__1));
		sin__ = cos__ * tan__;
		tau = sin__;
	L30:

		// apply the transformation to v and store the information
		// necessary to recover the givens rotation.

		v[*n] = sin__ * v[j] + cos__ * v[*n];
		v[j] = tau;

		// apply the transformation to s and extend the spike in w.

		l = jj;
		i__2 = *m;
		for (i__ = j; i__ <= i__2; ++i__) {
			temp = cos__ * s[l] - sin__ * w[i__];
			w[i__] = sin__ * s[l] + cos__ * w[i__];
			s[l] = temp;
			++l;
		}
	L50:
		;
	}
L70:

	// add the spike from the rank 1 update to w.

	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
		w[i__] += v[*n] * u[i__];
	}

	/*     eliminate the spike. */

	*sing = false;
	if (nm1 < 1) {
		goto L140;
	}
	i__1 = nm1;
	for (j = 1; j <= i__1; ++j) {
		if (w[j] == 0.0) {
			goto L120;
		}

		// determine a givens rotation which eliminates the
		// j-th element of the spike.

		if ((d__1 = s[jj], abs(d__1)) >= (d__2 = w[j], abs(d__2))) {
			goto L90;
		}
		cotan = s[jj] / w[j];
		/* Computing 2nd power */
		d__1 = cotan;
		sin__ = 0.5 / sqrt(0.25 + 0.25 * (d__1 * d__1));
		cos__ = sin__ * cotan;
		tau = 1.0;
		if (abs(cos__) * giant > 1.0) {
			tau = 1.0 / cos__;
		}
		goto L100;
	L90:
		tan__ = w[j] / s[jj];
		/* Computing 2nd power */
		d__1 = tan__;
		cos__ = 0.5 / sqrt(0.25 + 0.25 * (d__1 * d__1));
		sin__ = cos__ * tan__;
		tau = sin__;
	L100:

		// apply the transformation to s and reduce the spike in w.

		l = jj;
		i__2 = *m;
		for (i__ = j; i__ <= i__2; ++i__) {
			temp = cos__ * s[l] + sin__ * w[i__];
			w[i__] = -sin__ * s[l] + cos__ * w[i__];
			s[l] = temp;
			++l;
		}

		// store the information necessary to recover the */
		// givens rotation.

		w[j] = tau;
	L120:

		// test for 0.0 diagonal elements in the output s.

		if (s[jj] == 0.0) {
			*sing = true;
		}
		jj += *m - j + 1;
	}
L140:

	// move w back into the last column of the output s.

	l = jj;
	i__1 = *m;
	for (i__ = *n; i__ <= i__1; ++i__) {
		s[l] = w[i__];
		++l;
	}
	if (s[jj] == 0.0) {
		*sing = true;
	}
	return 0;
}


int rwupdt(int *n, double *r__, int *ldr, double *w, double *b, double *alpha, double *cos__, double *sin__)
{
	int r_dim1, r_offset, i__1, i__2;
	double d__1;

	static int i__, j, jm1;
	static double tan__, temp, rowj, cotan;

	--sin__;
	--cos__;
	--b;
	--w;
	r_dim1 = *ldr;
	r_offset = 1 + r_dim1;
	r__ -= r_offset;

	/* Function Body */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		rowj = w[j];
		jm1 = j - 1;

		/*        apply the previous transformations to */
		/*        r(i,j), i=1,2,...,j-1, and to w(j). */

		if (jm1 < 1) {
			goto L20;
		}
		i__2 = jm1;
		for (i__ = 1; i__ <= i__2; ++i__) {
			temp = cos__[i__] * r__[i__ + j * r_dim1] + sin__[i__] * rowj;
			rowj = -sin__[i__] * r__[i__ + j * r_dim1] + cos__[i__] * rowj;
			r__[i__ + j * r_dim1] = temp;
			/* L10: */
		}
	L20:

		/*        determine a givens rotation which eliminates w(j). */

		cos__[j] = 1.0;
		sin__[j] = 0.0;
		if (rowj == 0.0) {
			goto L50;
		}
		if ((d__1 = r__[j + j * r_dim1], abs(d__1)) >= abs(rowj)) {
			goto L30;
		}
		cotan = r__[j + j * r_dim1] / rowj;
		/* Computing 2nd power */
		d__1 = cotan;
		sin__[j] = 0.5 / sqrt(0.25 + 0.25 * (d__1 * d__1));
		cos__[j] = sin__[j] * cotan;
		goto L40;
	L30:
		tan__ = rowj / r__[j + j * r_dim1];
		/* Computing 2nd power */
		d__1 = tan__;
		cos__[j] = 0.5 / sqrt(0.25 + 0.25 * (d__1 * d__1));
		sin__[j] = cos__[j] * tan__;
	L40:

		/*        apply the current transformation to r(j,j), b(j), and alpha. */

		r__[j + j * r_dim1] = cos__[j] * r__[j + j * r_dim1] + sin__[j] *
			rowj;
		temp = cos__[j] * b[j] + sin__[j] * *alpha;
		*alpha = -sin__[j] * b[j] + cos__[j] * *alpha;
		b[j] = temp;
	L50:
		;
	}
	return 0;
}


int chkder(int *m, int *n, double *x, double *fvec, double *fjac, int *ldfjac, double *xp, double *fvecp, int *mode, double *err)
{
	static double factor = 100.;

	int fjac_dim1, fjac_offset, i__1, i__2;
	double d__1, d__2, d__3, d__4, d__5;
	static int i__, j;
	static double eps, epsf, temp, epsmch;
	static double epslog;

	--err;
	--fvecp;
	--fvec;
	--xp;
	--x;
	fjac_dim1 = *ldfjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;

	/*     epsmch is the machine precision. */

	epsmch = dpmpar(1);

	eps = sqrt(epsmch);

	if (*mode == 2) {
		goto L20;
	}

	// mode = 1

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		temp = eps * (d__1 = x[j], abs(d__1));
		if (temp == 0.0) {
			temp = eps;
		}
		xp[j] = x[j] + temp;
	}
	return 0;
L20:

	/*        mode = 2. */

	epsf = factor * epsmch;
	epslog = log10(eps);
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
		err[i__] = 0.0;
		/* L30: */
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		temp = (d__1 = x[j], abs(d__1));
		if (temp == 0.0) {
			temp = 1.0;
		}
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
			err[i__] += temp * fjac[i__ + j * fjac_dim1];
			/* L40: */
		}
		/* L50: */
	}
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
		temp = 1.0;
		if (fvec[i__] != 0.0 && fvecp[i__] != 0.0 && (d__2 = fvecp[i__] -
			fvec[i__], abs(d__2)) >= epsf * (d__1 = fvec[i__], abs(d__1)))
		{
			temp = eps * (d__3 = (fvecp[i__] - fvec[i__]) / eps - err[i__],
				abs(d__3)) / ((d__4 = fvec[i__], abs(d__4)) + (d__5 =
					fvecp[i__], abs(d__5)));
		}
		err[i__] = 1.0;
		if (temp > epsmch && temp < eps) {
			err[i__] = (log10(temp) - epslog) / epslog;
		}
		if (temp >= eps) {
			err[i__] = 0.0;
		}
	}
	return 0;
}


int dogleg(int *n, double *r__, int *lr, double *diag, double *qtb, double *delta, double *x, double *wa1, double *wa2)
{
	int i__1, i__2;
	double d__1, d__2, d__3, d__4;
	static int i__, j, k, l, jj, jp1;
	static double sum, temp, alpha, bnorm;
	static double gnorm, qnorm, epsmch;
	static double sgnorm;


	--wa2;
	--wa1;
	--x;
	--qtb;
	--diag;
	--r__;


	// epsmch is the machine precision. */

	epsmch = dpmpar(1);

	// first, calculate the gauss-newton direction. */

	jj = *n * (*n + 1) / 2 + 1;
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
		j = *n - k + 1;
		jp1 = j + 1;
		jj -= k;
		l = jj + 1;
		sum = 0.0;
		if (*n < jp1) {
			goto L20;
		}
		i__2 = *n;
		for (i__ = jp1; i__ <= i__2; ++i__) {
			sum += r__[l] * x[i__];
			++l;
		}
	L20:
		temp = r__[jj];
		if (temp != 0.0) {
			goto L40;
		}
		l = j;
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* Computing MAX */
			d__2 = temp, d__3 = (d__1 = r__[l], abs(d__1));
			temp = max(d__2, d__3);
			l = l + *n - i__;
		}
		temp = epsmch * temp;
		if (temp == 0.0) {
			temp = epsmch;
		}
	L40:
		x[j] = (qtb[j] - sum) / temp;
	}

	/*     test whether the gauss-newton direction is acceptable. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		wa1[j] = 0.0;
		wa2[j] = diag[j] * x[j];
	}
	qnorm = enorm(n, &wa2[1]);
	if (qnorm <= *delta) {
		return 0;
	}

	/*     the gauss-newton direction is not acceptable. */
	/*     next, calculate the scaled gradient direction. */

	l = 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		temp = qtb[j];
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
			wa1[i__] += r__[l] * temp;
			++l;
		}
		wa1[j] /= diag[j];
	}

	/*     calculate the norm of the scaled gradient and test for */
	/*     the special case in which the scaled gradient is 0.0. */

	gnorm = enorm(n, &wa1[1]);
	sgnorm = 0.0;
	alpha = *delta / qnorm;
	if (gnorm == 0.0) {
		goto L120;
	}

	/*     calculate the point along the scaled gradient */
	/*     at which the quadratic is minimized. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		wa1[j] = wa1[j] / gnorm / diag[j];
		/* L90: */
	}
	l = 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		sum = 0.0;
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
			sum += r__[l] * wa1[i__];
			++l;
			/* L100: */
		}
		wa2[j] = sum;
		/* L110: */
	}
	temp = enorm(n, &wa2[1]);
	sgnorm = gnorm / temp / temp;

	/*     test whether the scaled gradient direction is acceptable. */

	alpha = 0.0;
	if (sgnorm >= *delta) {
		goto L120;
	}

	/*     the scaled gradient direction is not acceptable. */
	/*     finally, calculate the point along the dogleg */
	/*     at which the quadratic is minimized. */

	bnorm = enorm(n, &qtb[1]);
	temp = bnorm / gnorm * (bnorm / qnorm) * (sgnorm / *delta);
	/* Computing 2nd power */
	d__1 = sgnorm / *delta;
	/* Computing 2nd power */
	d__2 = temp - *delta / qnorm;
	/* Computing 2nd power */
	d__3 = *delta / qnorm;
	/* Computing 2nd power */
	d__4 = sgnorm / *delta;
	temp = temp - *delta / qnorm * (d__1 * d__1) + sqrt(d__2 * d__2 + (1.0 -
		d__3 * d__3) * (1.0 - d__4 * d__4));
	/* Computing 2nd power */
	d__1 = sgnorm / *delta;
	alpha = *delta / qnorm * (1.0 - d__1 * d__1) / temp;
L120:

	/*     form appropriate convex combination of the gauss-newton */
	/*     direction and the scaled gradient direction. */

	temp = (1.0 - alpha) * min(sgnorm, *delta);
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		x[j] = temp * wa1[j] + alpha * x[j];
	}
	return 0;
}

double dpmpar(int id)
{
	static struct {
		double e_1[3];
		double fill_2[1];
	} equiv_2 = { 2.22044604926e-16, 2.22507385852e-308,1.79769313485e308 };

	double ret_val;

	/* Local variables */
#define dmach ((double *)&equiv_2)
#define minmag ((int *)&equiv_2 + 2)
#define maxmag ((int *)&equiv_2 + 4)
#define mcheps ((int *)&equiv_2)




	ret_val = dmach[(0 + (0 + (id - 1 << 3))) / 8];
	return ret_val;

}

#undef mcheps
#undef maxmag
#undef minmag
#undef dmach

double enorm(int *n, double *x)
{
	/* given an n - vector x, this function calculates the euclidean norm of x.

			The euclidean norm is computed by accumulating the sum of squares in three different sums. The sums of squares for the
		    small and large components are scaled so that no overflows occur. Non - destructive underflows are permitted. Underflows
		    and overflows do not occur in the computation of the unscaled sum of squares for the intermediate components.
		    The definitions of small, intermediate and large components depend on two constants, rdwarf and rgiant. The main
		    restrictions on these constants are that rdwarf^2 not underflow and rgiant^2 not overflow. The constants
		    given here are suitable for every known computer.

	*/
	static double rdwarf = 3.834e-20;
	static double rgiant = 1.304e19;

	int i__1;
	double ret_val, d__1;
	static int i__;
	static double s1, s2, s3, xabs, x1max, x3max, agiant, floatn;

	--x;
	s1 = 0.0;
	s2 = 0.0;
	s3 = 0.0;
	x1max = 0.0;
	x3max = 0.0;
	floatn = (double)(*n);
	agiant = rgiant / floatn;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		xabs = (d__1 = x[i__], abs(d__1));
		if (xabs > rdwarf && xabs < agiant) {
			goto L70;
		}
		if (xabs <= rdwarf) {
			goto L30;
		}

		/*              sum for large com0.1.0nts. */

		if (xabs <= x1max) {
			goto L10;
		}
		/* Computing 2nd power */
		d__1 = x1max / xabs;
		s1 = 1.0 + s1 * (d__1 * d__1);
		x1max = xabs;
		goto L20;
	L10:
		/* Computing 2nd power */
		d__1 = xabs / x1max;
		s1 += d__1 * d__1;
	L20:
		goto L60;
	L30:

		/*              sum for small com0.1.0nts. */

		if (xabs <= x3max) {
			goto L40;
		}
		/* Computing 2nd power */
		d__1 = x3max / xabs;
		s3 = 1.0 + s3 * (d__1 * d__1);
		x3max = xabs;
		goto L50;
	L40:
		if (xabs != 0.0) {
			/* Computing 2nd power */
			d__1 = xabs / x3max;
			s3 += d__1 * d__1;
		}
	L50:
	L60:
		goto L80;
	L70:

		/*           sum for intermediate com0.1.0nts. */

		/* Computing 2nd power */
		d__1 = xabs;
		s2 += d__1 * d__1;
	L80:
		/* L90: */
		;
	}

	/*     calculation of norm. */

	if (s1 == 0.0) {
		goto L100;
	}
	ret_val = x1max * sqrt(s1 + s2 / x1max / x1max);
	goto L130;
L100:
	if (s2 == 0.0) {
		goto L110;
	}
	if (s2 >= x3max) {
		ret_val = sqrt(s2 * (1.0 + x3max / s2 * (x3max * s3)));
	}
	if (s2 < x3max) {
		ret_val = sqrt(x3max * (s2 / x3max + x3max * s3));
	}
	goto L120;
L110:
	ret_val = x3max * sqrt(s3);
L120:
L130:
	return ret_val;
}

double enorm(Col<double> &x1)
{
	/* given an n - vector x, this function calculates the euclidean norm of x.

			The euclidean norm is computed by accumulating the sum of squares in three different sums. The sums of squares for the
			small and large components are scaled so that no overflows occur. Non - destructive underflows are permitted. Underflows
			and overflows do not occur in the computation of the unscaled sum of squares for the intermediate components.
			The definitions of small, intermediate and large components depend on two constants, rdwarf and rgiant. The main
			restrictions on these constants are that rdwarf^2 not underflow and rgiant^2 not overflow. The constants
			given here are suitable for every known computer.

	*/
	static double rdwarf = 3.834e-20;
	static double rgiant = 1.304e19;

	uword n = x1.n_elem;
	double *x = x1.memptr();
	double ret_val, d__1;
	static int i__;
	static double s1, s2, s3, xabs, x1max, x3max, agiant, floatn;
	--x;
	s1 = 0.0;
	s2 = 0.0;
	s3 = 0.0;
	x1max = 0.0;
	x3max = 0.0;
	floatn = (double)(n);
	agiant = rgiant / floatn;

	for (uword i__ = 1; i__ <= n; ++i__) {
		xabs = (d__1 = x[i__], abs(d__1));
		if (xabs > rdwarf && xabs < agiant) {
			goto L70;
		}
		if (xabs <= rdwarf) {
			goto L30;
		}

		// sum for large com0.1.0nts.

		if (xabs <= x1max) {
			goto L10;
		}
		// Computing 2nd power
		d__1 = x1max / xabs;
		s1 = 1.0 + s1 * (d__1 * d__1);
		x1max = xabs;
		goto L20;
	L10:
		// Computing 2nd power
		d__1 = xabs / x1max;
		s1 += d__1 * d__1;
	L20:
		goto L50;
	L30:

		// sum for small com0.1.0nts.

		if (xabs <= x3max) {
			goto L40;
		}
		/* Computing 2nd power */
		d__1 = x3max / xabs;
		s3 = 1.0 + s3 * (d__1 * d__1);
		x3max = xabs;
		goto L50;
	L40:
		if (xabs != 0.0) {
			// Computing 2nd power
			d__1 = xabs / x3max;
			s3 += d__1 * d__1;
		}
	L50:
		goto L80;
	L70:

		// sum for intermediate com0.1.0nts.

		// Computing 2nd power
		d__1 = xabs;
		s2 += d__1 * d__1;
	L80:
		;
	}

	// calculation of norm. 

	if (s1 == 0.0) {goto L100;}

	ret_val = x1max * sqrt(s1 + s2 / x1max / x1max);
	goto L120;
L100:
	if (s2 == 0.0) {
		goto L110;
	}
	if (s2 >= x3max) {
		ret_val = sqrt(s2 * (1.0 + x3max / s2 * (x3max * s3)));
	}
	if (s2 < x3max) {
		ret_val = sqrt(x3max * (s2 / x3max + x3max * s3));
	}
	goto L120;
L110:
	ret_val = x3max * sqrt(s3);
L120:
	return ret_val;
}
