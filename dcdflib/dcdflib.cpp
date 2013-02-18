// This is the main DLL file.

#pragma warning (default : 4412)

#include "dcdflib.h"

using namespace dcdflib;
using namespace System;

StatClass::StatClass()
{
	tol = 1.0e-10;
	atol = 1.0e-50;
	zero = 1.0e-300;
	inf = 1.0e300;
	one = 1.0;
	tent4 = 1.0e4;
}

// NOTE
// Calls to dpmpar have been replaced with constants.
// dpmpar(1) Precision: 2.22044604925031E-16
// dpmpar(2) Smallest Magniture: 4.94065645841247e-324; (double::Epsilon)
// dpmpar(3) Largest Magnitude: 1.79769313486232E+308 (double::MaxValue)


//****************************************************************************80
//
//  Purpose:
// 
//    ALGDIV computes ln ( Gamma ( B ) / Gamma ( A + B ) ) when 8 <= B.
//
//  Discussion:
//
//    In this algorithm, DEL(X) is the function defined by
//
//      ln ( Gamma(X) ) = ( X - 0.5 ) * ln ( X ) - X + 0.5 * ln ( 2 * PI ) 
//                      + DEL(X).
//
//  Parameters:
//
//    Input, double %A, *B, define the arguments.
//
//    Output, double ALGDIV, the value of ln(Gamma(B)/Gamma(A+B)).
//

double StatClass::algdiv ( double %a, double %b )
{
	const double c0 =  0.833333333333333e-01;
	const double c1 = -0.277777777760991e-02;
	const double c2 =  0.793650666825390e-03;
	const double c3 = -0.595202931351870e-03;
	const double c4 =  0.837308034031215e-03;
	const double c5 = -0.165322962780713e-02;
	double h = 0.0;
	double c = 0.0;
	double d = 0.0;
	double x = 0.0;

	if ( b <= a )
	{
		h = b / a;
		c = 1.0 / ( 1.0 + h );
		x = h / ( 1.0 + h );
		d = a + ( b - 0.5 );
	}
	else
	{
		h = a / b;
		c = h / ( 1.0 + h );
		x = 1.0 / ( 1.0 + h );
		d = b + ( a - 0.5 );
	}

	//
	//  SET SN = (1 - X**N)/(1 - X)
	//
	double x2 = x * x;
	double s3 = 1.0 + ( x + x2 );
	double s5 = 1.0 + ( x + x2 * s3 );
	double s7 = 1.0 + ( x + x2 * s5 );
	double s9 = 1.0 + ( x + x2 * s7 );
	double s11 = 1.0 + ( x + x2 * s9 );

	//
	//  SET W = DEL(B) - DEL(A + B)
	//
	double t = Math::Pow ( 1.0 / b, 2.0 );

	double w = (((( c5 * s11  * t 
		+ c4 * s9 ) * t
		+ c3 * s7 ) * t
		+ c2 * s5 ) * t
		+ c1 * s3 ) * t
		+ c0;

	w *= ( c / b );
	//
	//  Combine the results.
	//

	double u = d * alnrel ( (a / b) );
	double v = a * ( Math::Log ( b ) - 1.0 );

	// TODO This code doesn't make any sense to me
	// (w - v) - u == (w - u) - v, unless the distributive
	// properties of subtraction and addition have changed
	// see http://people.scs.fsu.edu/~burkardt/f_src/dcdflib/dcdflib.f90
	// if ( v < u )
	// {
	//    algdiv = (w - v) - u;
	// }
	// else
	// {
	//    algdiv = (w - u) - v;
	// }
	// return algdiv;

	return (w - u - v);
}

//****************************************************************************80
//
//  Purpose:
// 
//    ALNREL evaluates the function ln ( 1 + A ).
//
//  Modified:
//
//    17 November 2006
//
//  Reference:
//
//    Armido DiDinato, Alfred Morris,
//    Algorithm 708: 
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//
//  Parameters:
//
//    Input, double %A, the argument.
//
//    Output, double ALNREL, the value of ln ( 1 + A ).
//
double StatClass::alnrel ( double %a )
{
	const double p1 = -0.129418923021993e+01;
	const double p2 =  0.405303492862024e+00;
	const double p3 = -0.178874546012214e-01;
	const double q1 = -0.162752256355323e+01;
	const double q2 =  0.747811014037616e+00;
	const double q3 = -0.845104217945565e-01;

	if ( Math::Abs( a ) <= 0.375 )
	{
		double t = a / ( a + 2.0 );
		double t2 = t * t;
		double w = ((( p3 * t2 + p2 ) * t2 + p1) * t2 + 1.0) / ((( q3 * t2 + q2) * t2 + q1 ) * t2 + 1.0);
		return (2.0 * t * w);
	}
	else
	{
		return Math::Log ( 1.0 + a );
	}
}

//****************************************************************************80
//
//  Purpose:
// 
//    APSER computes the incomplete beta ratio I(SUB(1-X))(B,A).
//
//  Discussion:
//
//    APSER is used only for cases where
//
//      A <= min ( EPS, EPS * B ), 
//      B * X <= 1, and 
//      X <= 0.5.
//
//  Parameters:
//
//    Input, double %A, *B, *X, the parameters of teh
//    incomplete beta ratio.
//
//    Input, double %EPS, a tolerance.
//
//    Output, double APSER, the computed value of the
//    incomplete beta ratio.
//
double StatClass::apser ( double %a, double %b, double %x, double %eps )
{
	const double g = 0.577215664901533;
	double aj = 0.0;
	double c = 0.0;

	double bx = Math::Pow(b,x);
	double t = x-bx;

	if (Math::Pow(b,eps) > 2.e-2)
		c = Math::Log(bx) + g + t;
	else
		c = Math::Log(x) + psi(b) + g + t;


	double tol = Math::Pow(5.0, eps) * Math::Abs(c);
	double j = 1.0;
	double s = 0.0;

	do
	{
		j = j + 1.0;
		t *= (x - bx / j);
		aj = t / j;
		s = s + aj;
	}
	while ( Math::Abs(aj) > tol);

	return ( -( a * ( c + s )));
}

//****************************************************************************80
//
//  Purpose:
// 
//    BCORR evaluates DEL(A0) + DEL(B0) - DEL(A0 + B0).
//
//  Discussion:
//
//    The function DEL(A) is a remainder term that is used in the expression:
//
//      ln ( Gamma ( A ) ) = ( A - 0.5 ) * ln ( A ) 
//        - A + 0.5 * ln ( 2 * PI ) + DEL ( A ),
//
//    or, in other words, DEL ( A ) is defined as:
//
//      DEL ( A ) = ln ( Gamma ( A ) ) - ( A - 0.5 ) * ln ( A ) 
//        + A + 0.5 * ln ( 2 * PI ).
//
//  Parameters:
//
//    Input, double %A0, *B0, the arguments.
//    It is assumed that 8 <= A0 and 8 <= B0.
//
//    Output, double %BCORR, the value of the function.
//
double StatClass::bcorr ( double %a0, double %b0 )
{
	const double c0 =  0.833333333333333e-01;
	const double c1 = -0.277777777760991e-02;
	const double c2 =  0.793650666825390e-03;
	const double c3 = -0.595202931351870e-03;
	const double c4 =  0.837308034031215e-03;
	const double c5 = -0.165322962780713e-02;

	double a = Math::Min ( a0, b0 );
	double b = Math::Max ( a0, b0 );
	double h = a / b;
	double c = h / ( 1.0 + h );
	double x = 1.0 / ( 1.0 + h );
	double x2 = x * x;

	//
	//  SET SN = (1 - X**N)/(1 - X)
	//
	double s3 = 1.0 + ( x + x2 );
	double s5 = 1.0 + ( x + x2 * s3 );
	double s7 = 1.0 + ( x + x2 * s5 );
	double s9 = 1.0 + ( x + x2 * s7 );
	double s11 = 1.0 + ( x + x2 * s9 );

	//
	//  SET W = DEL(B) - DEL(A + B)
	//
	double t = Math::Pow ( 1.0 / b, 2.0 );

	double w = (((( c5 * s11  * t + c4 * s9 ) * t + c3 * s7 ) * t + c2 * s5 ) * t + c1 * s3 ) * t + c0;
	w *= ( c / b );

	//
	//  COMPUTE  DEL(A) + W
	//
	t = Math::Pow ( 1.0 / a, 2.0 );

	return ((((( c5 * t + c4 ) * t + c3 ) * t + c2 ) * t + c1 ) * t + c0 ) / a + w;
}

//****************************************************************************80**
//
//  Purpose:
//
//    BETA evaluates the beta function.
//
//  Modified:
//
//    03 December 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the arguments of the beta function.
//
//    Output, double BETA, the value of the beta function.
//
double StatClass::beta ( double a, double b )
{
	return ( Math::Exp ( beta_log ( a, b ) ) );
}

//*****************************************************************************
//
//  Purpose:
//
//    BETA_ASYM computes an asymptotic expansion for IX(A,B), for large A and B.
//
//  Parameters:
//
//    Input, double %A, *B, the parameters of the function.
//    A and B should be nonnegative.  It is assumed that both A and B
//    are greater than or equal to 15.
//
//    Input, double %LAMBDA, the value of ( A + B ) * Y - B.
//    It is assumed that 0 <= LAMBDA.
//
//    Input, double %EPS, the tolerance.
//
//*****************************************************************************
double StatClass::beta_asym ( double %a, double %b, double %lambda, double %eps )
{
	const double e0 = 1.12837916709551;
	const double e1 = .353553390593274;
	const int num = 20;
	//
	//  NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
	//            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
	//            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
	//     E0 = 2/SQRT(PI)
	//     E1 = 2**(-3/2)
	//
	int K3 = 1;

	array<double> ^a0 = gcnew array<double>(21);
	array<double> ^b0 = gcnew array<double>(21);
	array<double> ^c = gcnew array<double>(21);
	array<double> ^d = gcnew array<double>(21);

	double h = 0.0;
	double r0 = 0.0;
	double r1 = 0.0;
	double w0 = 0.0;

	if(a >= b)
	{
		h = b / a;
		r0 = 1.0 / (1.0+h);
		r1 = (b-a) / a;
		w0 = 1.0 / Math::Sqrt(b * (1.0 + h));
	}
	else
	{
		h = a / b;
		r0 = 1.0 /(1.0+h);
		r1 = (b-a) / b;
		w0 = 1.0 / Math::Sqrt(a * (1.0 + h));
	}

	double T1 = -(lambda / a);
	double T2 = lambda / b;

	double f = a * rlog1(T1) + b * rlog1(T2);
	double t = Math::Exp(-f);

	if (t == 0.0)
		return 0.0;

	double z0 = Math::Sqrt(f);
	double z = 0.5 * (z0 / e1);
	double z2 = f + f;

	a0[0] = 2.0/3.0*r1;
	c[0] = -(0.5*a0[0]);
	d[0] = -c[0];

	double j0 = 0.5 / e0 * error_fc ( K3, z0 );
	double j1 = e1;
	double sum = j0 + d[0] * w0 * j1;
	double s = 1.0;
	double h2 = h * h;
	double hn = 1.0;
	double w = w0;
	double znm1 = z;
	double zn = z2;

	for (int n = 2; n <= num; n += 2 )
	{
		hn = h2 * hn;
		a0[n-1] = 2.0 * r0 * (1.0 + h * hn) / ((double)n + 2.0);
		int np1 = n + 1;
		s += hn;
		a0[np1-1] = 2.0 * r1 * s / ((double)n + 3.0);
		for (int i = n; i <= np1; i++ )
		{
			double r = -(0.5 * ((double)i + 1.0));
			b0[0] = r * a0[0];
			for (int m = 2; m <= i; m++ )
			{
				double bsum = 0.0;
				int mm1 = m - 1;
				for (int j = 1; j <= mm1; j++ )
				{
					int mmj = m - j;
					bsum += (((double)j * r - (double)mmj) * a0[j - 1] * b0[mmj - 1]);
				}
				b0[m-1] = r * a0[m - 1] + bsum / (double)m;
			}
			c[i-1] = b0[i-1] / ((double)i + 1.0);
			double dsum = 0.0;
			int im1 = i - 1;
			for (int j = 1; j <= im1; j++ )
			{
				int imj = i-j;
				dsum += (d[imj - 1] * c[j - 1]);
			}
			d[i-1] = -(dsum + c[i - 1]);
		}
		j0 = e1 * znm1 + ((double)n - 1.0) * j0;
		j1 = e1 * zn + (double)n * j1;
		znm1 = z2 * znm1;
		zn = z2 * zn;
		w = w0 * w;
		double t0 = d[n-1] * w * j0;
		w = w0 * w;
		double t1 = d[np1-1] * w * j1;
		sum += (t0 + t1);

		if(Math::Abs(t0) + Math::Abs(t1) <= eps * sum)
			break;
	}

	double u = Math::Exp(-bcorr(a,b));
	return e0 * t * u * sum;;
}


//****************************************************************************80
//
//  Purpose:
//
//    BETA_FRAC evaluates a continued fraction expansion for IX(A,B).
//
//  Parameters:
//
//    Input, double *A, *B, the parameters of the function.
//    A and B should be nonnegative.  It is assumed that both A and
//    B are greater than 1.
//
//    Input, double *X, *Y.  X is the argument of the
//    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
//
//    Input, double *LAMBDA, the value of ( A + B ) * Y - B.
//
//    Input, double *EPS, a tolerance.
//
//    Output, double BETA_FRAC, the value of the continued
//    fraction approximation for IX(A,B).
//
//****************************************************************************80
double StatClass::beta_frac ( double %a, double %b, double %x, double %y, double %lambda, double %eps )
{
	double bfrac = beta_rcomp ( a, b, x, y );

	if ( bfrac == 0.0 ) 
	{
		return bfrac;
	}

	double c = 1.0 + lambda;
	double c0 = b / a;
	double c1 = 1.0 + 1.0 / a;
	double yp1 = y + 1.0;
	double n = 0.0;
	double p = 1.0;
	double s = a + 1.0;
	double an = 0.0;
	double bn = 1.0;
	double anp1 = 1.0;
	double bnp1 = c / c1;
	double r = c1 / c;

	//
	//  CONTINUED FRACTION CALCULATION
	//
	double r0 = 0.0;
	do
	{
		double t = ++n / a;
		double w = n * (b - n) * x;
		double e = a / s;
		double alpha = p * (p + c0) * e * e * (w * x);
		e = (1.0 + t) / (c1 + t + t);
		double beta = n + w / s + e * (c + n * yp1);
		p = 1.0 + t;
		s += 2.0;

		//
		//  UPDATE AN, BN, ANP1, AND BNP1
		//
		t = alpha * an + beta * anp1;
		an = anp1;
		anp1 = t;
		t = alpha * bn + beta * bnp1;
		bn = bnp1;
		bnp1 = t;
		r0 = r;
		r = anp1 / bnp1;

		//
		//  RESCALE AN, BN, ANP1, AND BNP1
		//
		an /= bnp1;
		bn /= bnp1;
		anp1 = r;
		bnp1 = 1.0;
	}
	while (Math::Abs(r - r0) > eps * r);

	//
	//  TERMINATION
	//
	return bfrac * r;
}

//****************************************************************************80
//
//  Purpose:
// 
//    BETA_GRAT evaluates an asymptotic expansion for IX(A,B).
//
//  Parameters:
//
//    Input, double %A, *B, the parameters of the function.
//    A and B should be nonnegative.  It is assumed that 15 <= A 
//    and B <= 1, and that B is less than A.
//
//    Input, double %X, *Y.  X is the argument of the
//    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
//
//    Input/output, double %W, a quantity to which the
//    result of the computation is to be added on output.
//
//    Input, double %EPS, a tolerance.
//
//    Output, int %IERR, an error flag, which is 0 if no error
//    was detected.
//
//****************************************************************************80
void StatClass::beta_grat ( double %a, double %b, double %x, double %y, double %w, double %eps,int %ierr )
{
	double lnx, p, q;
	int i,n,nm1;
	array<double> ^c = gcnew array<double>(30);
	array<double> ^d = gcnew array<double>(30);

	double bm1 = b - 1.0;
	double nu = a + 0.5 * bm1;

	if (y > 0.375)
		lnx = Math::Log(x);
	else
		lnx = alnrel(-y);

	double z = -(nu * lnx);

	if ( b * z == 0.0 )
	{
		ierr = 1;
		return;
	}

	//
	//  COMPUTATION OF THE EXPANSION
	//  SET R = EXP(-Z)*Z**B/GAMMA(B)
	//
	double r = b * (1.0 + gam1(b)) * Math::Exp(b * Math::Log(z));
	r *= (Math::Exp(a * lnx) * Math::Exp(0.5 * bm1 * lnx));

	double u = algdiv(b, a) + b * Math::Log(nu);
	u = r * Math::Exp( -u );

	if ( u == 0.0 )
	{
		ierr = 1;
		return;
	}

	gamma_rat1 ( b, z, r, p, q, eps );
	double v = Math::Pow(1.0 / nu, 2.0) / 4.0;
	double t2 = lnx * lnx / 4.0;
	double l = w / u;
	double j = q / r;
	double sum = j;
	double t = 1.0;
	double cn = 1.0;
	double n2 = 0.0;
	for ( n = 1; n <= 30; n++ )
	{
		double bp2n = b + n2;
		j = (bp2n * (bp2n + 1.0) * j + (z + bp2n + 1.0) * t) * v;
		n2 = n2 + 2.0;
		t *= t2;
		cn /= (n2 * (n2 + 1.0));
		c[n - 1] = cn;
		double s = 0.0;
		if ( n != 1 )
		{
			nm1 = n-1;
			double coef = b - (double)n;
			for ( i = 1; i <= nm1; i++ )
			{
				s = s + (coef * c[i - 1] * d[n - i - 1]);
				coef = coef + b;
			}
		}

		d[n-1] = bm1 * cn + s / (double)n;
		double dj = d[n-1] * j;
		sum += dj;
		if ( sum <= 0.0 )
		{
			ierr = 1;
			return;
		}
		if(Math::Abs(dj) <= eps * (sum + l))
			break;
	}
	//
	//  ADD THE RESULTS TO W
	//
	ierr = 0;
	w = w + (u * sum);
	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    BETA_INC evaluates the incomplete beta function IX(A,B).
//
//  Author:
//
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//
//  Parameters:
//
//    Input, double %A, *B, the parameters of the function.
//    A and B should be nonnegative.
//
//    Input, double %X, *Y.  X is the argument of the
//    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
//
//    Output, double %W, *W1, the values of IX(A,B) and
//    1-IX(A,B).
//
//    Output, int %IERR, the error flag.
//    0, no error was detected.
//    1, A or B is negative;
//    2, A = B = 0;
//    3, X < 0 or 1 < X;
//    4, Y < 0 or 1 < Y;
//    5, X + Y /= 1;
//    6, X = A = 0;
//    7, Y = B = 0.
//
//****************************************************************************80
void StatClass::beta_inc ( double %a, double %b, double %x, double %y, double %w, double %w1, int %ierr )
{
	int ierr1 = 0;
	double lambda = 0.0;

	//
	//  EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
	//  NUMBER FOR WHICH 1.0 + EPS .GT. 1.0
	//
	// Pseudo code to determine eps
	//{
	//	double eps = 1.0;
	//	S10:
	//	eps = eps/2.0;
	//	if (eps + 1.0 > 1.0) goto S10;
	//	eps *= 2.0;
	//}


	double eps = 2.22044604925031E-16;
	w = 0.0;
	w1 = 0.0;

	if (a < 0.0 || b < 0.0)
	{
		//  A or B is negative;
		ierr = 1;
		return;
	};

	if (a == 0.0 && b == 0.0)
	{
		// A = B = 0
		ierr = 2;
		return;
	}

	if (x < 0.0 || x > 1.0)
	{
		// X < 0 or 1 < X;
		ierr = 3;
		return;
	}

	if (y < 0.0 || y > 1.0)
	{
		// Y < 0 or 1 < Y
		ierr = 4;
		return;
	}

	double z = x + y - 1.0;

	if (Math::Abs(z) > 3.0 * eps)
	{
		// X + Y /= 1
		ierr = 5;
		return;
	}

	ierr = 0;

	// If we get to this point, we know that x + y == 1
	// a >= 0 && b >= 0
	// a != 0 && b != 0
	if (x == 0.0)
	{
		if (a == 0.0)
		{
			ierr = 6;
			return;
		}
		else
		{
			w = 0.0;
			w1 = 1.0;
			return;
		}
	}

	// If we get to this point, we know that x != 0
	if ((y == 0.0) || (a == 0.0))
	{
		if (b == 0.0)
		{
			ierr = 7;
			return;
		}
		else
		{
			w = 1.0;
			w1 = 0.0;
			return;
		}
	}

	// If we get to this point, we know that x !=0, y!= 0 and a != 0
	if (b == 0.0)
	{
		w = 0.0;
		w1 = 1.0;
		return;
	}


	eps = Math::Max(eps,1.e-15);

	if(Math::Max(a,b) < 1.e-3 * eps)
	{
		//
		//  PROCEDURE FOR A AND B .LT. 1.E-3*EPS
		//
		w = b/(a+b);
		w1 = a/(a+b);
		return;
	}

	bool swap = false;
	int n;
	double a0 = a;
	double b0 = b;
	double x0 = x;
	double y0 = y;

	if(Math::Min(a0,b0) > 1.0)
	{
		//
		//  PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1
		//
		if (a > b)
		{
			lambda = (a + b) * y - b;
		}
		else
		{
			lambda = a - (a + b) * x;
		}

		if (lambda < 0.0)
		{
			swap = true;
			a0 = b;
			b0 = a;
			x0 = y;
			y0 = x;
			lambda = Math::Abs(lambda);
		}

		if (b0 < 40.0 && b0 * x0 <= 0.7)
		{
			w = beta_pser(a0, b0, x0, eps);
			w1 = 1.0 - w;
		}
		else
		{
			if (b0 < 40.0)
			{
				n = ( int ) b0;
				b0 -= (double)n;
				if (b0 == 0.0)
				{
					n -= 1;
					b0 = 1.0;
				}
				w = beta_up ( b0, a0, y0, x0, n, eps );
				if (x0 <= 0.7)
				{
					w = w + beta_pser(a0, b0, x0, eps);
					w1 = 1.0 - w;
				}
				else
				{
					if(a0 <= 15.0)
					{
						n = 20;
						w = w + beta_up ( a0, b0, x0, y0, n, eps );
						a0 = a0 + (double)n;
					}
					// beta_grat ( a0, b0, x0, y0, w, 15.0*eps, ierr1 );
					beta_grat ( a0, b0, x0, y0, w, 5.0 * eps, ierr1 );
					w1 = 1.0 - w;
				}
			}
			else
			{
				if ((b0 <= 100.0) || (lambda > 0.03 * b0) || (a0 <= 100.0) || (lambda > 0.03 * a0))
				{
					//w = beta_frac ( a0, b0, x0, y0, lambda, 15.0*eps );
					w = beta_frac ( a0, b0, x0, y0, lambda, 5.0 * eps );
					w1 = w1 = 1.0 - w;
				}
				else
				{
					w = beta_asym ( a0, b0, lambda, 100.0 * eps );
					w1 = 1.0 - w;
				}
			}
		}
	}
	else
	{
		//
		//  PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1
		//
		if (x > 0.5)
		{
			swap = true;
			a0 = b;
			b0 = a;
			x0 = y;
			y0 = x;
		}

		if(b0 < Math::Min(eps, eps*a0))
		{
			w = fpser(a0, b0, x0, eps);
			w1 = 1.0 - w;
		}
		else
		{
			if ((a0 < Math::Min(eps, eps * b0)) && (b0 * x0 <= 1.0))
			{
				w1 = apser(a0, b0, x0, eps);
				w = 1.0 - w1;
			}
			else
			{
				if (Math::Max(a0, b0) <= 1.0)
				{
					if ((a0 >= Math::Min(0.2, b0)) || (Math::Pow(x0, a0) <= 0.9))
					{
						w = beta_pser(a0, b0, x0, eps);
						w1 = 1.0 - w;
					}
					else
					{
						if (x0 >= 0.3)
						{
							w1 = beta_pser(b0, a0, y0, eps);
							w = 1.0 - w1;
						}
						else
						{
							n = 20;
							w1 = beta_up ( b0, a0, y0, x0, n, eps );
							b0 = b0 + (double)n;
							beta_grat( b0, a0, y0, x0, w1, 15.0 * eps, ierr1);
							w = 1.0 - w1;
						}
					}
				}
				else
				{
					if ( (b0 <= 1.0) || ((x0 < 0.1) && (Math::Pow(x0 * b0, a0) <= 0.7)) )
					{
						w = beta_pser(a0, b0, x0, eps);
						w1 = 1.0 - w;
					}
					else
					{
						if (x0 >= 0.3)
						{
							w1 = beta_pser(b0, a0, y0, eps);
							w = 1.0 - w1;
						}
						else
						{
							if ( b0 <= 15.0 )
							{
								n = 20;
								w1 = beta_up ( b0, a0, y0, x0, n, eps );
								b0 = b0 + (double)n;
							}
							beta_grat( b0, a0, y0, x0, w1, 15.0 * eps, ierr1);
							w = 1.0 - w1;
						}
					}
				}
			}
		}
	}

	if (swap)
	{
		double t = w;
		w = w1;
		w1 = t;
	}
	return;
}

//******************************************************************************
//
//  Purpose: 
//
//    BETA_INC_VALUES returns some values of the incomplete Beta function.
//
//  Discussion:
//
//    The incomplete Beta function may be written
//
//      BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
//                      / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
//
//    Thus,
//
//      BETA_INC(A,B,0.0) = 0.0
//      BETA_INC(A,B,1.0) = 1.0
//
//    Note that in Mathematica, the expressions:
//
//      BETA[A,B]   = Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
//      BETA[X,A,B] = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
//
//    and thus, to evaluate the incomplete Beta function requires:
//
//      BETA_INC(A,B,X) = BETA[X,A,B] / BETA[A,B]
//
//  Modified:
//
//    09 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//    Karl Pearson,
//    Tables of the Incomplete Beta Function,
//    Cambridge University Press, 1968.
//
//  Parameters:
//
//    Input/output, int %N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double %A, *B, the parameters of the function.
//
//    Output, double %X, the argument of the function.
//
//    Output, double %FX, the value of the function.
//
//******************************************************************************
void StatClass::beta_inc_values ( int %n_data, double %a, double %b, double %x, double %fx )
{
	const int N_MAX = 30;

	array<double> ^a_vec =
	{
		0.5,  0.5,   0.5,  1.0, 
		1.0,  1.0,   1.0,  1.0, 
		2.0,  2.0,   2.0,  2.0, 
		2.0,  2.0,   2.0,  2.0, 
		2.0,  5.5,  10.0, 10.0, 
		10.0, 10.0, 20.0, 20.0, 
		20.0, 20.0, 20.0, 30.0, 
		30.0, 40.0 
	};

	array<double> ^b_vec =
	{ 
		0.5,  0.5,  0.5,  0.5, 
		0.5,  0.5,  0.5,  1.0, 
		2.0,  2.0,  2.0,  2.0, 
		2.0,  2.0,  2.0,  2.0, 
		2.0,  5.0,  0.5,  5.0, 
		5.0, 10.0,  5.0, 10.0, 
		10.0, 20.0, 20.0, 10.0, 
		10.0, 20.0 
	};

	array<double> ^fx_vec =
	{ 
		0.0637686, 0.2048328, 1.0000000, 0.0,       
		0.0050126, 0.0513167, 0.2928932, 0.5000000, 
		0.028,     0.104,     0.216,     0.352,     
		0.500,     0.648,     0.784,     0.896,     
		0.972,     0.4361909, 0.1516409, 0.0897827, 
		1.0000000, 0.5000000, 0.4598773, 0.2146816, 
		0.9507365, 0.5000000, 0.8979414, 0.2241297, 
		0.7586405, 0.7001783 
	};

	array<double> ^x_vec =
	{ 
		0.01, 0.10, 1.00, 0.0,  
		0.01, 0.10, 0.50, 0.50, 
		0.1,  0.2,  0.3,  0.4,  
		0.5,  0.6,  0.7,  0.8,  
		0.9,  0.50, 0.90, 0.50, 
		1.00, 0.50, 0.80, 0.60, 
		0.80, 0.50, 0.60, 0.70, 
		0.80, 0.70
	};

	// NOTE Use the ternary operator to set the value of n_data
	n_data = ( n_data < 0 ) ? 0 : n_data;

	if ( N_MAX <= n_data)
	{
		// n_data wraps around to 0 (no lookup values)
		n_data = 0;
		a = 0.0;
		b = 0.0;
		x = 0.0;
		fx = 0.0;
	}
	else
	{
		a = a_vec[n_data];
		b = b_vec[n_data];
		x = x_vec[n_data];
		fx = fx_vec[n_data];
		n_data++;
	}

}

//*****************************************************************************
//
//  Purpose:
//
//    BETA_LOG evaluates the logarithm of the beta function.
//
//  Reference:
//
//    Armido DiDinato and Alfred Morris,
//    Algorithm 708: 
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software, 
//    Volume 18, 1993, pages 360-373.
//
//  Parameters:
//
//    Input, double %A0, *B0, the parameters of the function.
//    A0 and B0 should be nonnegative.
//
//    Output, double %BETA_LOG, the value of the logarithm
//    of the Beta function.
//
//*****************************************************************************
double StatClass::beta_log ( double %a0, double %b0 )
{
	const double e = .918938533204673;
	double value, h, w, z;
	int i,n;

	double a = Math::Min(a0, b0);
	double b = Math::Max(a0, b0);

	if(a < 8.0)
	{
		if (a < 1.0)
		{
			//
			//  PROCEDURE WHEN A .LT. 1
			//
			if(b >= 8.0)
				value = gamma_log ( a ) + algdiv(a, b);
			else
				value = gamma_log ( a ) + ( gamma_log ( b ) - gamma_log ( a + b ) );
			return value;
		}
		//
		//  PROCEDURE WHEN 1 .LE. A .LT. 8
		//
		if (a <= 2.0)
		{
			if (b <= 2.0)
			{
				value = gamma_log(a) + gamma_log(b) - gsumln(a, b);
				return value;
			}
			else
			{
				w = 0.0;
				if (b < 8.0)
				{
					//
					//  REDUCTION OF B WHEN B .LT. 8
					//
					n = (int) (b - 1.0);
					z = 1.0;
					for ( i = 1; i <= n; i++ )
					{
						b -= 1.0;
						z *= (b / (a + b));
					}
					value = w + Math::Log(z) + (gamma_log(a) + ( gamma_log(b) - gsumln(a, b)));
					return value;
				}
				else
				{
					value = gamma_log(a) + algdiv(a, b);
					return value;
				}
			}
		}

		if (b > 1000.0)
		{
			//
			//  REDUCTION OF A WHEN B .GT. 1000
			//
			n = ( int ) ( a - 1.0 );
			w = 1.0;
			for ( i = 1; i <= n; i++ )
			{
				a -= 1.0;
				w *= (a / (1.0 + a / b));
			}
			value = Math::Log(w) - (double) n * Math::Log(b) + (gamma_log(a) + algdiv(a, b));
			return value;
		}
		else
		{
			//
			//  REDUCTION OF A WHEN B .LE. 1000
			//
			n = (int) ( a - 1.0 );
			w = 1.0;
			for ( i = 1; i <= n; i++ )
			{
				a -= 1.0;
				h = a / b;
				w *= (h / (1.0 + h));
			}
			w = Math::Log(w);
			if (b < 8.0)
			{
				//
				//  REDUCTION OF B WHEN B .LT. 8
				//
				n = ( int ) ( b - 1.0 );
				z = 1.0;
				for ( i = 1; i <= n; i++ )
				{
					b -= 1.0;
					z *= (b / (a + b));
				}
				value = w + Math::Log(z)+ (gamma_log(a) + (gamma_log(b) - gsumln(a, b)));
				return value;
			}
			else
			{
				value = w + gamma_log(a) + algdiv(a, b);
				return value;
			}
		}
	}
	else
	{
		//
		//  PROCEDURE WHEN A .GE. 8
		//
		w = bcorr(a, b);
		h = a / b;
		double c = h / (1.0 + h);
		double u = -((a - 0.5) * Math::Log(c));
		double v = b * alnrel(h);
		return -(0.5 * Math::Log(b)) + e + w - v - u;
	}
}

//****************************************************************************80
//
//  Purpose:
// 
//    BETA_PSER uses a power series expansion to evaluate IX(A,B)(X).
//
//  Discussion:
//
//    BETA_PSER is used when B <= 1 or B*X <= 0.7.
//
//  Parameters:
//
//    Input, double %A, *B, the parameters.
//
//    Input, double %X, the point where the function
//    is to be evaluated.
//
//    Input, double %EPS, the tolerance.
//
//    Output, double BETA_PSER, the approximate value of IX(A,B)(X).
//
//****************************************************************************80
double StatClass::beta_pser ( double %a, double %b, double %x, double %eps )
{
	// Removed from definition
	double bpser,a0,apb,b0,c,n,sum,t,tol,u,w,z;
	int i,m;

	bpser = 0.0;
	if(x == 0.0) return bpser;
	//
	//  COMPUTE THE FACTOR X**A/(A*BETA(A,B))
	//
	a0 = Math::Min(a,b);
	if (a0 >= 1.0)
	{
		z = a*Math::Log(x) - beta_log(a, b);
		bpser = Math::Exp(z) / a;
	}
	else
	{
		b0 = Math::Max(a, b);
		if(b0 < 8.0)
		{
			if(b0 <= 1.0)
			{
				//
				//  PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1
				//
				bpser = Math::Pow(x, a);
				if(bpser == 0.0) return bpser;
				apb = a + b;
				if(apb <= 1.0)
					z = 1.0 + gam1(apb);
				else
				{
					u = a + b - 1.0;
					z = (1.0 + gam1(u)) / apb;
				}
				c = (1.0 + gam1(a)) * (1.0 + gam1(b)) / z;
				bpser *= (c * (b / apb));
			}
			else
			{
				//
				//  PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8
				//
				u = gamma_ln1 ( a0 );
				m = ( int ) ( b0 - 1.0 );
				if (m >= 1)
				{
					c = 1.0;
					for ( i = 1; i <= m; i++ )
					{
						b0 -= 1.0;
						c *= (b0 / (a0 + b0));
					}
					u = Math::Log(c) + u;
				}
				z = a * Math::Log(x) - u;
				b0 -= 1.0;
				apb = a0 + b0;
				if(apb <= 1.0)
					t = 1.0 + gam1(apb);
				else
				{
					u = a0 + b0 - 1.;
					t = (1.0 + gam1(u)) / apb;
				}
				bpser = Math::Exp(z) * (a0 / a) * (1.0 + gam1(b0)) / t;
			}
		}
		else
		{
			//
			//  PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8
			//
			u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 );
			z = a * Math::Log(x) - u;
			bpser = a0 / a * Math::Exp(z);
		}
	}

	if(bpser == 0.0 || a <= 0.1 * eps) return bpser;
	//
	//  COMPUTE THE SERIES
	//
	sum = 0.0;
	n = 0.0;
	c = 1.0;
	tol = eps / a;
	do
	{
		n = n + 1.0;
		c *= ((0.5 + (0.5-b / n)) * x);
		w = c / (a + n);
		sum = sum + w;
	} while (Math::Abs(w) > tol);

	bpser *= (1.0 + a * sum);
	return bpser;
}

//****************************************************************************80
//
//  Purpose:
// 
//    BETA_RCOMP evaluates X**A * Y**B / Beta(A,B).
//
//  Parameters:
//
//    Input, double %A, *B, the parameters of the Beta function.
//    A and B should be nonnegative.
//
//    Input, double %X, *Y, define the numerator of the fraction.
//
//    Output, double BETA_RCOMP, the value of X**A * Y**B / Beta(A,B).
//
//****************************************************************************80
double StatClass::beta_rcomp ( double %a, double %b, double %x, double %y )
{
	double const Const = .398942280401433;
	double brcomp, a0, apb, b0, c, e, h, lambda, lnx, lny, t, u, v, x0, y0, z;
	int n;

	//
	//  CONST = 1/SQRT(2*PI)
	//
	brcomp = 0.0;
	if(x == 0.0 || y == 0.0) return brcomp;
	a0 = Math::Min(a,b);
	if (a0 < 8.0)
	{
		if (x <= 0.375)
		{
			lnx = Math::Log(x);
			lny = alnrel(-x);
		}
		else if (y <= 0.375)
		{
			lnx = alnrel(-y);
			lny = Math::Log(y);
		}
		else
		{
			lnx = Math::Log(x);
			lny = Math::Log(y);
		}
		z = a*lnx+b*lny;
		if (a0 >= 1.0)
		{
			z -= beta_log(a,b);
			brcomp = Math::Exp(z);
			return brcomp;
		}
		//
		//  PROCEDURE FOR A .LT. 1 OR B .LT. 1
		//
		b0 = Math::Max(a,b);
		if (b0 < 8.0)
		{
			if (b0 <= 1.0)
			{
				//
				//  ALGORITHM FOR B0 .LE. 1
				//
				brcomp = Math::Exp(z);

				if (brcomp == 0.0)
					return brcomp;

				apb = a + b;

				if (apb <= 1.0)
				{
					z = 1.0+gam1(apb);
				}
				else
				{
					u = a+b-1.;
					z = (1.0+gam1(u))/apb;
				}
				c = (1.0 + gam1(a)) * (1.0 + gam1(b)) / z;
				brcomp = brcomp * (a0 * c) / (1.0 + a0 / b0);
				return brcomp;
			}
			//
			//  ALGORITHM FOR 1 .LT. B0 .LT. 8
			//
			u = gamma_ln1 ( a0 );
			n = ( int ) ( b0 - 1.0 );
			if (n >= 1)
			{
				c = 1.0;
				for ( int i = 1; i <= n; i++ )
				{
					b0 -= 1.0;
					c *= (b0 / (a0 + b0));
				}
				u = Math::Log(c) + u;
			}
			z -= u;
			b0 -= 1.0;
			apb = a0 + b0;
			if (apb <= 1.0)
			{
				t = 1.0 + gam1(apb);
			}
			else
			{
				u = a0 + b0 - 1.0;
				t = (1.0 + gam1(u)) / apb;
			}
			brcomp = a0*Math::Exp(z)*(1.0+gam1(b0))/t;
			return brcomp;
		}
		//
		//  ALGORITHM FOR B0 .GE. 8
		//
		u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 );
		brcomp = a0 * Math::Exp(z - u);
		return brcomp;
	}
	//
	//  PROCEDURE FOR A .GE. 8 AND B .GE. 8
	//
	if (a <= b)
	{
		h = a / b;
		x0 = h / (1.0 + h);
		y0 = 1.0 / (1.0 + h);
		lambda = a - (a+b) * x;
	}
	else
	{
		h = b / a;
		x0 = 1.0 / (1.0 + h);
		y0 = h / (1.0 + h);
		lambda = (a + b) * y - b;
	}
	e = -(lambda / a);
	if (Math::Abs(e) <= 0.6)
	{
		u = rlog1(e);
	}
	else
	{
		u = e - Math::Log(x / x0);
	}
	e = lambda / b;

	if (Math::Abs(e) <= 0.6)
	{
		v = rlog1(e);
	}
	else
	{
		v = e - Math::Log(y / y0);
	}
	z = Math::Exp( - (a * u + b * v));
	brcomp = Const * Math::Sqrt(b * x0) * z * Math::Exp( -bcorr(a,b) );
	return brcomp;
}

//****************************************************************************80
//
//  Purpose:
// 
//    BETA_RCOMP1 evaluates Math::Exp(MU) * X**A * Y**B / Beta(A,B).
//
//  Parameters:
//
//    Input, int MU, ?
//
//    Input, double A, B, the parameters of the Beta function.
//    A and B should be nonnegative.
//
//    Input, double X, Y, ?
//
//    Output, double BETA_RCOMP1, the value of
//    Math::Exp(MU) * X**A * Y**B / Beta(A,B).
//
//****************************************************************************80
double StatClass::beta_rcomp1 ( int %mu, double %a, double %b, double %x, double %y )
{
	double const Const = .398942280401433e0;
	double brcmp1,a0,apb,b0,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z;

	//
	//     CONST = 1/SQRT(2*PI)
	//

	a0 = Math::Min(a,b);
	if (a0 < 8.0)
	{
		if (x <= 0.375)
		{
			lnx = Math::Log(x);
			lny = alnrel(-x);
		}
		else if (y <= 0.375)
		{
			lnx = alnrel(-y);
			lny = Math::Log(y);
		}
		else
		{
			lnx = Math::Log(x);
			lny = Math::Log(y);
		}
		z = a*lnx+b*lny;

		if (a0 >= 1.0)
		{
			z -= beta_log(a,b);
			brcmp1 = esum(mu,z);
			return brcmp1;
		}
		//
		//   PROCEDURE FOR A .LT. 1 OR B .LT. 1
		//
		b0 = Math::Max(a,b);
		if (b0 < 8.0)
		{
			if (b0 <= 1.0)
			{
				//
				//  ALGORITHM FOR B0 .LE. 1
				//
				brcmp1 = esum(mu,z);
				if (brcmp1 == 0.0)
					return brcmp1;

				apb = a+b;

				if (apb <= 1.0)
				{
					z = 1.0+gam1(apb);
				}
				else
				{
					u = a+b-1.e0;
					z = (1.0+gam1(u))/apb;
				}
				c = (1.0+gam1(a))*(1.0+gam1(b))/z;
				brcmp1 = brcmp1*(a0*c)/(1.0+a0/b0);
				return brcmp1;
			}
			//
			//  ALGORITHM FOR 1 .LT. B0 .LT. 8
			//
			u = gamma_ln1 ( a0 );
			int n = ( int ) ( b0 - 1.0 );

			if (n >= 1)
			{
				c = 1.0;
				for ( int i = 1; i <= n; i++ )
				{
					b0 -= 1.0;
					c *= (b0/(a0+b0));
				}
				u = Math::Log(c)+u;
			}
			z -= u;
			b0 -= 1.0;
			apb = a0 + b0;
			if (apb <= 1.0)
			{
				t = 1.0+gam1(apb);
			}
			else
			{
				u = a0 + b0 - 1.0;
				t = (1.0 + gam1(u)) / apb;
			}
			brcmp1 = a0 * esum(mu,z) * (1.0 + gam1(b0)) / t;
			return brcmp1;
		}
		else
		{
			//
			//  ALGORITHM FOR B0 .GE. 8
			//
			u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 );
			brcmp1 = a0 * esum(mu,z-u);
			return brcmp1;
		}
	}
	//
	//    PROCEDURE FOR A .GE. 8 AND B .GE. 8
	//
	if (a <= b)
	{
		h = a / b;
		x0 = h / (1.0 + h);
		y0 = 1.0 / (1.0 + h);
		lambda = a - (a + b) * x;
	}
	else
	{
		h = b / a;
		x0 = 1.0 / (1.0 + h);
		y0 = h / (1.0 + h);
		lambda = (a + b) * y - b;
	}
	e = -(lambda / a);
	if (Math::Abs(e) <= 0.6)
	{
		u = rlog1(e);
	}
	else
	{
		u = e-Math::Log(x/x0);
	}
	e = lambda / b;
	if (Math::Abs(e) <= 0.6)
	{
		v = rlog1(e);
	}
	else
	{
		v = e-Math::Log(y/y0);
	}
	z = esum(mu, -(a * u + b * v));
	brcmp1 = Const * Math::Sqrt(b * x0) * z * Math::Exp(-bcorr(a, b));
	return brcmp1;
}

//****************************************************************************80
//
//  Purpose:
// 
//    BETA_UP evaluates IX(A,B) - IX(A+N,B) where N is a positive integer.
//
//  Parameters:
//
//    Input, double %A, *B, the parameters of the function.
//    A and B should be nonnegative.
//
//    Input, double %X, *Y, ?
//
//    Input, int %N, ?
//
//    Input, double %EPS, the tolerance.
//
//    Output, double BETA_UP, the value of IX(A,B) - IX(A+N,B).
//
//****************************************************************************80
double StatClass::beta_up ( double %a, double %b, double %x, double %y, int %n, double %eps )
{
	int K1 = 1;
	int K2 = 0;
	//
	//  OBTAIN THE SCALING FACTOR EXP(-MU) AND
	//  EXP(MU)*(X**A*Y**B/BETA(A,B))/A
	//
	double t = 0.0;
	double apb = a + b;
	double ap1 = a + 1.0;
	int mu = 0, k = 0;
	double d = 1.0;
	double r = 0.0;

	if (!((n == 1) || (a < 1.0) || (apb < 1.1 * ap1)))
	{
		mu = ( int ) Math::Abs ( exparg(K1) );
		k = ( int ) exparg ( K2 );
		if (k < mu) mu = k;
		t = mu;
		d = Math::Exp(-t);
	}
	double bup = beta_rcomp1 ( mu, a, b, x, y ) / a;
	if (n == 1 || bup == 0.0)
		return bup;
	int nm1 = n-1;
	double w = d;
	//
	//  LET K BE THE INDEX OF THE MAXIMUM TERM
	//
	k = 0;
	if (b > 1.0)
	{
		// Next two statments were moved from else clause and modified
		// to get rid of goto
		r = (b - 1.0) * x / y - a;
		if (!((y > 1.0e-4) && (r < 1.0)))
		{

			if (y <= 1.0e-4)
			{
				k = nm1;
			}
			else
			{
				t = ( double ) nm1;
				k = nm1;
				if ( r < t ) k = ( int ) r;
			}
			//
			//          ADD THE INCREASING TERMS OF THE SERIES
			//
			for (int  i = 1; i <= k; i++ )
			{
				double l = i - 1;
				d = (apb + l)/(ap1 + l) * x * d;
				w = w + d;
			}
			if (k == nm1)
			{
				bup *= w;
				return bup;
			}
		}
	}
	//
	//          ADD THE REMAINING TERMS OF THE SERIES
	//
	int kp1 = k + 1;
	for (int i = kp1; i <= nm1; i++ )
	{
		double l = i - 1;
		d = (apb + l)/(ap1 + l) * x * d;
		w = w + d;
		if (d <= eps * w) 
			break;
	}
	//
	//  TERMINATE THE PROCEDURE
	//
	bup *= w;
	return bup;
}

//****************************************************************************80***
//
//  Purpose: 
//
//    BINOMIAL_CDF_VALUES returns some values of the binomial CDF.
//
//  Discussion:
//
//    CDF(X)(A,B) is the probability of at most X successes in A trials,
//    given that the probability of success on a single trial is B.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int %N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int %A, double %B, the parameters of the function.
//
//    Output, int %X, the argument of the function.
//
//    Output, double %FX, the value of the function.
//
//****************************************************************************80***
void StatClass::binomial_cdf_values ( int %n_data, int %a, double %b, int %x, double %fx )
{
	array<const int> ^a_vec =
	{
		2,  2,  2,  2, 
		2,  4,  4,  4, 
		4, 10, 10, 10, 
		10, 10, 10, 10, 
		10
	};

	array<const double> ^b_vec =
	{
		0.05, 0.05, 0.05, 0.50, 
		0.50, 0.25, 0.25, 0.25, 
		0.25, 0.05, 0.10, 0.15, 
		0.20, 0.25, 0.30, 0.40, 
		0.50
	};

	array<const double> ^fx_vec =
	{
		0.9025, 0.9975, 1.0000, 0.2500, 
		0.7500, 0.3164, 0.7383, 0.9492, 
		0.9961, 0.9999, 0.9984, 0.9901, 
		0.9672, 0.9219, 0.8497, 0.6331, 
		0.3770
	};

	array<const int> ^x_vec =
	{
		0, 1, 2, 0, 
		1, 0, 1, 2, 
		3, 4, 4, 4, 
		4, 4, 4, 4, 
		4
	};

	const int N_MAX = 17;

	if ( n_data < 0 )
	{
		n_data = 0;
	}

	n_data++;

	if ( N_MAX < n_data )
	{
		n_data = 0;
		a = 0;
		b = 0.0;
		x = 0;
		fx = 0.0;
	}
	else
	{
		a = a_vec[n_data-1];
		b = b_vec[n_data-1];
		x = x_vec[n_data-1];
		fx = fx_vec[n_data-1];
	}
	return;
}

//*****************************************************************************
//
//  Purpose:
// 
//    CDFBET evaluates the CDF of the Beta Distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the beta distribution 
//    given the others.
//
//    The value P of the cumulative distribution function is calculated 
//    directly by code associated with the reference.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The beta density is proportional to t%(A-1) * (1-t)%(B-1).
//
//  Modified:
//
//    09 June 2004
//
//  Reference:
//
//    Armido DiDinato and Alfred Morris,
//    Algorithm 708: 
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software, 
//    Volume 18, 1993, pages 360-373.
//
//  Parameters:
//
//    Input, int %WHICH, indicates which of the next four argument
//    values is to be calculated from the others.
//    1: Calculate P and Q from X, Y, A and B;
//    2: Calculate X and Y from P, Q, A and B;
//    3: Calculate A from P, Q, X, Y and B;
//    4: Calculate B from P, Q, X, Y and A.
//
//    Input/output, double %P, the integral from 0 to X of the
//    chi-square distribution.  Input range: [0, 1].
//
//    Input/output, double %Q, equals 1-P.  Input range: [0, 1].
//
//    Input/output, double %X, the upper limit of integration 
//    of the beta density.  If it is an input value, it should lie in
//    the range [0,1].  If it is an output value, it will be searched for
//    in the range [0,1].
//
//    Input/output, double %Y, equal to 1-X.  If it is an input
//    value, it should lie in the range [0,1].  If it is an output value,
//    it will be searched for in the range [0,1].
//
//    Input/output, double %A, the first parameter of the beta
//    density.  If it is an input value, it should lie in the range
//    (0, +infinity).  If it is an output value, it will be searched
//    for in the range [1D-300,1D300].
//
//    Input/output, double %B, the second parameter of the beta
//    density.  If it is an input value, it should lie in the range
//    (0, +infinity).  If it is an output value, it will be searched
//    for in the range [1D-300,1D300].
//
//    Output, int %STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +4, if X + Y /= 1.
//
//    Output, double %BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
//*****************************************************************************
void StatClass::cdfbet ( int %which, double %p, double %q, double %x, double %y, double %a, double %b, int %status, double %bound )
{
	double K2 = 0.0;
	double K3 = 1.0;
	double K8 = 0.5;
	double K9 = 5.0;
	double fx,xhi,xlo,cum,ccum,xy,pq;
	unsigned long qhi,qleft,qporq = false;
	double T4,T5,T6,T7,T10,T11,T12,T13,T14,T15;

	status = 0;
	bound = 0.0;
	//
	//     Check arguments
	//
	if(which < 1 || which > 4)
	{
		if(which < 1)
			bound = 1.0;
		else
			bound = 4.0e0;
		status = -1;
		return;
	}

	if (which != 1)
	{
		//
		//     P
		//
		if(p < 0.0 || p > 1.0)
		{
			if(p < 0.0)
				bound = 0.0;
			else
				bound = 1.0;

			status = -2;
			return;
		}

		//
		//     Q
		//
		if (q < 0.0 || q > 1.0)
		{
			if (q < 0.0)
				bound = 0.0;
			else
				bound = 1.0;

			status = -3;
			return;
		}

	}

	if (which != 2)
	{
		//
		//     X
		//
		if (x < 0.0 || x > 1.0)
		{
			if (x < 0.0)
				bound = 0.0;
			else
				bound = 1.0;

			status = -4;
			return;
		}

		//
		//     Y
		//
		if (y < 0.0 || y > 1.0)
		{
			if (y < 0.0)
				bound = 0.0;
			else
				bound = 1.0;

			status = -5;
			return;
		}

	}

	if (which != 3)
	{
		//
		//     A
		//
		if (a <= 0.0)
		{
			bound = 0.0;
			status = -6;
			return;
		}
	}

	if (which != 4)
	{
		//
		//     B
		//
		if (b <= 0.0)
		{
			bound = 0.0;
			status = -7;
			return;
		}
	}

	if (which != 1)
	{
		//
		//     P + Q
		//
		pq = p+q;

		// REMINDER dpmpar declarations removed
		if (Math::Abs(pq-0.5e0-0.5e0) > 3.0 * 2.22044604925031E-16 )
		{
			if (pq < 0.0)
				bound = 0.0;
			else
				bound = 1.0;

			status = 3;
			return;
		}
	}

	if (which != 2)
	{
		//
		//     X + Y
		//
		xy = x+y;

		// REMINDER dpmpar declarations removed
		if (Math::Abs(xy-0.5e0-0.5e0) > 3.0 * 2.22044604925031E-16 )
		{
			if (xy < 0.0)
				bound = 0.0;
			else
				bound = 1.0;

			status = 4;
			return;
		}
	}

	if (which != 1)
		qporq = p <= q;
	//
	//     Select the minimum of P or Q
	//     Calculate ANSWERS
	//
	switch (which)
	{
	case 1:
		//
		//     Calculating P and Q
		//
		cumbet(x,y,a,b,p,q);
		status = 0;
		break;
	case 2:
		//
		//     Calculating X and Y
		//
		T4 = atol;
		T5 = tol;
		dstzr(K2,K3,T4,T5);

		if (qporq)
		{
			status = 0;
			dzror(status, y, fx, xlo, xhi, qleft, qhi);
			x = one-y;

			while (status == 1)
			{
				cumbet(x,y,a,b,cum,ccum);
				fx = ccum-q;
				dzror(status,y,fx,xlo,xhi,qleft,qhi);
				x = one-y;
			}
		}
		else
		{
			status = 0;
			dzror(status,x,fx,xlo,xhi,qleft,qhi);
			y = one-x;

			while (status == 1)
			{
				cumbet(x,y,a,b,cum,ccum);
				fx = cum-p;
				dzror(status,x,fx,xlo,xhi,qleft,qhi);
				y = one-x;
			}
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = 0.0;
			}
			else
			{
				status = 2;
				bound = 1.0;
			}
		}
		break;
	case 3:
		//
		//     Computing A
		//
		a = 5.0;
		T6 = zero;
		T7 = inf;
		T10 = atol;
		T11 = tol;
		dstinv(T6,T7,K8,K8,K9,T10,T11);
		status = 0;
		dinvr(status,a,fx,qleft,qhi);

		while (status == 1)
		{
			cumbet(x,y,a,b,cum,ccum);
			if (qporq)
				fx = cum-p;
			else
				fx = ccum-q;
			dinvr(status,a,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = zero;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
		break;
	case 4:
		//
		//     Computing B
		//
		b = 5.0;
		T12 = zero;
		T13 = inf;
		T14 = atol;
		T15 = tol;
		dstinv(T12,T13,K8,K8,K9,T14,T15);
		status = 0;
		dinvr(status,b,fx,qleft,qhi);
		while (status == 1)
		{
			cumbet(x,y,a,b,cum,ccum);
			if (qporq)
				fx = cum-p;
			else
				fx = ccum-q;

			dinvr(status,b,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = zero;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}
	return;
}

//****************************************************************************
//
//  Purpose:
// 
//    CDFBIN evaluates the CDF of the Binomial distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the binomial distribution 
//    given the others.
//
//    The value P of the cumulative distribution function is calculated 
//    directly.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    P is the probablility of S or fewer successes in XN binomial trials,
//    each trial having an individual probability of success of PR.  
//
//  Modified:
//
//    09 June 2004
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.5.24.
//
//  Parameters:
//
//    Input, int %WHICH, indicates which of argument values is to 
//    be calculated from the others.
//    1: Calculate P and Q from S, XN, PR and OMPR;
//    2: Calculate S from P, Q, XN, PR and OMPR;
//    3: Calculate XN from P, Q, S, PR and OMPR;
//    4: Calculate PR and OMPR from P, Q, S and XN.
//
//    Input/output, double %P, the cumulation, from 0 to S,
//    of the binomial distribution.  If P is an input value, it should 
//    lie in the range [0,1].
//
//    Input/output, double %Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double %S, the number of successes observed.
//    Whether this is an input or output value, it should lie in the
//    range [0,XN].
//
//    Input/output, double %XN, the number of binomial trials.
//    If this is an input value it should lie in the range: (0, +infinity).
//    If it is an output value it will be searched for in the
//    range [1.0D-300, 1.0D+300].
//
//    Input/output, double %PR, the probability of success in each 
//    binomial trial.  Whether this is an input or output value, it should
//    lie in the range: [0,1].
//
//    Input/output, double %OMPR, equal to 1-PR.  Whether this is an
//    input or output value, it should lie in the range [0,1].  Also, it should
//    be the case that PR + OMPR = 1.
//
//    Output, int %STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +4, if PR + OMPR /= 1.
//
//    Output, double %BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//



void StatClass::cdfbin ( int %which, double %p, double %q, double %s, double %xn,
						double %pr, double %ompr, int %status, double %bound )

{
	double K2 = 0.0;
	double K3 = 0.5e0;
	double K4 = 5.0e0;
	double K11 = 1.0;
	double fx,xhi,xlo,cum,ccum,pq,prompr;
	unsigned long qhi,qleft,qporq = 0;
	double T5,T6,T7,T8,T9,T10,T12,T13;

	status = 0;
	bound = 0.0;
	//
	//     Check arguments
	//
	if (which < 1 && which > 4)
	{
		if (which < 1)
			bound = 1.0;
		else
			bound = 4.0e0;

		status = -1;
		return;

	}

	if (which != 1)
	{
		//
		//     P
		//
		if (p < 0.0 || p > 1.0)
		{
			if (p < 0.0)
				bound = 0.0;
			else
				bound = 1.0;

			status = -2;
			return;

		}
	}

	if (which != 1)
	{
		//
		//     Q
		//
		if (q < 0.0 || q > 1.0)
		{
			if (q < 0.0)
				bound = 0.0;
			else
				bound = 1.0;

			status = -3;
			return;
		}
	}

	if (which != 3)
	{
		//
		//     XN
		//
		if (xn <= 0.0)
		{
			bound = 0.0;
			status = -5;
			return;
		}
	}

	if (which != 2)
	{
		//
		//     S
		//

		if (s < 0.0 || which != 3 && s > xn)
		{
			if (s < 0.0)
				bound = 0.0;
			else
				bound = xn;

			status = -4;
			return;
		}
	}

	if (which != 4)
	{
		//
		//     PR
		//
		if (pr < 0.0 || pr > 1.0)
		{
			if (pr < 0.0)
				bound = 0.0;
			else
				bound = 1.0;

			status = -6;
			return;
		}
	}

	if (which != 4)
	{
		//
		//     OMPR
		//

		if (ompr < 0.0 || ompr > 1.0)
		{
			if (ompr < 0.0)
				bound = 0.0;
			else
				bound = 1.0;

			status = -7;
			return;
		}
	}

	if (which != 1)
	{
		//
		//     P + Q
		//
		pq = p + q;
		if (Math::Abs(pq - 1.0) > 3.0 * 2.22044604925031E-16 )
		{
			if (pq < 0.0)
				bound = 0.0;
			else
				bound = 1.0;

			status = 3;
			return;
		}
	}

	if (which != 4)
	{
		//
		//     PR + OMPR
		//
		prompr = pr + ompr;
		if (Math::Abs(prompr - 1.0) > 3.0 * 2.22044604925031E-16 )
		{
			if (prompr < 0.0)
				bound = 0.0;
			else
				bound = 1.0;

			status = 4;
			return;
		}
	}

	if (which != 1)
		qporq = p <= q;
	//
	//     Select the minimum of P or Q
	//     Calculate ANSWERS
	//
	if (1 == which)
	{
		//
		//     Calculating P
		//
		cumbin(s,xn,pr,ompr,p,q);
		status = 0;
	}
	else if (2 == which)
	{
		//
		//     Calculating S
		//
		s = 5.0e0;
		T5 = atol;
		T6 = tol;
		dstinv(K2,xn,K3,K3,K4,T5,T6);
		status = 0;
		dinvr(status,s,fx,qleft,qhi);

		while  (status ==1)
		{
			cumbin(s,xn,pr,ompr,cum,ccum);

			if (qporq)
				fx = cum-p;
			else
				fx = ccum-q;

			dinvr(status,s,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = 0.0;
			}
			else
			{
				status = 2;
				bound = xn;
			}
		}
	}
	else if (3 == which)
	{
		//
		//     Calculating XN
		//
		xn = 5.0e0;
		T7 = zero;
		T8 = inf;
		T9 = atol;
		T10 = tol;
		dstinv(T7,T8,K3,K3,K4,T9,T10);
		status = 0;
		dinvr(status,xn,fx,qleft,qhi);

		while (status == 1)
		{
			cumbin(s,xn,pr,ompr,cum,ccum);
			if (qporq)
				fx = cum-p;
			else
				fx = ccum-q;

			dinvr(status,xn,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = zero;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}
	else if (4 == which)
	{
		//
		//     Calculating PR and OMPR
		//
		T12 = atol;
		T13 = tol;
		dstzr(K2,K11,T12,T13);
		if (qporq)
		{
			status = 0;
			dzror(status,pr,fx,xlo,xhi,qleft,qhi);
			ompr = one-pr;

			while (status == 1)
			{
				cumbin(s,xn,pr,ompr,cum,ccum);
				fx = cum-p;
				dzror(status,pr,fx,xlo,xhi,qleft,qhi);
				ompr = one-pr;
			}
		}
		else
		{
			status = 0;
			dzror(status,ompr,fx,xlo,xhi,qleft,qhi);
			pr = one-ompr;

			while (status == 1)
			{
				cumbin(s,xn,pr,ompr,cum,ccum);
				fx = ccum-q;
				dzror(status,ompr,fx,xlo,xhi,qleft,qhi);
				pr = one-ompr;
			}
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = 0.0;
			}
			else
			{
				status = 2;
				bound = 1.0;
			}
		}
	}
	return;
}

//****************************************************************************
//
//  Purpose:
// 
//    CDFCHI evaluates the CDF of the chi square distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the chi square distribution 
//    given the others.
//
//    The value P of the cumulative distribution function is calculated 
//    directly.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The CDF of the chi square distribution can be evaluated 
//    within Mathematica by commands such as:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF [ ChiSquareDistribution [ DF ], X ]
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.4.19.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input, int %WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from X and DF;
//    2: Calculate X from P, Q and DF;
//    3: Calculate DF from P, Q and X.
//
//    Input/output, double %P, the integral from 0 to X of 
//    the chi-square distribution.  If this is an input value, it should
//    lie in the range [0,1].
//
//    Input/output, double %Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double %X, the upper limit of integration 
//    of the chi-square distribution.  If this is an input 
//    value, it should lie in the range: [0, +infinity).  If it is an output
//    value, it will be searched for in the range: [0,1.0D+300].
//
//    Input/output, double %DF, the degrees of freedom of the
//    chi-square distribution.  If this is an input value, it should lie
//    in the range: (0, +infinity).  If it is an output value, it will be
//    searched for in the range: [ 1.0D-300, 1.0D+300].
//
//    Output, int %STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +10, an error was returned from CUMGAM.
//
//    Output, double %BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
//****************************************************************************
void StatClass::cdfchi ( int %which, double %p, double %q, double %x, double %df,
						int %status, double %bound )
{

	double K2 = 0.0;
	double K4 = 0.5e0;
	double K5 = 5.0e0;
	double fx,cum,ccum,pq,porq = 0.0;
	unsigned long qhi,qleft,qporq =0;
	double T3,T6,T7,T8,T9,T10,T11;

	status = 0;
	bound = 0.0;
	//
	//     Check arguments
	//
	if((which < 1 || which > 3))
	{
		if(which < 1)
			bound = 1.0;
		else
			bound = 3.0e0;
		status = -1;
		return;
	}

	if (which != 1)
	{
		//
		//     P
		//
		if((p < 0.0) || (p > 1.0))
		{
			if(p < 0.0)
				bound = 0.0;
			else
				bound = 1.0;
			status = -2;
			return;
		}

		//
		//     Q
		//
		if ((q <= 0.0) || (q > 1.0))
		{
			if (q <= 0.0)
				bound = 0.0;
			else
				bound = 1.0;
			status = -3;
			return;
		}
	}

	if (which != 2)
	{
		//
		//     X
		//
		if (x < 0.0)
		{
			bound = 0.0;
			status = -4;
			return;
		}
	}

	if (which != 3)
	{
		//
		//     DF
		//
		if (df <= 0.0)
		{
			bound = 0.0;
			status = -5;
			return;
		}
	}


	if (which != 1)
	{
		//
		//     P + Q
		//
		pq = p+q;
		if (Math::Abs(pq - 1.0) > 3.0 * 2.22044604925031E-16 )
		{
			if (pq < 0.0)
				bound = 0.0;
			else
				bound = 1.0;

			status = 3;
			return;
		}

		//
		//     Select the minimum of P or Q
		//
		qporq = p <= q;
		if (qporq)
			porq = p;
		else
			porq = q;
	}

	//
	//     Calculate ANSWERS
	//
	if(1 == which) {
		//
		//     Calculating P and Q
		//
		status = 0;
		cumchi(x,df,p,q);
		if(porq > 1.5e0) {
			status = 10;
			return;
		}
	}
	else if (2 == which)
	{
		//
		//     Calculating X
		//
		x = 5.0e0;
		T3 = inf;
		T6 = atol;
		T7 = tol;
		dstinv(K2,T3,K4,K4,K5,T6,T7);
		status = 0;
		dinvr(status,x,fx,qleft,qhi);

		while (1 == status)
		{
			cumchi(x,df,cum,ccum);
			if (qporq)
				fx = cum-p;
			else
				fx = ccum-q;

			if (fx+porq > 1.5e0)
			{
				status = 10;
				return;
			}
			else
				dinvr(status,x,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = 0.0;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}
	else if(3 == which)
	{
		//
		//  Calculating DF
		//
		df = 5.0e0;
		T8 = zero;
		T9 = inf;
		T10 = atol;
		T11 = tol;
		dstinv(T8,T9,K4,K4,K5,T10,T11);
		status = 0;
		dinvr(status,df,fx,qleft,qhi);

		while (1 == status)
		{
			cumchi(x,df,cum,ccum);
			if (qporq)
				fx = cum-p;
			else
				fx = ccum-q;

			if (fx+porq > 1.5e0)
			{
				status = 10;
				return;
			}
			dinvr(status,df,fx,qleft,qhi);
		}

		if (-1 == status)
		{
			if (qleft)
			{
				status = 1;
				bound = zero;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}
	return;
}

//*****************************************************************************
//
//  Purpose:
// 
//    CDFCHN evaluates the CDF of the Noncentral Chi-Square.
//
//  Discussion:
//
//    This routine calculates any one parameter of the noncentral chi-square
//    distribution given values for the others.
//
//    The value P of the cumulative distribution function is calculated 
//    directly.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The computation time required for this routine is proportional
//    to the noncentrality parameter (PNONC).  Very large values of
//    this parameter can consume immense computer resources.  This is
//    why the search range is bounded by 10,000.
//
//    The CDF of the noncentral chi square distribution can be evaluated 
//    within Mathematica by commands such as:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.5.25.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input, int %WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from X, DF and PNONC;
//    2: Calculate X from P, DF and PNONC;
//    3: Calculate DF from P, X and PNONC;
//    4: Calculate PNONC from P, X and DF.
//
//    Input/output, double %P, the integral from 0 to X of 
//    the noncentral chi-square distribution.  If this is an input
//    value, it should lie in the range: [0, 1.0-1.0D-16).
//
//    Input/output, double %Q, is generally not used by this 
//    subroutine and is only included for similarity with other routines.
//    However, if P is to be computed, then a value will also be computed
//    for Q.
//
//    Input, double %X, the upper limit of integration of the 
//    noncentral chi-square distribution.  If this is an input value, it
//    should lie in the range: [0, +infinity).  If it is an output value,
//    it will be sought in the range: [0,1.0D+300].
//
//    Input/output, double %DF, the number of degrees of freedom 
//    of the noncentral chi-square distribution.  If this is an input value,
//    it should lie in the range: (0, +infinity).  If it is an output value,
//    it will be searched for in the range: [ 1.0D-300, 1.0D+300].
//
//    Input/output, double %PNONC, the noncentrality parameter of 
//    the noncentral chi-square distribution.  If this is an input value, it
//    should lie in the range: [0, +infinity).  If it is an output value,
//    it will be searched for in the range: [0,1.0D+4]
//
//    Output, int %STATUS, reports on the calculation.
//    0, if calculation completed correctly;
//    -I, if input parameter number I is out of range;
//    1, if the answer appears to be lower than the lowest search bound;
//    2, if the answer appears to be higher than the greatest search bound.
//
//    Output, double %BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//*****************************************************************************
void StatClass::cdfchn ( int %which, double %p, double %q, double %x, double %df,
						double %pnonc, int %status, double %bound )
{
	double K1 = 0.0;
	double K3 = 0.5;
	double K4 = 5.0;
	double fx,cum,ccum;
	unsigned long qhi,qleft;
	double T2,T5,T6,T7,T8,T9,T10,T11,T12,T13;

	status = 0;
	bound = 0.0;

	//
	//     Check arguments
	//
	if (which < 1)
	{
		bound = 1.0;
		status = -1;
		return;
	}
	else if (which > 4)
	{
		bound = 4.0;
		status = -1;
		return;
	}

	if (which != 1)
	{
		//
		//     P
		//
		if (p < 0.0)
		{
			bound = 0.0;
			status = -2;
			return;
		}
		else if (p > one)
		{
			bound = one;
			status = -2;
			return;
		}
	}

	if (which != 2)
	{
		//
		//     X
		//
		if (x < 0.0)
		{
			bound = 0.0;
			status = -4;
			return;
		}
	}

	if (which != 3)
	{
		//
		//     DF
		//
		if (df <= 0.0)
		{
			bound = 0.0;
			status = -5;
			return;
		}
	}

	if (which != 4)
	{
		//
		//     PNONC
		//
		if (pnonc < 0.0)
		{
			bound = 0.0;
			status = -6;
			return;
		}
	}

	//
	//     Calculate ANSWERS
	//
	if (1 == which)
	{
		//
		//     Calculating P and Q
		//
		cumchn(x,df,pnonc,p,q);
		status = 0;
	}
	else if (2 == which)
	{
		//
		//     Calculating X
		//
		x = 5.0;
		T2 = inf;
		T5 = atol;
		T6 = tol;
		dstinv(K1,T2,K3,K3,K4,T5,T6);
		status = 0;
		dinvr(status,x,fx,qleft,qhi);

		while (status == 1)
		{
			cumchn(x,df,pnonc,cum,ccum);
			fx = cum-p;
			dinvr(status,x,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = 0.0;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}
	else if(3 == which)
	{
		//
		//     Calculating DF
		//
		df = 5.0;
		T7 = zero;
		T8 = inf;
		T9 = atol;
		T10 = tol;
		dstinv(T7,T8,K3,K3,K4,T9,T10);
		status = 0;
		dinvr(status,df,fx,qleft,qhi);

		while (status == 1)
		{
			cumchn(x,df,pnonc,cum,ccum);
			fx = cum-p;
			dinvr(status,df,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = zero;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}
	else if(4 == which)
	{
		//
		//     Calculating PNONC
		//
		pnonc = 5.0;
		T11 = tent4;
		T12 = atol;
		T13 = tol;
		dstinv(K1,T11,K3,K3,K4,T12,T13);
		status = 0;
		dinvr(status,pnonc,fx,qleft,qhi);

		while (status == 1)
		{
			cumchn(x,df,pnonc,cum,ccum);
			fx = cum-p;
			dinvr(status,pnonc,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = zero;
			}
			else
			{
				status = 2;
				bound = tent4;
			}
		}
	}

	return;
}

//*****************************************************************************
//
//  Purpose:
// 
//    CDFF evaluates the CDF of the F distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the F distribution 
//    given the others.
//
//    The value P of the cumulative distribution function is calculated 
//    directly.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The value of the cumulative F distribution is not necessarily
//    monotone in either degree of freedom.  There thus may be two
//    values that provide a given CDF value.  This routine assumes
//    monotonicity and will find an arbitrary one of the two values.
//
//  Modified:
//
//    14 April 2007
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.6.2.
//
//  Parameters:
//
//    Input, int %WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from F, DFN and DFD;
//    2: Calculate F from P, Q, DFN and DFD;
//    3: Calculate DFN from P, Q, F and DFD;
//    4: Calculate DFD from P, Q, F and DFN.
//
//    Input/output, double %P, the integral from 0 to F of 
//    the F-density.  If it is an input value, it should lie in the
//    range [0,1].
//
//    Input/output, double %Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double %F, the upper limit of integration 
//    of the F-density.  If this is an input value, it should lie in the
//    range [0, +infinity).  If it is an output value, it will be searched
//    for in the range [0,1.0D+300].
//
//    Input/output, double %DFN, the number of degrees of 
//    freedom of the numerator sum of squares.  If this is an input value,
//    it should lie in the range: (0, +infinity).  If it is an output value,
//    it will be searched for in the range: [ 1.0D-300, 1.0D+300].
//
//    Input/output, double %DFD, the number of degrees of freedom 
//    of the denominator sum of squares.  If this is an input value, it should
//    lie in the range: (0, +infinity).  If it is an output value, it will
//    be searched for in the  range: [ 1.0D-300, 1.0D+300].
//
//    Output, int %STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1.
//
//    Output, double %BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
//*****************************************************************************

void StatClass::cdff ( int %which, double %p, double %q, double %f, double %dfn,
					  double %dfd, int %status, double %bound )

{
	double K2 = 0.0;
	double K4 = 0.5e0;
	double K5 = 5.0e0;
	double pq,fx,cum,ccum;
	unsigned long qhi,qleft,qporq = 0;
	double T3,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15;

	status = 0;
	bound = 0.0;

	//
	//  Check arguments
	//
	if (which < 1)
	{
		bound = 1.0;
		status = -1;
		return;
	}
	else if (which > 4)
	{
		bound = 4.0;
		status = -1;
		return;
	}

	if (which != 1)
	{
		//
		//     P
		//
		if (p < 0.0)
		{
			bound = 0.0;
			status = -2;
			return;
		}
		else if (p > 1.0)
		{
			bound = 1.0;
			status = -2;
			return;
		}

		//
		//     Q
		//
		if (q <= 0.0)
		{
			bound = 0.0;
			status = -3;
			return;
		}
		else if ( q > 1.0)
		{
			bound = 1.0;
			status = -3;
			return;
		}
	}

	if (which != 2)
	{
		//
		//     F
		//
		if (f < 0.0)
		{
			bound = 0.0;
			status = -4;
			return;
		}
	}

	if (which != 3)
	{
		//
		//     DFN
		//
		if (dfn <= 0.0)
		{
			bound = 0.0;
			status = -5;
			return;
		}
	}

	if (which != 4)
	{
		//
		//     DFD
		//
		if (dfd <= 0.0)
		{
			bound = 0.0;
			status = -6;
			return;
		}
	}

	if (which == 1)
	{
		//
		//     P + Q
		//
		pq = p+q;

		if (Math::Abs(pq - 1.0) < 3.0 * 2.22044604925031E-16 )
		{
			if (pq < 0.0)
				bound = 0.0;
			else
				bound = 1.0;

			status = 3;
			return;
		}
	}

	if (which != 1)
		qporq = p <= q;
	//
	//     Select the minimum of P or Q
	//     Calculate ANSWERS
	//

	if(1 == which)
	{
		//
		//     Calculating P
		//
		cumf(f,dfn,dfd,p,q);
		status = 0;
	}
	else if(2 == which)
	{
		//
		//     Calculating F
		//
		f = 5.0e0;
		T3 = inf;
		T6 = atol;
		T7 = tol;
		dstinv(K2,T3,K4,K4,K5,T6,T7);
		status = 0;
		dinvr(status,f,fx,qleft,qhi);

		while (status == 1)
		{
			cumf(f,dfn,dfd,cum,ccum);
			if (qporq)
				fx = cum-p;
			else
				fx = ccum-q;

			dinvr(status,f,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = 0.0;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}
	//
	//  Calculate DFN.
	//
	//  Note that, in the original calculation, the lower bound for DFN was 0.
	//  Using DFN = 0 causes an error in CUMF when it calls BETA_INC.
	//  The lower bound was set to the more reasonable value of 1.
	//  JVB, 14 April 2007.
	//
	else if ( 3 == which ) 
	{

		T8 = 1.0;
		T9 = inf;
		T10 = atol;
		T11 = tol;
		dstinv ( T8, T9, K4, K4, K5, T10, T11 );

		status = 0;
		dfn = 5.0;
		fx = 0.0;

		dinvr ( status, dfn, fx, qleft, qhi );

		while ( status == 1 )
		{
			cumf ( f, dfn, dfd, cum, ccum );

			if ( p <= q ) 
			{
				fx = cum - p;
			}
			else
			{
				fx = ccum - q;
			}
			dinvr ( status, dfn, fx, qleft, qhi );
		}

		if ( status == -1 )
		{
			if ( qleft )
			{
				status = 1;
				bound = 1.0;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}
	//
	//  Calculate DFD.
	//
	//  Note that, in the original calculation, the lower bound for DFD was 0.
	//  Using DFD = 0 causes an error in CUMF when it calls BETA_INC.
	//  The lower bound was set to the more reasonable value of 1.
	//  JVB, 14 April 2007.
	//
	//
	else if ( 4 == which )
	{

		T12 = 1.0;
		T13 = inf;
		T14 = atol;
		T15 = tol;
		dstinv ( T12, T13, K4, K4, K5, T14, T15 );

		status = 0;
		dfd = 5.0;
		fx = 0.0;
		dinvr ( status, dfd, fx, qleft, qhi );

		while ( status == 1 )
		{
			cumf ( f, dfn, dfd, cum, ccum );

			if ( p <= q ) 
			{
				fx = cum - p;
			}
			else
			{
				fx = ccum - q;
			}
			dinvr ( status, dfd, fx, qleft, qhi );
		}

		if ( status == -1 )
		{
			if ( qleft )
			{
				status = 1;
				bound = 1.0;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}

	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    CDFFNC evaluates the CDF of the Noncentral F distribution.
//
//  Discussion:
//
//    This routine originally used 1.0E+300 as the upper bound for the 
//    interval in which many of the missing parameters are to be sought.
//    Since the underlying rootfinder routine needs to evaluate the 
//    function at this point, it is no surprise that the program was
//    experiencing overflows.  A less extravagant upper bound
//    is being tried for now!
//
//
//    This routine calculates any one parameter of the Noncentral F distribution 
//    given the others.
//
//    The value P of the cumulative distribution function is calculated 
//    directly.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The computation time required for this routine is proportional
//    to the noncentrality parameter PNONC.  Very large values of
//    this parameter can consume immense computer resources.  This is
//    why the search range is bounded by 10,000.
//
//    The value of the cumulative noncentral F distribution is not
//    necessarily monotone in either degree of freedom.  There thus
//    may be two values that provide a given CDF value.  This routine
//    assumes monotonicity and will find an arbitrary one of the two
//    values.
//
//    The CDF of the noncentral F distribution can be evaluated 
//    within Mathematica by commands such as:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF [ NoncentralFRatioDistribution [ DFN, DFD, PNONC ], X ]
//
//  Modified:
//
//    15 June 2004
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.6.20.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input, int %WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from F, DFN, DFD and PNONC;
//    2: Calculate F from P, Q, DFN, DFD and PNONC;
//    3: Calculate DFN from P, Q, F, DFD and PNONC;
//    4: Calculate DFD from P, Q, F, DFN and PNONC;
//    5: Calculate PNONC from P, Q, F, DFN and DFD.
//
//    Input/output, double %P, the integral from 0 to F of 
//    the noncentral F-density.  If P is an input value it should
//    lie in the range [0,1) (Not including 1!).
//
//    Dummy, double %Q, is not used by this subroutine,
//    and is only included for similarity with the other routines.
//    Its input value is not checked.  If P is to be computed, the
//    Q is set to 1 - P.
//
//    Input/output, double %F, the upper limit of integration 
//    of the noncentral F-density.  If this is an input value, it should
//    lie in the range: [0, +infinity).  If it is an output value, it
//    will be searched for in the range: [0,1.0D+30].
//
//    Input/output, double %DFN, the number of degrees of freedom 
//    of the numerator sum of squares.  If this is an input value, it should
//    lie in the range: (0, +infinity).  If it is an output value, it will
//    be searched for in the range: [ 1.0, 1.0D+30].
//
//    Input/output, double %DFD, the number of degrees of freedom 
//    of the denominator sum of squares.  If this is an input value, it should
//    be in range: (0, +infinity).  If it is an output value, it will be
//    searched for in the range [1.0, 1.0D+30].
//
//    Input/output, double %PNONC, the noncentrality parameter
//    If this is an input value, it should be nonnegative.
//    If it is an output value, it will be searched for in the range: [0,1.0D+4].
//
//    Output, int %STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1.
//
//    Output, double %BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
//****************************************************************************80

void StatClass::cdffnc ( int %which, double %p, double %q, double %f, double %dfn,
						double %dfd, double %phonc, int %status, double %bound )

{
	double K1 = 0.0;
	double K3 = 0.5e0;
	double K4 = 5.0e0;
	double fx,cum,ccum;
	unsigned long qhi,qleft;
	double T2,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17;

	status = 0;
	bound = 0.0;
	//
	//     Check arguments
	//
	if (which < 1)
		bound = 1.0;
	else if (which > 5)
		bound = 5.0;

	if (bound != 0.0)
	{
		status = -1;
		return;
	}

	if (which != 1)
	{
		//
		//     P
		//
		if (p < 0.0)
		{
			bound = 0.0;
			status = -2;
			return;
		}
		else if (p > one)
		{
			bound = one;
			status = -2;
			return;
		}
	}

	if (which != 2)
	{
		//
		//     F
		//
		if (f < 0.0)
		{
			bound = 0.0;
			status = -4;
			return;
		}
	}

	if (which != 3)
	{
		//
		//     DFN
		//
		if (dfn <= 0.0)
		{
			bound = 0.0;
			status = -5;
			return;
		}
	}

	if (which != 4)
	{
		//
		//     DFD
		//
		if (dfd <= 0.0)
		{
			bound = 0.0;
			status = -6;
			return;
		}
	}

	if (which != 5)
	{
		//
		//     PHONC
		//
		if (phonc < 0.0)
		{
			bound = 0.0;
			status = -7;
			return;
		}
	}

	//
	//     Calculate ANSWERS
	//
	if(1 == which)
	{
		//
		//     Calculating P
		//
		cumfnc(f,dfn,dfd,phonc,p,q);
		status = 0;
	}
	else if(2 == which)
	{
		//
		//     Calculating F
		//
		f = 5.0e0;
		T2 = inf;
		T5 = atol;
		T6 = tol;
		dstinv(K1,T2,K3,K3,K4,T5,T6);
		status = 0;
		dinvr(status,f,fx,qleft,qhi);

		while (status == 1)
		{
			cumfnc(f,dfn,dfd,phonc,cum,ccum);
			fx = cum-p;
			dinvr(status,f,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = 0.0;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}
	else if(3 == which)
	{
		//
		//     Calculating DFN
		//
		dfn = 5.0e0;
		T7 = zero;
		T8 = inf;
		T9 = atol;
		T10 = tol;
		dstinv(T7,T8,K3,K3,K4,T9,T10);
		status = 0;
		dinvr(status,dfn,fx,qleft,qhi);

		while (status == 1)
		{
			cumfnc(f,dfn,dfd,phonc,cum,ccum);
			fx = cum-p;
			dinvr(status,dfn,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = zero;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}
	else if (4 == which)
	{
		//
		//     Calculating DFD
		//
		dfd = 5.0e0;
		T11 = zero;
		T12 = inf;
		T13 = atol;
		T14 = tol;
		dstinv(T11,T12,K3,K3,K4,T13,T14);
		status = 0;
		dinvr(status,dfd,fx,qleft,qhi);

		while (status == 1)
		{
			cumfnc(f,dfn,dfd,phonc,cum,ccum);
			fx = cum-p;
			dinvr(status,dfd,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = zero;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}
	else if(5 == which)
	{
		//
		//     Calculating PHONC
		//
		phonc = 5.0e0;
		T15 = tent4;
		T16 = atol;
		T17 = tol;
		dstinv(K1,T15,K3,K3,K4,T16,T17);
		status = 0;
		dinvr(status,phonc,fx,qleft,qhi);

		while (status == 1)
		{
			cumfnc(f,dfn,dfd,phonc,cum,ccum);
			fx = cum-p;
			dinvr(status,phonc,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = 0.0;
			}
			else
			{
				status = 2;
				bound = tent4;
			}
		}
	}

	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    CDFGAM evaluates the CDF of the Gamma Distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the Gamma distribution 
//    given the others.
//
//    The cumulative distribution function P is calculated directly.
//
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The gamma density is proportional to T**(SHAPE - 1) * EXP(- SCALE * T)
//
//  Reference:
//
//    Armido DiDinato and Alfred Morris, 
//    Computation of the incomplete gamma function ratios and their inverse,
//    ACM Transactions on Mathematical Software,
//    Volume 12, 1986, pages 377-393.
//
//  Parameters:
//
//    Input, int %WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from X, SHAPE and SCALE;
//    2: Calculate X from P, Q, SHAPE and SCALE;
//    3: Calculate SHAPE from P, Q, X and SCALE;
//    4: Calculate SCALE from P, Q, X and SHAPE.
//
//    Input/output, double %P, the integral from 0 to X of the 
//    Gamma density.  If this is an input value, it should lie in the 
//    range: [0,1].
//
//    Input/output, double %Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double %X, the upper limit of integration of 
//    the Gamma density.  If this is an input value, it should lie in the
//    range: [0, +infinity).  If it is an output value, it will lie in
//    the range: [0,1E300].
//
//    Input/output, double %SHAPE, the shape parameter of the 
//    Gamma density.  If this is an input value, it should lie in the range: 
//    (0, +infinity).  If it is an output value, it will be searched for
//    in the range: [1.0D-300,1.0D+300].
//
//    Input/output, double %SCALE, the scale parameter of the 
//    Gamma density.  If this is an input value, it should lie in the range
//    (0, +infinity).  If it is an output value, it will be searched for
//    in the range: (1.0D-300,1.0D+300].
//
//    Output, int %STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +10, if the Gamma or inverse Gamma routine cannot compute the answer.  
//    This usually happens only for X and SHAPE very large (more than 1.0D+10.
//
//    Output, double %BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
//****************************************************************************80
void StatClass::cdfgam ( int %which, double %p, double %q, double %x, double %shape,
						double %scale, int %status, double %bound )
{

	double K5 = 0.5e0;
	double K6 = 5.0e0;
	double xx,fx,xscale,cum,ccum,pq,porq = 0;
	int ierr;
	unsigned long qhi,qleft,qporq = 0;
	double T2,T3,T4,T7,T8,T9;

	status = 0;
	bound = 0.0;
	//
	//     Check arguments
	//
	if(which < 1)
	{
		bound = 1.0;
		status = -1;
		return;
	}
	else if (which > 4)
	{
		bound = 4.0e0;
		status = -1;
		return;
	}

	if(which != 1)
	{
		//
		//     P
		//
		if(p < 0.0)
		{
			bound = 0.0;
			status = -2;
			return;
		}
		else if (p > 1.0)
		{
			bound = 1.0;
			status = -2;
			return;
		}
		//
		//     Q
		//
		if(q < 0.0)
		{
			bound = 0.0;
			status = -3;
			return;
		}
		else if (q > 1.0)
		{
			bound = 1.0;
			status = -3;
			return;
		}
		//
		//     P + Q
		//
		pq = p + q;
		if (Math::Abs(pq - 1.0) > 3.0 * 2.22044604925031E-16)
		{
			if (pq < 0.0)
			{
				bound = 0.0;
				status = 3;
				return;
			}
			else if (pq > 1.0)
			{
				bound = 1.0;
				status = 3;
				return;
			}
		}
		//
		//     Select the minimum of P or Q
		//
		qporq = p <= q;
		porq = (qporq) ? p : q;
	}

	if (which != 2)
	{
		//
		//     X
		//
		if(x < 0.0)
		{
			bound = 0.0;
			status = -4;
			return;
		}
	}

	if (which != 3)
	{
		//
		//     SHAPE
		//
		if (shape <= 0.0)
		{
			bound = 0.0;
			status = -5;
			return;
		}
	}
	if (which != 4)
	{
		//
		//     SCALE
		//
		if (scale <= 0.0)
		{
			bound = 0.0;
			status = -6;
			return;
		}
	}

	//
	//     Calculate ANSWERS
	//
	if(1 == which)
	{
		//
		//     Calculating P
		//
		status = 0;
		xscale = x*scale;
		cumgam(xscale,shape,p,q);
		if (porq > 1.5e0)
			status = 10;
	}
	else if(2 == which)
	{
		//
		//     Computing X
		//
		T2 = -1.0;
		gamma_inc_inv ( shape, xx, T2, p, q, ierr );
		if(ierr < 0.0) {
			status = 10;
			return;
		}
		else  {
			x = xx/ scale;
			status = 0;
		}
	}
	else if(3 == which)
	{
		//
		//     Computing SHAPE
		//
		shape = 5.0e0;
		xscale = x*scale;
		T3 = zero;
		T4 = inf;
		T7 = atol;
		T8 = tol;
		dstinv(T3,T4,K5,K5,K6,T7,T8);
		status = 0;
		dinvr(status,shape,fx,qleft,qhi);

		while (1 == status)
		{
			cumgam(xscale,shape,cum,ccum);
			if (qporq)
				fx = cum-p;
			else
				fx = ccum-q;

			if((qporq && cum > 1.5e0) || (!qporq && ccum > 1.5e0))
			{
				status = 10;
				return;
			}

			dinvr(status,shape,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = zero;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}
	else if(4 == which)
	{
		//
		//     Computing SCALE
		//
		T9 = -1.0;
		gamma_inc_inv ( shape, xx, T9, p, q, ierr );
		if(ierr < 0.0)
		{
			status = 10;
			return;
		}
		else
		{
			scale = xx/ x;
			status = 0;
		}
	}

	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    CDFNBN evaluates the CDF of the Negative Binomial distribution
//
//  Discussion:
//
//    This routine calculates any one parameter of the negative binomial
//    distribution given values for the others.
//
//    The cumulative negative binomial distribution returns the
//    probability that there will be F or fewer failures before the
//    S-th success in binomial trials each of which has probability of
//    success PR.
//
//    The individual term of the negative binomial is the probability of
//    F failures before S successes and is
//    Choose( F, S+F-1 ) * PR%(S) * (1-PR)%F
//
//    Computation of other parameters involve a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.5.26.
//
//  Parameters:
//
//    Input, int WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from F, S, PR and OMPR;
//    2: Calculate F from P, Q, S, PR and OMPR;
//    3: Calculate S from P, Q, F, PR and OMPR;
//    4: Calculate PR and OMPR from P, Q, F and S.
//
//    Input/output, double P, the cumulation from 0 to F of 
//    the negative binomial distribution.  If P is an input value, it
//    should lie in the range [0,1].
//
//    Input/output, double Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double F, the upper limit of cumulation of 
//    the binomial distribution.  There are F or fewer failures before 
//    the S-th success.  If this is an input value, it may lie in the
//    range [0,+infinity), and if it is an output value, it will be searched
//    for in the range [0,1.0D+300].
//
//    Input/output, double S, the number of successes.
//    If this is an input value, it should lie in the range: [0, +infinity).
//    If it is an output value, it will be searched for in the range: 
//    [0, 1.0D+300].
//
//    Input/output, double PR, the probability of success in each 
//    binomial trial.  Whether an input or output value, it should lie in the
//    range [0,1].
//
//    Input/output, double OMPR, the value of (1-PR).  Whether an
//    input or output value, it should lie in the range [0,1].  
//
//    Output, int STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +4, if PR + OMPR /= 1.
//
//    Output, double BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
//****************************************************************************80
void StatClass::cdfnbn ( int %which, double %p, double %q, double %s, double %xn,
						double %pr, double %ompr, int %status, double %bound )
{
	double K2 = 0.0;
	double K4 = 0.5e0;
	double K5 = 5.0e0;
	double K11 = 1.0;
	double fx,xhi,xlo,pq,prompr,cum,ccum;
	unsigned long qhi,qleft,qporq = 0;
	double T3,T6,T7,T8,T9,T10,T12,T13;

	status = 0;
	bound = 0.0;
	//
	//     Check arguments
	//
	// Which is between 1 and 4, otherwise exit
	if (which < 1)
	{
		bound = 1.0;
		status = -1;
		return;
	}
	else if (which > 4)
	{
		bound = 4.0e0;
		status = -1;
		return;
	}

	if (which != 1)
	{
		//
		//     P
		//
		if (p < 0.0)
		{
			bound = 0.0;
			status = -2;
			return;
		}
		else if (p > 1.0)
		{
			bound = 1.0;
			status = -2;
			return;
		}

		//
		//     Q
		//
		if (q <= 0.0)
		{
			bound = 0.0;
			status = -3;
			return;
		}
		else if ( q > 1.0)
		{
			bound = 1.0;
			status = -3;
			return;
		}
		//
		//     P + Q
		//
		pq = p+q;
		if (Math::Abs(pq - 1.0) > 3.0 * 2.22044604925031E-16)
		{
			if(pq < 0.0)
			{
				bound = 0.0;
				status = 3;
				return;
			}
			else if (pq > 1.0)
			{
				bound = 1.0;
				status = 3;
				return;
			}
		}

		qporq = p <= q;

	}

	if (which != 2)
		//
		//     S
		//
		if (s < 0.0)
		{
			bound = 0.0;
			status = -4;
			return;
		}

		if (which != 3)
			//
			//     XN
			//
			if (xn < 0.0)
			{
				bound = 0.0;
				status = -5;
				return;
			}

			if (which != 4)
			{
				//
				//     PR
				//
				if (pr < 0.0)
				{
					bound = 0.0;
					status = -6;
					return;
				}
				else if (pr > 1.0)
				{
					bound = 1.0;
					status = -6;
					return;
				}

				//
				//     OMPR
				//
				if (ompr < 0.0)
				{
					bound = 0.0;
					status = -7;
					return;
				}
				else if (ompr > 1.0)
				{
					bound = 1.0;
					status = -7;
					return;
				}

				//
				//     PR + OMPR
				//
				prompr = pr+ompr;
				if (Math::Abs(prompr - 1.0) > 3.0 * 2.22044604925031E-16)
				{
					if (prompr < 0.0)
					{
						bound = 0.0;
						status = 4;
						return;
					}
					else if (prompr > 1.0)
					{
						bound = 1.0;
						status = 4;
						return;
					}
				}
			}

			//
			//     Select the minimum of P or Q
			//     Calculate ANSWERS
			//
			if (1 == which)
			{
				//
				//     Calculating P
				//
				cumnbn(s,xn,pr,ompr,p,q);
				status = 0;
			}
			else if (2 == which)
			{
				//
				//     Calculating S
				//
				s = 5.0e0;
				T3 = inf;
				T6 = atol;
				T7 = tol;
				dstinv(K2,T3,K4,K4,K5,T6,T7);
				status = 0;
				dinvr(status,s,fx,qleft,qhi);

				while (status == 1)
				{
					cumnbn(s,xn,pr,ompr,cum,ccum);
					if (qporq)
						fx = cum-p;
					else 
						fx = ccum-q;
					dinvr(status,s,fx,qleft,qhi);
				}

				if (status == -1)
				{
					if (qleft)
					{
						status = 1;
						bound = 0.0;
					}
					else
					{
						status = 2;
						bound = inf;
					}
				}
			}
			else if (3 == which)
			{
				//
				//     Calculating XN
				//
				xn = 5.0e0;
				T8 = inf;
				T9 = atol;
				T10 = tol;
				dstinv(K2,T8,K4,K4,K5,T9,T10);
				status = 0;
				dinvr(status,xn,fx,qleft,qhi);

				while (1 == status)
				{
					cumnbn(s,xn,pr,ompr,cum,ccum);
					if (qporq)
						fx = cum-p;
					else
						fx = ccum-q;
					dinvr(status,xn,fx,qleft,qhi);
				}

				if (status == -1)
				{
					if (qleft)
					{
						status = 1;
						bound = 0.0;
					}
					else
					{
						status = 2;
						bound = inf;
					}
				}
			}
			else if (4 == which)
			{
				//
				//     Calculating PR and OMPR
				//
				T12 = atol;
				T13 = tol;
				dstzr(K2,K11,T12,T13);
				if (qporq)
				{
					status = 0;
					dzror(status,pr,fx,xlo,xhi,qleft,qhi);
					ompr = one-pr;

					while (1 == status)
					{
						cumnbn(s,xn,pr,ompr,cum,ccum);
						fx = cum-p;
						dzror(status,pr,fx,xlo,xhi,qleft,qhi);
						ompr = one-pr;
					}
				}
				else
				{
					status = 0;
					dzror(status,ompr,fx,xlo,xhi,qleft,qhi);
					pr = one-ompr;
					while (1 == status)
					{
						cumnbn(s,xn,pr,ompr,cum,ccum);
						fx = ccum-q;
						dzror(status,ompr,fx,xlo,xhi,qleft,qhi);
						pr = one-ompr;
					}
				}

				if (status == -1)
				{
					if (qleft)
					{
						status = 1;
						bound = 0.0;
					}
					else
					{
						status = 2;
						bound = 1.0;
					}
				}
			}
			return;
}

//*****************************************************************************
//
//  Purpose:
// 
//    CDFNOR evaluates the CDF of the Normal distribution.
//
//  Discussion:
//
//    A slightly modified version of ANORM from SPECFUN
//    is used to calculate the cumulative standard normal distribution.
//
//    The rational functions from pages 90-95 of Kennedy and Gentle
//    are used as starting values to Newton's Iterations which 
//    compute the inverse standard normal.  Therefore no searches are
//    necessary for any parameter.
//
//    For X < -15, the asymptotic expansion for the normal is used  as
//    the starting value in finding the inverse standard normal.
//
//    The normal density is proportional to
//    Math::Exp( - 0.5D+00 * (( X - MEAN)/SD)**2)
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.2.12.
//
//    William Cody,
//    Algorithm 715: SPECFUN - A Portable FORTRAN Package of
//      Special Function Routines and Test Drivers,
//    ACM Transactions on Mathematical Software,
//    Volume 19, pages 22-32, 1993.
//
//    Kennedy and Gentle,
//    Statistical Computing,
//    Marcel Dekker, NY, 1980,
//    QA276.4  K46
//
//  Parameters:
//
//    Input, int %WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from X, MEAN and SD;
//    2: Calculate X from P, Q, MEAN and SD;
//    3: Calculate MEAN from P, Q, X and SD;
//    4: Calculate SD from P, Q, X and MEAN.
//
//    Input/output, double %P, the integral from -infinity to X 
//    of the Normal density.  If this is an input or output value, it will
//    lie in the range [0,1].
//
//    Input/output, double %Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double %X, the upper limit of integration of 
//    the Normal density.
//
//    Input/output, double %MEAN, the mean of the Normal density.
//
//    Input/output, double %SD, the standard deviation of the 
//    Normal density.  If this is an input value, it should lie in the
//    range (0,+infinity).
//
//    Output, int %STATUS, the status of the calculation.
//    0, if calculation completed correctly;
//    -I, if input parameter number I is out of range;
//    1, if answer appears to be lower than lowest search bound;
//    2, if answer appears to be higher than greatest search bound;
//    3, if P + Q /= 1.
//
//    Output, double %BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
//*****************************************************************************

void StatClass::cdfnor ( int %which, double %p, double %q, double %x, double %mean,	double %sd, int %status, double %bound )
{
	double z,pq;

	status = 0;
	bound = 0.0;
	//
	//     Check arguments
	//

	if (which < 1)
	{
		bound = 1.0;
		status = -1;
		return;
	}
	else if (which > 4)
	{
		bound = 4.0e0;
		status = -1;
		return;
	}

	if (which != 1)
	{
		//
		//     P
		//
		if (p <= 0.0)
		{
			bound = 0.0;
			status = -2;
			return;
		}
		else if ( p > 1.0)
		{
			bound = 1.0;
			status = -2;
			return;
		}

		//
		//     Q
		//
		if (q <= 0.0)
		{
			bound = 0.0;
			status = -3;
			return;
		}
		else if (q > 1.0)
		{
			bound = 1.0;
			status = -3;
			return;
		}

		//
		//     P + Q
		//
		pq = p+q;
		if (Math::Abs(pq - 1.0) > 3.0 * 2.22044604925031E-16)
		{
			if (pq < 0.0)
			{
				bound = 0.0;
				status = 3;
				return;
			}
			else if (pq > 1)
			{
				bound = 1.0;
				status = 3;
				return;
			}
		}
	}

	if (which != 4)
	{
		//
		//     SD
		//
		if (sd <= 0.0)
		{
			bound = 0.0;
			status = -6;
			return;
		}
	}

	//
	//     Calculate ANSWERS
	//
	if(1 == which) {
		//
		//     Computing P
		//
		z = (x-mean)/ sd;
		cumnor(z,p,q);
	}
	else if(2 == which) {
		//
		//     Computing X
		//
		z = dinvnr(p,q);
		x = sd*z+mean;
	}
	else if(3 == which) {
		//
		//     Computing the MEAN
		//
		z = dinvnr(p,q);
		mean = x-sd*z;
	}
	else if(4 == which) {
		//
		//     Computing SD
		//
		z = dinvnr(p,q);
		sd = (x-mean)/z;
	}
	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    CDFPOI evaluates the CDF of the Poisson distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the Poisson distribution 
//    given the others.
//
//    The value P of the cumulative distribution function is calculated 
//    directly.
//
//    Computation of other parameters involve a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.4.21.
//
//  Parameters:
//
//    Input, int %WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from S and XLAM;
//    2: Calculate A from P, Q and XLAM;
//    3: Calculate XLAM from P, Q and S.
//
//    Input/output, double %P, the cumulation from 0 to S of the 
//    Poisson density.  Whether this is an input or output value, it will
//    lie in the range [0,1].
//
//    Input/output, double %Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double %S, the upper limit of cumulation of
//    the Poisson CDF.  If this is an input value, it should lie in
//    the range: [0, +infinity).  If it is an output value, it will be
//    searched for in the range: [0,1.0D+300].
//
//    Input/output, double %XLAM, the mean of the Poisson 
//    distribution.  If this is an input value, it should lie in the range
//    [0, +infinity).  If it is an output value, it will be searched for
//    in the range: [0,1E300].
//
//    Output, int %STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1.
//
//    Output, double %BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
//****************************************************************************80

void StatClass::cdfpoi ( int %which, double %p, double %q, double %s, double %xlam,	int %status, double %bound )
{

	double K2 = 0.0;
	double K4 = 0.5e0;
	double K5 = 5.0e0;
	double fx,cum,ccum,pq;
	unsigned long qhi,qleft,qporq = 0;
	double T3,T6,T7,T8,T9,T10;

	status = 0;
	bound = 0.0;
	//
	//     Check arguments
	//
	if (which < 1)
	{
		bound = 1.0;
		status = -1;
		return;
	}
	else if (which > 3)
	{
		bound = 3.0e0;
		status = -1;
		return;
	}

	if (which != 1)
	{
		//
		//     P
		//
		if (p < 0.0)
		{
			bound = 0.0;
			status = -2;
			return;
		}
		else if ( p > 1.0)
		{
			bound = 1.0;
			status = -2;
			return;
		}

		//
		//     Q
		//
		if (q < 0.0)
		{
			bound = 0.0;
			status = -3;
			return;
		}
		else if (q > 1.0)
		{
			bound = 1.0;
			status = -3;
			return;
		}

		pq = p+q;

		if (Math::Abs(pq - 1.0) <= 3.0 * 2.22044604925031E-16)
		{
			if (pq < 0.0)
			{
				bound = 0.0;
				status = 3;
				return;
			}
			else if (pq > 1.0)
			{
				bound = 1.0;
				status = 3;
				return;
			}
		}

		qporq = p <= q;
	}

	if (which != 2)
	{
		//
		//     S
		//
		if (s < 0.0)
		{
			bound = 0.0;
			status = -4;
			return;
		}
	}

	if (which != 3)
	{
		//
		//     XLAM
		//
		if (xlam < 0.0)
		{
			bound = 0.0;
			status = -5;
			return;
		}
	}

	//
	//     Select the minimum of P or Q
	//     Calculate ANSWERS
	//
	if (1 == which)
	{
		//
		//     Calculating P
		//
		cumpoi(s,xlam,p,q);
		status = 0;
	}
	else if(2 == which)
	{
		//
		//     Calculating S
		//
		s = 5.0e0;
		T3 = inf;
		T6 = atol;
		T7 = tol;
		dstinv(K2,T3,K4,K4,K5,T6,T7);
		status = 0;
		dinvr(status,s,fx,qleft,qhi);

		while (status == 1)
		{
			cumpoi(s,xlam,cum,ccum);
			if (qporq)
				fx = cum-p;
			else
				fx = ccum-q;

			dinvr(status,s,fx,qleft,qhi);

		}

		if  (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = 0.0;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}
	else if(3 == which)
	{
		//
		//     Calculating XLAM
		//
		xlam = 5.0e0;
		T8 = inf;
		T9 = atol;
		T10 = tol;
		dstinv(K2,T8,K4,K4,K5,T9,T10);
		status = 0;
		dinvr(status,xlam,fx,qleft,qhi);

		while (status == 1)
		{
			cumpoi(s,xlam,cum,ccum);
			if (qporq)
				fx = cum-p;
			else
				fx = ccum-q;

			dinvr(status,xlam,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = 0.0;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}

	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    CDFT evaluates the CDF of the T distribution.
//
//  Discussion:
//
//    This routine calculates any one parameter of the T distribution 
//    given the others.
//
//    The value P of the cumulative distribution function is calculated 
//    directly.
//
//    Computation of other parameters involve a seach for a value that
//    produces the desired value of P.   The search relies on the
//    monotonicity of P with respect to the other parameters.
//
//    The original version of this routine allowed the search interval
//    to extend from -1.0E+300 to +1.0E+300, which is fine until you 
//    try to evaluate a function at such a point!
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.5.27.
//
//  Parameters:
//
//    Input, int %WHICH, indicates which argument is to be calculated
//    from the others.
//    1 : Calculate P and Q from T and DF;
//    2 : Calculate T from P, Q and DF;
//    3 : Calculate DF from P, Q and T.
//
//    Input/output, double %P, the integral from -infinity to T of 
//    the T-density.  Whether an input or output value, this will lie in the
//    range [0,1].
//
//    Input/output, double %Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//
//    Input/output, double %T, the upper limit of integration of 
//    the T-density.  If this is an input value, it may have any value.
//    It it is an output value, it will be searched for in the range
//    [ -1.0D+30, 1.0D+30 ].
//
//    Input/output, double %DF, the number of degrees of freedom
//    of the T distribution.  If this is an input value, it should lie
//    in the range: (0 , +infinity).  If it is an output value, it will be
//    searched for in the range: [1, 1.0D+10].
//
//    Output, int %STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1.
//
//    Output, double %BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
//****************************************************************************80
void StatClass::cdft ( int %which, double %p, double %q, double %t, double %df,
					  int %status, double %bound )

{
	double K4 = 0.5e0;
	double K5 = 5.0e0;
	double fx,cum,ccum,pq;
	unsigned long qhi,qleft,qporq = 0;
	double T2,T3,T6,T7,T8,T9,T10,T11;

	status = 0;
	bound = 0.0;
	//
	//     Check arguments
	//
	if (which < 1)
	{
		bound = 1.0;
		status = -1;
		return;
	}

	if (which > 3)
	{
		bound = 3.0e0;
		status = -1;
		return;
	}

	if (which != 1)
	{
		//
		//     P
		//
		if (p <= 0.0)
		{
			bound = 0.0;
			status = -2;
			return;
		}

		if (p > 1.0)
		{
			bound = 1.0;
			status = -2;
			return;
		}

		//
		//     Q
		//
		if (q <= 0.0)
		{
			bound = 0.0;
			status = -3;
			return;
		}

		if (q > 1.0)
		{
			bound = 1.0;
			status = -3;
			return;
		}

		//
		//     P + Q
		//
		pq = p+q;

		if (Math::Abs(pq - 1.0) > 3.0 * 2.22044604925031E-16)
		{

			if (pq <= 0.0)
			{
				bound = 0.0;
				status = 3;
				return;
			}

			if (pq > 1.0)
			{
				bound = 1.0;
				status = 3;
				return;
			}
		}

		qporq = p <= q;

	}

	if (which != 3)
	{
		//
		//     DF
		//
		if (df <= 0.0)
		{
			bound = 0.0;
			status = -5;
			return;
		}
	}

	//
	//     Select the minimum of P or Q
	//     Calculate ANSWERS
	//
	if (1 == which)
	{
		//
		//     Computing P and Q
		//
		cumt(t,df,p,q);
		status = 0;
	}
	else if(2 == which)
	{
		//
		//     Computing T
		//     .. Get initial approximation for T
		//
		t = dt1(p,q,df);
		T2 = -inf;
		T3 = inf;
		T6 = atol;
		T7 = tol;
		dstinv(T2,T3,K4,K4,K5,T6,T7);
		status = 0;
		dinvr(status,t,fx,qleft,qhi);

		while (status == 1)
		{
			cumt(t,df,cum,ccum);
			if (qporq)
				fx = cum-p;
			else
				fx = ccum-q;
			dinvr(status,t,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if (qleft)
			{
				status = 1;
				bound = -inf;
			}
			else
			{
				status = 2;
				bound = inf;
			}
		}
	}
	else if(3 == which) 
	{
		//
		//     Computing DF
		//
		df = 5.0e0;
		T8 = zero;
		T9 = 1.0e10;
		T10 = atol;
		T11 = tol;
		dstinv(T8,T9,K4,K4,K5,T10,T11);
		status = 0;
		dinvr(status,df,fx,qleft,qhi);

		while (status == 1)
		{
			cumt(t,df,cum,ccum);

			if (qporq)
				fx = cum-p;
			else
				fx = ccum-q;

			dinvr(status,df,fx,qleft,qhi);
		}

		if (status == -1)
		{
			if(qleft)
			{
				status = 1;
				bound = zero;
			}
			else
			{
				status = 2;
				bound = 1.0e10;
			}
		}
	}
}

//****************************************************************************80**
//
//  Purpose:
//
//    CHI_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
//
//  Discussion:
//
//    The CDF of the noncentral chi square distribution can be evaluated
//    within Mathematica by commands such as:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF [ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
//
//  Modified:
//
//    12 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input/output, int %N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double %X, the argument of the function.
//
//    Output, double %LAMBDA, the noncentrality parameter.
//
//    Output, int %DF, the number of degrees of freedom.
//
//    Output, double %CDF, the noncentral chi CDF.
//
//****************************************************************************80**
void StatClass::chi_noncentral_cdf_values ( int %n_data, double %x, double %lambda, int %df, double %cdf )
{
	const int N_MAX = 27;

	array<double> ^cdf_vec = {
		0.839944, 0.695906, 0.535088, 
		0.764784, 0.620644, 0.469167,
		0.307088, 0.220382, 0.150025, 
		0.307116E-02, 0.176398E-02, 0.981679E-03, 
		0.165175E-01, 0.202342E-03, 0.498448E-06, 
		0.151325E-01, 0.209041E-02, 0.246502E-03, 
		0.263684E-01, 0.185798E-01, 0.130574E-01, 
		0.583804E-01, 0.424978E-01, 0.308214E-01, 
		0.105788, 0.794084E-01, 0.593201E-01 
	};
	array<int> ^df_vec = {
		1,   2,   3, 
		1,   2,   3, 
		1,   2,   3,
		1,   2,   3, 
		60,  80, 100, 
		1,   2,   3, 
		10,  10,  10, 
		10,  10,  10, 
		10,  10,  10 
	};
	array<double> ^lambda_vec = { 
		0.5,  0.5,  0.5, 
		1.0,  1.0,  1.0, 
		5.0,  5.0,  5.0, 
		20.0, 20.0, 20.0, 
		30.0, 30.0, 30.0, 
		5.0,  5.0,  5.0, 
		2.0,  3.0,  4.0, 
		2.0,  3.0,  4.0, 
		2.0,  3.0,  4.0 
	};
	array<double> ^x_vec = {
		3.000,  3.000,  3.000, 
		3.000,  3.000,  3.000, 
		3.000,  3.000,  3.000, 
		3.000,  3.000,  3.000, 
		60.000, 60.000, 60.000, 
		0.050,  0.050,  0.050, 
		4.000,  4.000,  4.000, 
		5.000,  5.000,  5.000, 
		6.000,  6.000,  6.000 
	};

	if ( n_data < 0 )
	{
		n_data = 0;
	}

	n_data = n_data + 1;

	if ( N_MAX < n_data )
	{
		n_data = 0;
		x = 0.0;
		lambda = 0.0;
		df = 0;
		cdf = 0.0;
	}
	else
	{
		x = x_vec[n_data-1];
		lambda = lambda_vec[n_data-1];
		df = df_vec[n_data-1];
		cdf = cdf_vec[n_data-1];
	}

	return;
}

//****************************************************************************80***
//
//  Purpose: 
//
//    CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
//
//  Discussion:
//
//    The value of CHI_CDF ( DF, X ) can be evaluated in Mathematica by 
//    commands like:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[ChiSquareDistribution[DF], X ]
//
//  Modified:
//
//    11 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input/output, int %N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int %A, the parameter of the function.
//
//    Output, double %X, the argument of the function.
//
//    Output, double %FX, the value of the function.
//
//****************************************************************************80***
void StatClass::chi_square_cdf_values ( int %n_data, int %a, double %x, double %fx )
{
	const int N_MAX = 21;

	array<int> ^a_vec = { 
		1,  2,  1,  2, 
		1,  2,  3,  4, 
		1,  2,  3,  4, 
		5,  3,  3,  3, 
		3,  3, 10, 10,
		10 
	};
	array<double> ^fx_vec = { 
		0.0796557, 0.00498752, 0.112463,    0.00995017, 
		0.472911,  0.181269,   0.0597575,   0.0175231, 
		0.682689,  0.393469,   0.198748,    0.090204, 
		0.0374342, 0.427593,   0.608375,    0.738536, 
		0.828203,  0.88839,    0.000172116, 0.00365985, 
		0.0185759 
	};
	array<double> ^x_vec = { 
		0.01, 0.01, 0.02, 0.02, 
		0.40, 0.40, 0.40, 0.40, 
		1.00, 1.00, 1.00, 1.00, 
		1.00, 2.00, 3.00, 4.00, 
		5.00, 6.00, 1.00, 2.00, 
		3.00 
	};

	if ( n_data < 0 )
	{
		n_data = 0;
	}

	n_data = n_data + 1;

	if ( N_MAX < n_data )
	{
		n_data = 0;
		a = 0;
		x = 0.0;
		fx = 0.0;
	}
	else
	{
		a = a_vec[n_data-1];
		x = x_vec[n_data-1];
		fx = fx_vec[n_data-1];
	}
	return;
}


//****************************************************************************80
// 
//  Purpose:
// 
//    CUMBET evaluates the cumulative incomplete beta distribution.
//
//  Discussion:
//
//    This routine calculates the CDF to X of the incomplete beta distribution
//    with parameters A and B.  This is the integral from 0 to x
//    of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)
//
//  Modified:
//
//    14 March 2006
//
//  Reference:
//
//    A R Didonato and Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios. 
//    ACM Transactions on Mathematical Software, 
//    Volume 18, Number 3, September 1992, pages 360-373.
//
//  Parameters:
//
//    Input, double %X, the upper limit of integration.
//
//    Input, double %Y, the value of 1-X.
//
//    Input, double %A, *B, the parameters of the distribution.
//
//    Output, double %CUM, *CCUM, the values of the cumulative
//    density function and complementary cumulative density function.
//
//****************************************************************************80
void StatClass::cumbet ( double %x, double %y, double %a, double %b, double %cum, double %ccum )
{
	int ierr;

	if ( x <= 0.0 )
	{
		cum = 0.0;
		ccum = 1.0;
	}
	else if ( y <= 0.0 ) 
	{
		cum = 1.0;
		ccum = 0.0;
	}
	else
	{
		beta_inc ( a, b, x, y, cum, ccum, ierr );
	}
	return;
}

//****************************************************************************80
// 
//  Purpose:
// 
//    CUMBIN evaluates the cumulative binomial distribution.
//
//  Discussion:
//
//    This routine returns the probability of 0 to S successes in XN binomial
//    trials, each of which has a probability of success, PR.
//
//  Modified:
//
//    14 March 2006
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.5.24.
//
//  Parameters:
//
//    Input, double %S, the upper limit of summation.
//
//    Input, double %XN, the number of trials.
//
//    Input, double %PR, the probability of success in one trial.
//
//    Input, double %OMPR, equals ( 1 - PR ).
//
//    Output, double %CUM, the cumulative binomial distribution.
//
//    Output, double %CCUM, the complement of the cumulative 
//    binomial distribution.
//
//****************************************************************************80

void StatClass::cumbin ( double %s, double %xn, double %pr, double %ompr, double %cum, double %ccum )
{
	double T1,T2;

	if ( s < xn )
	{
		T1 = s + 1.0;
		T2 = xn - s;
		cumbet ( pr, ompr, T1, T2, ccum, cum );
	}
	else
	{
		cum = 1.0;
		ccum = 0.0;
	}
	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    CUMCHI evaluates the cumulative chi-square distribution.
//
//  Parameters:
//
//    Input, double %X, the upper limit of integration.
//
//    Input, double %DF, the degrees of freedom of the
//    chi-square distribution.
//
//    Output, double %CUM, the cumulative chi-square distribution.
//
//    Output, double %CCUM, the complement of the cumulative 
//    chi-square distribution.
//
//****************************************************************************80
void StatClass::cumchi ( double %x, double %df, double %cum, double %ccum )
{
	double a;
	double xx;

	a = df * 0.5;
	xx = x * 0.5;
	cumgam ( xx, a, cum, ccum );
	return;
}

//****************************************************************************80
// 
//  Purpose:
// 
//    CUMCHN evaluates the cumulative noncentral chi-square distribution.
//
//  Discussion:
//
//    Calculates the cumulative noncentral chi-square
//    distribution, i.e., the probability that a random variable
//    which follows the noncentral chi-square distribution, with
//    noncentrality parameter PNONC and continuous degrees of
//    freedom DF, is less than or equal to X.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.4.25.
//
//  Parameters:
//
//    Input, double %X, the upper limit of integration.
//
//    Input, double %DF, the number of degrees of freedom.
//
//    Input, double %PNONC, the noncentrality parameter of 
//    the noncentral chi-square distribution.
//
//    Output, double %CUM, *CCUM, the CDF and complementary
//    CDF of the noncentral chi-square distribution.
//
//  Local Parameters:
//
//    Local, double EPS, the convergence criterion.  The sum 
//    stops when a term is less than EPS*SUM.
//
//    Local, int NTIRED, the maximum number of terms to be evaluated
//    in each sum.
//
//    Local, bool QCONV, is TRUE if convergence was achieved, that is,
//    the program did not stop on NTIRED criterion.
//
//****************************************************************************80

void StatClass::cumchn ( double %x, double %df, double %pnonc, double %cum,	double %ccum )
{
	double eps = 1.0e-5;
	int ntired = 1000;
	double adj,centaj,centwt,chid2,dfd2,lcntaj,lcntwt,lfact,pcent,pterm,sum,
		sumadj,term,wt,xnonc;
	int i,icent,iterb,iterf;
	double T1,T2,T3;

	if (x <= 0.0)
	{
		cum = 0.0;
		ccum = 1.0;
		return;
	}

	if (pnonc <= 1.0e-10)
	{
		//
		//     When non-centrality parameter is (essentially) zero,
		//     use cumulative chi-square distribution
		//
		cumchi(x,df,cum,ccum);
		return;
	}

	xnonc = pnonc/2.0e0;
	//
	//     The following code calculates the weight, chi-square, and
	//     adjustment term for the central term in the infinite series.
	//     The central term is the one in which the poisson weight is
	//     greatest.  The adjustment term is the amount that must
	//     be subtracted from the chi-square to move up two degrees
	//     of freedom.
	//
	icent = fifidint(xnonc);
	if(icent == 0) icent = 1;
	chid2 = x/2.0e0;
	//
	//     Calculate central weight term
	//
	T1 = (double)(icent+1);
	lfact = gamma_log ( T1 );
	lcntwt = -xnonc+(double)icent*Math::Log(xnonc)-lfact;
	centwt = Math::Exp(lcntwt);
	//
	//     Calculate central chi-square
	//
	T2 = (df+2.0e0*(double)(icent));
	cumchi(x,T2,pcent,ccum);
	//
	//     Calculate central adjustment term
	//
	dfd2 = (df+2.0e0*(double)(icent))/2.0e0;
	T3 = 1.0+dfd2;
	lfact = gamma_log ( T3 );
	lcntaj = dfd2*Math::Log(chid2)-chid2-lfact;
	centaj = Math::Exp(lcntaj);
	sum = centwt*pcent;
	//
	//     Sum backwards from the central term towards zero.
	//     Quit whenever either
	//     (1) the zero term is reached, or
	//     (2) the term gets small relative to the sum, or
	//     (3) More than NTIRED terms are totaled.
	//
	iterb = 0;
	sumadj = 0.0;
	adj = centaj;
	wt = centwt;
	i = icent;
	do
	{
		dfd2 = (df+2.0e0*(double)(i))/2.0e0;
		//
		//     Adjust chi-square for two fewer degrees of freedom.
		//     The adjusted value ends up in PTERM.
		//
		adj = adj*dfd2/chid2;
		sumadj = sumadj + adj;
		pterm = pcent+sumadj;
		//
		//     Adjust poisson weight for J decreased by one
		//
		wt *= ((double)i/xnonc);
		term = wt*pterm;
		sum = sum + term;
		i -= 1;
		iterb = iterb + 1;
	}
	while (!((int)((iterb) > ntired) || (int)(sum < 1.0e-20 || (term) < eps*sum) || i == 0));

	iterf = 0;
	//
	//     Now sum forward from the central term towards infinity.
	//     Quit when either
	//     (1) the term gets small relative to the sum, or
	//     (2) More than NTIRED terms are totaled.
	//
	sumadj = adj = centaj;
	wt = centwt;
	i = icent;

	do
	{

		//
		//     Update weights for next higher J
		//
		wt *= (xnonc/(double)(i+1));
		//
		//     Calculate PTERM and add term to sum
		//
		pterm = pcent-sumadj;
		term = wt*pterm;
		sum = sum + term;
		//
		//  Update adjustment term for DF for next iteration
		//
		i = i + 1;
		dfd2 = (df+2.0e0*(double)(i))/2.0e0;
		adj = adj*chid2/dfd2;
		sumadj = sum + adj;
		iterf = iterf + 1;
	}
	while (!((int)((iterf) > ntired) || (int)(sum < 1.0e-20 || (term) < eps*sum)));

	cum = sum;
	ccum = 0.5e0+(0.5e0-cum);
	return;
}

//****************************************************************************80
// 
//  Purpose:
// 
//    CUMF evaluates the cumulative F distribution.
//
//  Discussion:
//
//    CUMF computes the integral from 0 to F of the F density with DFN
//    numerator and DFD denominator degrees of freedom.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.5.28.
//
//  Parameters:
//
//    Input, double %F, the upper limit of integration.
//
//    Input, double %DFN, *DFD, the number of degrees of
//    freedom for the numerator and denominator.
//
//    Output, double %CUM, *CCUM, the value of the F CDF and
//    the complementary F CDF.
//
//****************************************************************************80
void StatClass::cumf ( double %f, double %dfn, double %dfd, double %cum, double %ccum )
{
	double dsum,prod,xx,yy;
	int ierr;
	double T1,T2;

	if (f <= 0.0)
	{
		cum = 0.0;
		ccum = 1.0;
		return;
	}

	prod = dfn*f;
	//
	//     XX is such that the incomplete beta with parameters
	//     DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM
	//     YY is 1 - XX
	//     Calculate the smaller of XX and YY accurately
	//
	dsum = dfd+prod;
	xx = dfd/dsum;

	if ( xx > 0.5e0 )
	{
		yy = prod/dsum;
		xx = 1.0-yy;
	}
	else  
	{
		yy = 1.0-xx;
	}

	T1 = dfd*0.5e0;
	T2 = dfn*0.5e0;
	beta_inc ( T1, T2, xx, yy, ccum, cum, ierr );
	return;
}

//****************************************************************************80
// 
//  Purpose:
//
//    CUMFNC evaluates the cumulative noncentral F distribution.
//
//  Discussion:
//
//    This routine computes the noncentral F distribution with DFN and DFD
//    degrees of freedom and noncentrality parameter PNONC.
//
//    The series is calculated backward and forward from J = LAMBDA/2
//    (this is the term with the largest Poisson weight) until
//    the convergence criterion is met.
//
//    The sum continues until a succeeding term is less than EPS
//    times the sum (or the sum is less than 1.0e-20).  EPS is
//    set to 1.0e-4 in a data statement which can be changed.
//
//
//    The original version of this routine allowed the input values
//    of DFN and DFD to be negative (nonsensical) or zero (which
//    caused numerical overflow.)  I have forced both these values
//    to be at least 1.
//
//  Modified:
//
//    15 June 2004
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.5.16, 26.6.17, 26.6.18, 26.6.20.
//
//  Parameters:
//
//    Input, double %F, the upper limit of integration.
//
//    Input, double %DFN, *DFD, the number of degrees of freedom
//    in the numerator and denominator.  Both DFN and DFD must be positive,
//    and normally would be integers.  This routine requires that they
//    be no less than 1.
//
//    Input, double %PNONC, the noncentrality parameter.
//
//    Output, double %CUM, *CCUM, the noncentral F CDF and 
//    complementary CDF.
//
//****************************************************************************80
void StatClass::cumfnc ( double %f, double %dfn, double %dfd, double %pnonc, double %cum, double %ccum )
{
	double eps = 1.0e-4;
	double dsum,dummy,prod,xx,yy,adn,aup,b,betdn,betup,centwt,dnterm,sum,
		upterm,xmult,xnonc;
	int i,icent,ierr;
	double T1,T2,T3,T4,T5,T6;

	if (f <= 0.0)
	{
		cum = 0.0;
		ccum = 1.0;
		return;
	}

	if (pnonc < 1.0e-10)
	{
		//
		//  Handle case in which the non-centrality parameter is
		//  (essentially) zero.
		//
		cumf(f,dfn,dfd,cum,ccum);
		return;
	}

	xnonc = pnonc/2.0e0;
	//
	//  Calculate the central term of the poisson weighting factor.
	//
	icent = ( int ) xnonc;
	if(icent == 0) icent = 1;
	//
	//  Compute central weight term
	//
	T1 = (double)(icent+1);
	centwt = Math::Exp(-xnonc+(double)icent*Math::Log(xnonc)- gamma_log ( T1 ) );
	//
	//  Compute central incomplete beta term
	//  Assure that minimum of arg to beta and 1 - arg is computed
	//  accurately.
	//
	prod = dfn*f;
	dsum = dfd+prod;
	yy = dfd/dsum;
	if(yy > 0.5e0) {
		xx = prod/dsum;
		yy = 1.0-xx;
	}
	else  xx = 1.0-yy;
	T2 = dfn*0.5e0+(double)icent;
	T3 = dfd*0.5e0;
	beta_inc ( T2, T3, xx, yy, betdn, dummy, ierr );
	adn = dfn/2.0e0+(double)icent;
	aup = adn;
	b = dfd/2.0e0;
	betup = betdn;
	sum = centwt*betdn;
	//
	//  Now sum terms backward from icent until convergence or all done
	//
	xmult = centwt;
	i = icent;
	T4 = adn+b;
	T5 = adn+1.0;
	dnterm = Math::Exp( gamma_log ( T4 ) - gamma_log ( T5 )
		- gamma_log ( b ) + adn * Math::Log ( xx ) + b * Math::Log(yy));

	do
	{
		xmult *= ((double)i/xnonc);
		i -= 1;
		adn -= 1.0;
		dnterm = (adn+1.0)/((adn+b)*xx)*dnterm;
		betdn += dnterm;
		sum += (xmult*betdn);
	}
	while (!((int)(sum < 1.0e-20 || (xmult*betdn) < eps*sum) || i <= 0));

	i = icent+1;
	//
	//  Now sum forwards until convergence
	//
	xmult = centwt;
	if (aup-1.0+b == 0)
		upterm = Math::Exp(-gamma_log ( aup ) - gamma_log ( b ) + (aup-1.0)*Math::Log(xx) +	b*Math::Log(yy));
	else
	{
		T6 = aup-1.0+b;
		upterm = Math::Exp( gamma_log ( T6 ) - gamma_log ( aup ) 
			- gamma_log ( b ) + (aup-1.0)*Math::Log(xx)+b*
			Math::Log(yy));
	}

	do
	{
		xmult *= (xnonc/(double)i);
		i += 1;
		aup += 1.0;
		upterm = (aup+b-2.0e0)*xx/(aup-1.0)*upterm;
		betup -= upterm;
		sum += (xmult*betup);
	}
	while(!((int)(sum < 1.0e-20 || (xmult*betup) < eps*sum)));

	cum = sum;
	ccum = 0.5e0+(0.5e0-cum);
	return;
}

//****************************************************************************80
// 
//  Purpose:
// 
//    CUMGAM evaluates the cumulative incomplete gamma distribution.
//
//  Discussion:
//
//    This routine computes the cumulative distribution function of the 
//    incomplete gamma distribution, i.e., the integral from 0 to X of
//
//      (1/GAM(A))*EXP(-T)*T**(A-1) DT
//
//    where GAM(A) is the complete gamma function of A, i.e.,
//
//      GAM(A) = integral from 0 to infinity of EXP(-T)*T**(A-1) DT
//
//  Parameters:
//
//    Input, double %X, the upper limit of integration.
//
//    Input, double %A, the shape parameter of the incomplete
//    Gamma distribution.
//
//    Output, double %CUM, *CCUM, the incomplete Gamma CDF and
//    complementary CDF.
//
//****************************************************************************80
void StatClass::cumgam ( double %x, double %a, double %cum, double %ccum )
{
	int K1 = 0;

	if (x <= 0.0)
	{
		cum = 0.0;
		ccum = 1.0;
		return;
	}

	cum = ccum = 0;
	gamma_inc ( a, x, cum, ccum, K1 );

	return;
}

//****************************************************************************80
// 
//  Purpose:
// 
//    CUMNBN evaluates the cumulative negative binomial distribution.
//
//  Discussion:
//
//    This routine returns the probability that there will be F or 
//    fewer failures before there are S successes, with each binomial 
//    trial having a probability of success PR.
//
//    Prob(# failures = F | S successes, PR)  =
//                        ( S + F - 1 )
//                        (            ) * PR%S * (1-PR)%F
//                        (      F     )
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.5.26.
//
//  Parameters:
//
//    Input, double %F, the number of failures.
//
//    Input, double %S, the number of successes.
//
//    Input, double %PR, *OMPR, the probability of success on
//    each binomial trial, and the value of (1-PR).
//
//    Output, double %CUM, *CCUM, the negative binomial CDF,
//    and the complementary CDF.
//
//****************************************************************************80
void StatClass::cumnbn ( double %s, double %xn, double %pr, double %ompr, double %cum, double %ccum )
{
	double T1;

	T1 = s + 1.0;
	cumbet(pr,ompr,xn,T1,cum,ccum);
	return;
}

//****************************************************************************80
// 
//  Purpose:
// 
//    CUMNOR computes the cumulative normal distribution.
//
//  Discussion:
//
//    This function evaluates the normal distribution function:
//
//                              / x
//                     1       |       -t*t/2
//          P(x) = ----------- |      e       dt
//                 Math::Sqrt(2 pi)  |
//                             /-oo
//
//    This transportable program uses rational functions that
//    theoretically approximate the normal distribution function to
//    at least 18 significant decimal digits.  The accuracy achieved
//    depends on the arithmetic system, the compiler, the intrinsic
//    functions, and proper selection of the machine-dependent
//    constants.
//
//  Author: 
//
//    William Cody
//    Mathematics and Computer Science Division
//    Argonne National Laboratory
//    Argonne, IL 60439
//
//  Reference:
//
//    William Cody,
//    Rational Chebyshev approximations for the error function,
//    Mathematics of Computation, 
//    1969, pages 631-637.
//
//    William Cody, 
//    Algorithm 715: 
//    SPECFUN - A Portable FORTRAN Package of Special Function Routines 
//      and Test Drivers,
//    ACM Transactions on Mathematical Software,
//    Volume 19, 1993, pages 22-32.
//
//  Parameters:
//
//    Input, double %ARG, the upper limit of integration.
//
//    Output, double %CUM, *CCUM, the Normal density CDF and
//    complementary CDF.
//
//  Local Parameters:
//
//    Local, double EPS, the argument below which anorm(x) 
//    may be represented by 0.5D+00 and above which  x*x  will not underflow.
//    A conservative value is the largest machine number X
//    such that   1.0D+00 + X = 1.0D+00   to machine precision.
//
//****************************************************************************80

void StatClass::cumnor ( double %arg, double %result, double %ccum )
{
	array<double> ^a = {
		2.2352520354606839287e00,1.6102823106855587881e02,1.0676894854603709582e03,
		1.8154981253343561249e04,6.5682337918207449113e-2
	};
	array<double> ^b = {
		4.7202581904688241870e01,9.7609855173777669322e02,1.0260932208618978205e04,
		4.5507789335026729956e04
	};
	array<double> ^c = {
		3.9894151208813466764e-1,8.8831497943883759412e00,9.3506656132177855979e01,
		5.9727027639480026226e02,2.4945375852903726711e03,6.8481904505362823326e03,
		1.1602651437647350124e04,9.8427148383839780218e03,1.0765576773720192317e-8
	};
	array<double> ^d = {
		2.2266688044328115691e01,2.3538790178262499861e02,1.5193775994075548050e03,
		6.4855582982667607550e03,1.8615571640885098091e04,3.4900952721145977266e04,
		3.8912003286093271411e04,1.9685429676859990727e04
	};
	array<double> ^p = {
		2.1589853405795699e-1,1.274011611602473639e-1,2.2235277870649807e-2,
		1.421619193227893466e-3,2.9112874951168792e-5,2.307344176494017303e-2
	};
	double one = 1.0;
	array<double> ^q = {
		1.28426009614491121e00,4.68238212480865118e-1,6.59881378689285515e-2,
		3.78239633202758244e-3,7.29751555083966205e-5
	};
	double sixten = 1.60e0;
	double sqrpi = 3.9894228040143267794e-1;
	double thrsh = 0.66291e0;
	double root32 = 5.656854248e0;
	double zero = 0.0;
	int i;
	double del,eps,temp,x,xden,xnum,y,xsq,min;
	//
	//  Machine dependent constants
	//
	eps = 2.22044604925031E-16 * 0.5;
	min = double::Epsilon;
	x = arg;
	y = Math::Abs(x);
	if(y <= thrsh) {
		//
		//  Evaluate  anorm  for  |X| <= 0.66291
		//
		xsq = zero;
		if(y > eps) xsq = x*x;
		xnum = a[4]*xsq;
		xden = xsq;
		for ( i = 0; i < 3; i++ )
		{
			xnum = (xnum+a[i])*xsq;
			xden = (xden+b[i])*xsq;
		}
		result = x*(xnum+a[3])/(xden+b[3]);
		temp = result;
		result = 0.5e0+temp;
		ccum = 0.5e0-temp;
	}
	//
	//  Evaluate  anorm  for 0.66291 <= |X| <= Math::Sqrt(32)
	//
	else if (y <= root32)
	{
		xnum = c[8] * y;
		xden = y;
		for ( i = 0; i < 7; i++ )
		{
			xnum = (xnum+c[i]) * y;
			xden = (xden+d[i]) * y;
		}
		result = (xnum+c[7])/(xden+d[7]);
		xsq = Math::Truncate(y*sixten)/sixten;
		del = (y-xsq)*(y+xsq);
		result = Math::Exp(-(xsq*xsq*0.5e0))*Math::Exp(-(del*0.5e0))*result;
		ccum = one-result;
		if(x > zero)
		{
			temp = result;
			result = ccum;
			ccum = temp;
		}
	}
	//
	//  Evaluate  anorm  for |X| > Math::Sqrt(32)
	//
	else  {
		result = zero;
		xsq = one/(x*x);
		xnum = p[5]*xsq;
		xden = xsq;
		for ( i = 0; i < 4; i++ )
		{
			xnum = (xnum+p[i])*xsq;
			xden = (xden+q[i])*xsq;
		}
		result = xsq*(xnum+p[4])/(xden+q[4]);
		result = (sqrpi-result)/y;
		xsq = Math::Truncate(x*sixten)/sixten;
		del = (x-xsq)*(x+xsq);
		result = Math::Exp(-(xsq*xsq*0.5e0))*Math::Exp(-(del*0.5e0))*result;
		ccum = one-result;
		if(x > zero) {
			temp = result;
			result = ccum;
			ccum = temp;
		}
	}
	if(result < min) result = 0.0;
	//
	//  Fix up for negative argument, erf, etc.
	//
	if(ccum < min) ccum = 0.0;
}

//****************************************************************************80
// 
//  Purpose:
// 
//    CUMPOI evaluates the cumulative Poisson distribution.
//
//  Discussion:
//
//    CUMPOI returns the probability of S or fewer events in a Poisson
//    distribution with mean XLAM.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    Formula 26.4.21.
//
//  Parameters:
//
//    Input, double %S, the upper limit of cumulation of the 
//    Poisson density function.
//
//    Input, double %XLAM, the mean of the Poisson distribution.
//
//    Output, double %CUM, *CCUM, the Poisson density CDF and 
//    complementary CDF.
//
//****************************************************************************80
void StatClass::cumpoi ( double %s, double %xlam, double %cum, double %ccum )
{
	double chi,df;

	df = 2.0 * (s + 1.0);
	chi = 2.0 * xlam;
	cumchi(chi,df,ccum,cum);
	return;
}

//****************************************************************************80
// 
//  Purpose:
// 
//    CUMT evaluates the cumulative T distribution.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    Formula 26.5.27.
//
//  Parameters:
//
//    Input, double %T, the upper limit of integration.
//
//    Input, double %DF, the number of degrees of freedom of 
//    the T distribution.
//
//    Output, double %CUM, *CCUM, the T distribution CDF and
//    complementary CDF.
//
//****************************************************************************80
void StatClass::cumt ( double %t, double %df, double %cum, double %ccum )
{
	double a;
	double dfptt;
	double K2 = 0.5e0;
	double oma;
	double T1;
	double tt;
	double xx;
	double yy;

	tt = (t) * (t);
	dfptt = ( df ) + tt;
	xx = df / dfptt;
	yy = tt / dfptt;
	T1 = 0.5e0 * ( df );
	cumbet ( xx, yy, T1, K2, a, oma );

	if ( t <= 0.0 )
	{
		cum = 0.5e0 * a;
		ccum = oma + ( cum );
	}
	else
	{
		ccum = 0.5e0 * a;
		cum = oma + ( ccum );
	}
	return;
}


//****************************************************************************80
// 
//  Purpose:
// 
//    DBETRM computes the Sterling remainder for the complete beta function.
//
//  Discussion:
//
//    Math::Log(Beta(A,B)) = Lgamma(A) + Lgamma(B) - Lgamma(A+B)
//    where Lgamma is the Math::Log of the (complete) gamma function
//
//    Let ZZ be approximation obtained if each Math::Log gamma is approximated
//    by Sterling's formula, i.e.,
//    Sterling(Z) = LOG( SQRT( 2*PI ) ) + ( Z-0.5D+00 ) * LOG( Z ) - Z
//
//    The Sterling remainder is Math::Log(Beta(A,B)) - ZZ.
//
//  Parameters:
//
//    Input, double %A, *B, the parameters of the Beta function.
//
//    Output, double DBETRM, the Sterling remainder.
//
//****************************************************************************80
double StatClass::dbetrm ( double %a, double %b )
{
	double dbetrm,T1,T2,T3;
	//
	//     Try to sum from smallest to largest
	//
	T1 = a+b;
	dbetrm = -dstrem(T1);
	T2 = Math::Max(a,b);
	dbetrm += dstrem(T2);
	T3 = Math::Min(a,b);
	dbetrm += dstrem(T3);
	return dbetrm;
}


//****************************************************************************80
// 
//  Purpose:
// 
//    DEXPM1 evaluates the function EXP(X) - 1.
//
//  Reference:
//
//    Armido DiDinato and Alfred Morris,
//    Algorithm 708: 
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software, 
//    Volume 18, 1993, pages 360-373.
//
//  Parameters:
//
//    Input, double %X, the value at which Math::Exp(X)-1 is desired.
//
//    Output, double DEXPM1, the value of Math::Exp(X)-1.
//
//****************************************************************************80
double StatClass::dexpm1 ( double %x )
{
	double p1 = .914041914819518e-09;
	double p2 = .238082361044469e-01;
	double q1 = -.499999999085958e+00;
	double q2 = .107141568980644e+00;
	double q3 = -.119041179760821e-01;
	double q4 = .595130811860248e-03;
	double dexpm1;
	double w;

	if ( Math::Abs(x) <= 0.15e0 )
	{
		dexpm1 =   x * ( ( ( 
			p2   * x 
			+ p1 ) * x 
			+ 1.0 )
			/(((( 
			q4   * x 
			+ q3 ) * x 
			+ q2 ) * x 
			+ q1 ) * x 
			+ 1.0 ) );
	}
	else if ( x <= 0.0 )
	{
		w = Math::Exp(x);
		dexpm1 = w-0.5e0-0.5e0;
	}
	else
	{
		w = Math::Exp(x);
		dexpm1 = w*(0.5e0+(0.5e0-1.0/w));
	}

	return dexpm1;
}

//****************************************************************************80
// 
//  Purpose:
// 
//    DINVNR computes the inverse of the normal distribution.
//
//  Discussion:
//
//    Returns X such that CUMNOR(X)  =   P,  i.e., the  integral from -
//    infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
//
//    The rational function on page 95 of Kennedy and Gentle is used as a start
//    value for the Newton method of finding roots.
//
//  Reference:
//
//    Kennedy and Gentle, 
//    Statistical Computing,
//    Marcel Dekker, NY, 1980,
//    QA276.4  K46
//
//  Parameters:
//
//    Input, double %P, *Q, the probability, and the complementary
//    probability.
//
//    Output, double DINVNR, the argument X for which the
//    Normal CDF has the value P.
//
//****************************************************************************80
double StatClass::dinvnr ( double %p, double %q )
{
	const int maxit = 100;
	const double r2pi = 0.3989422804014326e0;
	double dinvnr,strtx,xcur,cum,ccum,pp,dx;

	//
	//     FIND MINIMUM OF P AND Q
	//
	bool qporq = p <= q;
	if (qporq)
		pp = p;
	else
		pp = q;

	//
	//     INITIALIZATION STEP
	//
	strtx = stvaln(pp);
	xcur = strtx;
	//
	//     NEWTON INTERATIONS
	//
	for (int i = 1; i <= maxit; i++ )
	{
		cumnor(xcur,cum,ccum);
		dx = (cum-pp)/(r2pi*Math::Exp(-0.5e0*(xcur)*(xcur)));
		xcur -= dx;
		if(Math::Abs(dx/xcur) < (1.0e-13))
		{
			//
			//     IF WE GET HERE, NEWTON HAS SUCCEDED
			//
			dinvnr = xcur;
			if (!qporq)
				dinvnr = -dinvnr;
			return dinvnr;
		}
	}

	dinvnr = strtx;
	//
	//     IF WE GET HERE, NEWTON HAS FAILED
	//
	if (!qporq)
		dinvnr = -dinvnr;
	return dinvnr;

}

//****************************************************************************80
//
//  Purpose:
// 
//    DINVR bounds the zero of the function and invokes DZROR.
//
//  Discussion:
//
//    This routine seeks to find bounds on a root of the function and 
//    invokes ZROR to perform the zero finding.  STINVR must have been 
//    called before this routine in order to set its parameters.
//
//  Reference:
//
//    J C P Bus and T J Dekker,
//    Two Efficient Algorithms with Guaranteed Convergence for 
//      Finding a Zero of a Function,
//    ACM Transactions on Mathematical Software,
//    Volume 1, Number 4, pages 330-345, 1975.
//
//  Parameters:
//
//    Input/output, integer STATUS.  At the beginning of a zero finding 
//    problem, STATUS should be set to 0 and INVR invoked.  The value
//    of parameters other than X will be ignored on this call.
//    If INVR needs the function to be evaluated, it will set STATUS to 1 
//    and return.  The value of the function should be set in FX and INVR 
//    again called without changing any of its other parameters.
//    If INVR finishes without error, it returns with STATUS 0, and X an
//    approximate root of F(X).
//    If INVR cannot bound the function, it returns a negative STATUS and 
//    sets QLEFT and QHI.
//
//    Output, double precision X, the value at which F(X) is to be evaluated.
//
//    Input, double precision FX, the value of F(X) calculated by the user
//    on the previous call, when INVR returned with STATUS = 1.
//
//    Output, logical QLEFT, is defined only if QMFINV returns FALSE.  In that
//    case, QLEFT is TRUE if the stepping search terminated unsucessfully 
//    at SMALL, and FALSE if the search terminated unsucessfully at BIG.
//
//    Output, logical QHI, is defined only if QMFINV returns FALSE.  In that
//    case, it is TRUE if Y < F(X) at the termination of the search and FALSE
//    if F(X) < Y.
//
//****************************************************************************80
void StatClass::dinvr ( int %status, double %x, double %fx,
					   unsigned long %qleft, unsigned long %qhi )

{
	double nothing = 0.0;
	E0000(0,status,x,fx,qleft,qhi,nothing,nothing,nothing,nothing,nothing,nothing,nothing);
}

//****************************************************************************80
//
//  Purpose:
// 
//    DLANOR evaluates the logarithm of the asymptotic Normal CDF.
//
//  Discussion:
//
//    This routine computes the logarithm of the cumulative normal distribution
//    from abs ( x ) to infinity for  5 <= abs ( X ).
//
//    The relative error at X = 5 is about 0.5D-5.
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions 
//    1966, Formula 26.2.12.
//
//  Parameters:
//
//    Input, double %X, the value at which the Normal CDF is to be
//    evaluated.  It is assumed that 5 <= abs ( X ).
//
//    Output, double DLANOR, the logarithm of the asymptotic
//    Normal CDF.
//
//****************************************************************************80
double StatClass::dlanor ( double %x )

{
	const double dlsqpi = 0.91893853320467274177e0;

	array<double> ^coef = {
		-1.0,3.0,-15.0,105.0,-945.0,10395.0,-135135.0,2027025.0,
		-34459425.0,654729075.0,-13749310575.0,316234143225.0
	};
	int K1 = 12;
	double dlanor,approx,correc,xx,xx2,T2;

	xx = Math::Abs(x);
	if ( xx < 5.0 )
		throw gcnew ArgumentException("Argument less than 5 in DLANOR", "x");
	approx = -dlsqpi-0.5*xx*xx-Math::Log(xx);
	xx2 = xx*xx;
	T2 = 1.0/xx2;
	correc = eval_pol ( coef, K1, T2 ) / xx2;
	correc = alnrel ( correc );
	dlanor = approx+correc;
	return dlanor;
}

//*****************************************************************************
//
//  Purpose:
//
//    DSTINV seeks a value X such that F(X) = Y.
//
//  Discussion:
//
//      Double Precision - SeT INverse finder - Reverse Communication
//                              Function
//     Concise Description - Given a monotone function F finds X
//     such that F(X) = Y.  Uses Reverse communication -- see invr.
//     This routine sets quantities needed by INVR.
//          More Precise Description of INVR -
//     F must be a monotone function, the results of QMFINV are
//     otherwise undefined.  QINCR must be .TRUE. if F is non-
//     decreasing and .FALSE. if F is non-increasing.
//     QMFINV will return .TRUE. if and only if F(SMALL) and
//     F(BIG) bracket Y, i. e.,
//          QINCR is .TRUE. and F(SMALL).LE.Y.LE.F(BIG) or
//          QINCR is .FALSE. and F(BIG).LE.Y.LE.F(SMALL)
//     if QMFINV returns .TRUE., then the X returned satisfies
//     the following condition.  let
//               TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
//     then if QINCR is .TRUE.,
//          F(X-TOL(X)) .LE. Y .LE. F(X+TOL(X))
//     and if QINCR is .FALSE.
//          F(X-TOL(X)) .GE. Y .GE. F(X+TOL(X))
//                              Arguments
//     SMALL --> The left endpoint of the interval to be
//          searched for a solution.
//                    SMALL is DOUBLE PRECISION
//     BIG --> The right endpoint of the interval to be
//          searched for a solution.
//                    BIG is DOUBLE PRECISION
//     ABSSTP, RELSTP --> The initial step size in the search
//          is MAX(ABSSTP,RELSTP*ABS(X)). See algorithm.
//                    ABSSTP is DOUBLE PRECISION
//                    RELSTP is DOUBLE PRECISION
//     STPMUL --> When a step doesn't bound the zero, the step
//                size is multiplied by STPMUL and another step
//                taken.  A popular value is 2.0
//                    DOUBLE PRECISION STPMUL
//     ABSTOL, RELTOL --> Two numbers that determine the accuracy
//          of the solution.  See function for a precise definition.
//                    ABSTOL is DOUBLE PRECISION
//                    RELTOL is DOUBLE PRECISION
//                              Method
//     Compares F(X) with Y for the input value of X then uses QINCR
//     to determine whether to step left or right to bound the
//     desired x.  the initial step size is
//          MAX(ABSSTP,RELSTP*ABS(S)) for the input value of X.
//     Iteratively steps right or left until it bounds X.
//     At each step which doesn't bound X, the step size is doubled.
//     The routine is careful never to step beyond SMALL or BIG.  If
//     it hasn't bounded X at SMALL or BIG, QMFINV returns .FALSE.
//     after setting QLEFT and QHI.
//     If X is successfully bounded then Algorithm R of the paper
//     'Two Efficient Algorithms with Guaranteed Convergence for
//     Finding a Zero of a Function' by J. C. P. Bus and
//     T. J. Dekker in ACM Transactions on Mathematical
//     Software, Volume 1, No. 4 page 330 (DEC. '75) is employed
//     to find the zero of the function F(X)-Y. This is routine
//     QRZERO.
//
//*****************************************************************************
void StatClass::dstinv (double %zsmall, double %zbig, double %zabsst,
						double %zrelst, double %zstpmu, double %zabsto, double %zrelto )
{
	int nnothing = 0;
	double nothing = 0.0;
	unsigned long lnothing = 0;
	E0000( 1, nnothing, nothing, nothing, lnothing, lnothing,
		zabsst,	zabsto, zbig, zrelst, zrelto, zsmall, zstpmu);
}
//*****************************************************************************

//****************************************************************************80
//
//  Purpose:
// 
//    DSTREM computes the Sterling remainder ln ( Gamma ( Z ) ) - Sterling ( Z ).
//
//  Discussion:
//
//    This routine returns 
//
//      ln ( Gamma ( Z ) ) - Sterling ( Z )  
//
//    where Sterling(Z) is Sterling's approximation to ln ( Gamma ( Z ) ).
//
//    Sterling(Z) = ln ( Math::Sqrt ( 2 * PI ) ) + ( Z - 0.5 ) * ln ( Z ) - Z
//
//    If 6 <= Z, the routine uses 9 terms of a series in Bernoulli numbers,
//    with values calculated using Maple.
//
//    Otherwise, the difference is computed explicitly.
//
//  Modified:
//
//    14 June 2004
//
//  Parameters:
//
//    Input, double %Z, the value at which the Sterling 
//    remainder is to be calculated.  Z must be positive.
//
//    Output, double DSTREM, the Sterling remainder.
//
double StatClass::dstrem ( double %z )
{
	array<double> ^coef =
	{
		0.0,0.0833333333333333333333333333333,
		-0.00277777777777777777777777777778,
		0.000793650793650793650793650793651,
		-0.000595238095238095238095238095238,
		0.000841750841750841750841750841751,
		-0.00191752691752691752691752691753,
		0.00641025641025641025641025641026,
		-0.0295506535947712418300653594771,
		0.179644372368830573164938490016
	};

	double dstrem = 0.0;
	//
	//    For information, here are the next 11 coefficients of the
	//    remainder term in Sterling's formula
	//            -1.39243221690590111642743221691
	//            13.4028640441683919944789510007
	//            -156.848284626002017306365132452
	//            2193.10333333333333333333333333
	//            -36108.7712537249893571732652192
	//            691472.268851313067108395250776
	//            -0.152382215394074161922833649589D8
	//            0.382900751391414141414141414141D9
	//            -0.108822660357843910890151491655D11
	//            0.347320283765002252252252252252D12
	//            -0.123696021422692744542517103493D14
	//
	if (z <= 0.0)
		throw gcnew ArgumentException("dstrem: Argument less than zero.", "z");

	if (z > 6.0e0)
	{
		int K1 = 10;
		double T2 = 1.0 / (z * z);
		dstrem = eval_pol ( coef, K1, T2 )*z;
	}
	else
	{
		const double hln2pi = 0.91893853320467274178;
		double sterl = hln2pi + (z - 0.5) * Math::Log(z)-z;
		dstrem = gamma_log ( z ) - sterl;
	}

	return dstrem;
}

//*****************************************************************************
//
//  Purpose:
// 
//    DSTXR sets quantities needed by the zero finder.
//
//  Discussion:
//
//     Double precision SeT ZeRo finder - Reverse communication version
//                              Function
//     Sets quantities needed by ZROR.  The function of ZROR
//     and the quantities set is given here.
//     Concise Description - Given a function F
//     find XLO such that F(XLO) = 0.
//          More Precise Description -
//     Input condition. F is a double function of a single
//     double argument and XLO and XHI are such that
//          F(XLO)*F(XHI)  .LE.  0.0
//     If the input condition is met, QRZERO returns .TRUE.
//     and output values of XLO and XHI satisfy the following
//          F(XLO)*F(XHI)  .LE. 0.
//          ABS(F(XLO)  .LE. ABS(F(XHI)
//          ABS(XLO-XHI)  .LE. TOL(X)
//     where
//          TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
//     If this algorithm does not find XLO and XHI satisfying
//     these conditions then QRZERO returns .FALSE.  This
//     implies that the input condition was not met.
//                              Arguments
//     XLO --> The left endpoint of the interval to be
//           searched for a solution.
//                    XLO is DOUBLE PRECISION
//     XHI --> The right endpoint of the interval to be
//           for a solution.
//                    XHI is DOUBLE PRECISION
//     ABSTOL, RELTOL --> Two numbers that determine the accuracy
//                      of the solution.  See function for a
//                      precise definition.
//                    ABSTOL is DOUBLE PRECISION
//                    RELTOL is DOUBLE PRECISION
//                              Method
//     Algorithm R of the paper 'Two Efficient Algorithms with
//     Guaranteed Convergence for Finding a Zero of a Function'
//     by J. C. P. Bus and T. J. Dekker in ACM Transactions on
//     Mathematical Software, Volume 1, no. 4 page 330
//     (Dec. '75) is employed to find the zero of F(X)-Y.
//
//*****************************************************************************
void StatClass::dstzr ( double %zxlo, double %zxhi, double %zabstl, double %zreltl )
{
	int nnothing = 0;
	double nothing = 0.0;
	unsigned long lnothing = 0;
	E0001(1,nnothing,nothing,nothing,nothing,nothing,lnothing,lnothing,zabstl,zreltl,zxhi,zxlo);
}

//****************************************************************************80
//
//  Purpose:
// 
//    DT1 computes an approximate inverse of the cumulative T distribution.
//
//  Discussion:
//
//    Returns the inverse of the T distribution function, i.e.,
//    the integral from 0 to INVT of the T density is P. This is an
//    initial approximation.
//
//  Parameters:
//
//    Input, double %P, *Q, the value whose inverse from the 
//    T distribution CDF is desired, and the value (1-P).
//
//    Input, double %DF, the number of degrees of freedom of the 
//    T distribution.
//
//    Output, double DT1, the approximate value of X for which
//    the T density CDF with DF degrees of freedom has value P.
//
//****************************************************************************80
double StatClass::dt1 ( double %p, double %q, double %df )
{
	array<double, 2> ^coef = 
	{
		{   1.0,	    1.0,	   0.0,	  0.0,	 0.0},
		{   3.0,	   16.0,	   5.0,	  0.0,	 0.0},
		{ -15.0,	   17.0,	  19.0,	  3.0,	 0.0},
		{-945.0,	-1920.0,	1482.0,	776.0,	79.0}
	};
	array<double> ^denom = 
	{
		4.0,	96.0,	384.0,	92160.0
	};

	array<int> ^ideg =
	{
		2,3,4,5
	};

	double x = Math::Abs(dinvnr(p,q));
	double xx = x*x;
	double sum = x;
	double denpow = 1.0;

	for (int i = 0; i < 4; i++ )
	{
		array<double> ^cf = {coef[i,0], coef[i,1], coef[i,2], coef[i,3], coef[i,4]};
		double term = eval_pol ( cf, ideg[i], xx ) * x;
		denpow *= df;
		sum += ( term / (denpow*denom[i]) );
	}

	if (p >= 0.5)
		return sum;
	else
		return -sum;
}

//****************************************************************************80
// 
//  Purpose:
// 
//    DZROR seeks the zero of a function using reverse communication.
//
//  Discussion:
//
//     Performs the zero finding.  STZROR must have been called before
//     this routine in order to set its parameters.
// 
// 
//                              Arguments
// 
// 
//     STATUS <--> At the beginning of a zero finding problem, STATUS
//                 should be set to 0 and ZROR invoked.  (The value
//                 of other parameters will be ignored on this call.)
// 
//                 When ZROR needs the function evaluated, it will set
//                 STATUS to 1 and return.  The value of the function
//                 should be set in FX and ZROR again called without
//                 changing any of its other parameters.
// 
//                 When ZROR has finished without error, it will return
//                 with STATUS 0.  In that case (XLO,XHI) bound the answe
// 
//                 If ZROR finds an error (which implies that F(XLO)-Y an
//                 F(XHI)-Y have the same sign, it returns STATUS -1.  In
//                 this case, XLO and XHI are undefined.
//                         INTEGER STATUS
// 
//     X <-- The value of X at which F(X) is to be evaluated.
//                         DOUBLE PRECISION X
// 
//     FX --> The value of F(X) calculated when ZROR returns with
//            STATUS = 1.
//                         DOUBLE PRECISION FX
// 
//     XLO <-- When ZROR returns with STATUS = 0, XLO bounds the
//             inverval in X containing the solution below.
//                         DOUBLE PRECISION XLO
// 
//     XHI <-- When ZROR returns with STATUS = 0, XHI bounds the
//             inverval in X containing the solution above.
//                         DOUBLE PRECISION XHI
// 
//     QLEFT <-- .TRUE. if the stepping search terminated unsucessfully
//                at XLO.  If it is .FALSE. the search terminated
//                unsucessfully at XHI.
//                    QLEFT is LOGICAL
// 
//     QHI <-- .TRUE. if F(X) .GT. Y at the termination of the
//              search and .FALSE. if F(X) .LT. Y at the
//              termination of the search.
//                    QHI is LOGICAL
// 
//
//****************************************************************************80
void StatClass::dzror ( int %status, double %x, double %fx, double %xlo, double %xhi, unsigned long %qleft, unsigned long %qhi )
{
	double dnothing = 0.0;
	E0001(0, status, x, fx, xlo, xhi, qleft, qhi, dnothing, dnothing, dnothing, dnothing);
}

//*****************************************************************************
//
//  Purpose:
//
//    E0000 is a reverse-communication zero bounder.
//
//*****************************************************************************
void StatClass::E0000 ( int IENTRY,
					   int %status,
					   double %x,
					   double %fx,
					   unsigned long %qleft,
					   unsigned long %qhi,
					   double %zabsst,
					   double %zabsto,
					   double %zbig,
					   double %zrelst,
					   double %zrelto,
					   double %zsmall,
					   double %zstpmu )
{
	switch(IENTRY)
	{
	case 0:
		if (status > 0)
		{
			switch(E0000_lastState)
			{
			case 1:
				E0000_fsmall = fx;
				x = E0000_big;
				//
				//     GET-FUNCTION-VALUE
				//
				E0000_lastState = 2;
				//
				//     TO GET-FUNCTION-VALUE
				//
				status = 1;
				return;

			case 2: 
				E0000_fbig = fx;
				E0000_qincr = E0000_fbig > E0000_fsmall;

				if (E0000_qincr)
				{
					if (E0000_fsmall > 0.0)
					{
						status = -1;
						qleft = qhi = 1;
						return;
					}

					if (E0000_fbig < 0.0)
					{
						status = -1;
						qleft = qhi = 0;
						return;
					}

					x = E0000_xsave;
					E0000_step = Math::Max(E0000_absstp,E0000_relstp*Math::Abs(x));
					//
					//      YY = F(X) - Y
					//     GET-FUNCTION-VALUE
					//
					E0000_lastState = 3;
					//
					//     TO GET-FUNCTION-VALUE
					//
					status = 1;
					return;
				}

				if (E0000_fsmall < 0.0)
				{
					status = -1;
					qleft = 1;
					qhi = 0;
					return;
				}

				if (E0000_fbig > 0.0)
				{
					status = -1;
					qleft = 0;
					qhi = 1;
					return;
				}

				x = E0000_xsave;
				E0000_step = Math::Max(E0000_absstp,E0000_relstp*Math::Abs(x));
				//
				//      YY = F(X) - Y
				//     GET-FUNCTION-VALUE
				//
				E0000_lastState = 3;
				//
				//     TO GET-FUNCTION-VALUE
				//
				status = 1;
				return;

			case 3: 

				E0000_yy = fx;

				if (E0000_yy == 0.0)
				{
					status = 0;
					E0000_qok = 1;
					return;
				}

				E0000_qup = E0000_qincr && E0000_yy < 0.0 || !E0000_qincr && E0000_yy > 0.0;

				if (!E0000_qup)
				{
					//
					//     HANDLE CASE IN WHICH WE MUST STEP LOWER
					//
					E0000_xub = E0000_xsave;
					E0000_xlb = Math::Max(E0000_xub-E0000_step,E0000_small);
					//
					//      YY = F(XLB) - Y
					//
					x = E0000_xlb;
					//
					//     GET-FUNCTION-VALUE
					//
					E0000_lastState = 5;
					//
					//     TO GET-FUNCTION-VALUE
					//
					status = 1;
					return;
				}
				else
				{
					//
					//     HANDLE CASE IN WHICH WE MUST STEP HIGHER
					//
					E0000_xlb = E0000_xsave;
					E0000_xub = Math::Min(E0000_xlb+E0000_step,E0000_big);
					//
					//      YY = F(XUB) - Y
					//
					x = E0000_xub;
					//
					//     GET-FUNCTION-VALUE
					//
					E0000_lastState = 4;
					//
					//     TO GET-FUNCTION-VALUE
					//
					status = 1;
					return;
				}

			case 4: 
				E0000_yy = fx;
				E0000_qbdd = E0000_qincr && E0000_yy >= 0.0 || !E0000_qincr && E0000_yy <= 0.0;
				E0000_qlim = E0000_xub >= E0000_big;
				E0000_qcond = E0000_qbdd || E0000_qlim;

				if (!E0000_qcond)
				{
					E0000_step = E0000_stpmul*E0000_step;
					E0000_xlb = E0000_xub;
					E0000_xub = Math::Min(E0000_xlb+E0000_step,E0000_big);
					//
					//      YY = F(XUB) - Y
					//
					x = E0000_xub;
					//
					//     GET-FUNCTION-VALUE
					//
					E0000_lastState = 4;
					//
					//     TO GET-FUNCTION-VALUE
					//
					status = 1;
					return;
				}

				if (E0000_qlim || !E0000_qbdd)
				{
					status = -1;
					qleft = 0;
					qhi = !E0000_qincr;
					x = E0000_big;
					return;
				}

				dstzr(E0000_xlb,E0000_xub,E0000_abstol,E0000_reltol);
				//
				//  IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
				//
				status = 0;
				dzror ( status, x, fx, E0000_xlo, E0000_xhi, E0000_qdum1, E0000_qdum2 );
				if (status != 1)
				{
					x = E0000_xlo;
					status = 0;
					return;
				}
				else
				{
					//
					//     GET-FUNCTION-VALUE
					//
					E0000_lastState = 6;
					//
					//     TO GET-FUNCTION-VALUE
					//
					status = 1;
					return;
				}

			case 5:
				E0000_yy = fx;
				E0000_qbdd = E0000_qincr && E0000_yy <= 0.0 || !E0000_qincr && E0000_yy >= 0.0;
				E0000_qlim = E0000_xlb <= E0000_small;
				E0000_qcond = E0000_qbdd || E0000_qlim;

				if (!E0000_qcond)
				{
					E0000_step = E0000_stpmul*E0000_step;
					E0000_xub = E0000_xlb;
					E0000_xlb = Math::Max(E0000_xub-E0000_step,E0000_small);

					//
					//      YY = F(XLB) - Y
					//
					x = E0000_xlb;
					//
					//     GET-FUNCTION-VALUE
					//
					E0000_lastState = 5;
					//
					//     TO GET-FUNCTION-VALUE
					//
					status = 1;
					return;
				}

				if (E0000_qlim && !E0000_qbdd)
				{
					status = -1;
					qleft = 1;
					qhi = E0000_qincr;
					x = E0000_small;
					return;
				}

				dstzr(E0000_xlb,E0000_xub,E0000_abstol,E0000_reltol);
				//
				//  IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
				//
				status = 0;

				dzror ( status, x, fx, E0000_xlo, E0000_xhi, E0000_qdum1, E0000_qdum2 );
				if(!(status == 1))
				{
					x = E0000_xlo;
					status = 0;
					return;
				}

				//
				//     GET-FUNCTION-VALUE
				//
				E0000_lastState = 6;
				//
				//     TO GET-FUNCTION-VALUE
				//
				status = 1;
				return;

			case 6:
				if (status != 1)
				{
					x = E0000_xlo;
					status = 0;
					return;
				}
				dzror ( status, x, fx, E0000_xlo, E0000_xhi, E0000_qdum1, E0000_qdum2 );
				if(!(status == 1))
				{
					x = E0000_xlo;
					status = 0;
					return;
				}

				//
				//     GET-FUNCTION-VALUE
				//
				E0000_lastState = 6;
				//
				//     TO GET-FUNCTION-VALUE
				//
				status = 1;
				return;

			default: 
				throw gcnew ArgumentOutOfRangeException("E0000_lastState","Must be in the range of 1-6 inclusive.");
			}
		}
		else
		{
			E0000_qcond = !(int)((E0000_small) <= (x) && (x) <= (E0000_big));
			if (E0000_qcond)
				throw gcnew ApplicationException ("Arguments E0000_small, X, E0000_big not monotone in INVR");
			E0000_xsave = x;
			//
			//     See that SMALL and BIG bound the zero and set QINCR
			//
			x = E0000_small;
			//
			//     GET-FUNCTION-VALUE
			//
			E0000_lastState = 1;
			//
			//     TO GET-FUNCTION-VALUE
			//
			status = 1;
			return;
		}

	case 1:
		{
			E0000_small = zsmall;
			E0000_big = zbig;
			E0000_absstp = zabsst;
			E0000_relstp = zrelst;
			E0000_stpmul = zstpmu;
			E0000_abstol = zabsto;
			E0000_reltol = zrelto;
			return;
		}

	default: 
		throw gcnew System::ArgumentOutOfRangeException("IENTRY", "Must have a value 0 or 1");
	}
}

//****************************************************************************80
//
//  Purpose:
//
//    E00001 is a reverse-communication zero finder.
//
//****************************************************************************80
void StatClass::E0001 ( int IENTRY,
					   int %status,
					   double %x,
					   double %fx,
					   double %xlo,
					   double %xhi,
					   unsigned long %qleft,
					   unsigned long %qhi,
					   double %zabstl,
					   double %zreltl,
					   double %zxhi,
					   double %zxlo )
{
	switch(IENTRY)
	{
	case 0:
		if (status <= 0)
		{
			xlo = E0001_xxlo;
			xhi = E0001_xxhi;
			E0001_b = x = xlo;
			//
			//     GET-FUNCTION-VALUE
			//
			E0001_lastState = 1;
			//
			//     TO GET-FUNCTION-VALUE
			//
			status = 1;
			return;
		}

		switch(E0001_lastState)
		{
		case 1:
			E0001_fb = fx;
			xlo = xhi;
			E0001_a = x = xlo;
			//
			//     GET-FUNCTION-VALUE
			//
			E0001_lastState = 2;
			//
			//     TO GET-FUNCTION-VALUE
			//
			status = 1;
			return;

		case 2:
			//
			//     Check that F(ZXLO) < 0 < F(ZXHI)  or
			//                F(ZXLO) > 0 > F(ZXHI)
			//
			if (E0001_fb < 0.0)
			{
				if (fx < 0.0)
				{
					status = -1;
					qleft = fx < E0001_fb;
					qhi = 0;
					return;
				}
			}

			if (E0001_fb > 0.0)
			{
				if (fx > 0.0)
				{
					status = -1;
					qleft = fx > E0001_fb;
					qhi = 1;
					return;
				}
			}

			E0001_fa = fx;
			E0001_first = 1;
			E0001_c = E0001_a;
			E0001_fc = E0001_fa;
			E0001_ext = 0;

			break;

		case 3:

			E0001_fb = fx;

			if (E0001_fc*E0001_fb >= 0.0)
			{
				E0001_c = E0001_a;
				E0001_fc = E0001_fa;
				E0001_ext = 0;
			}
			else if (E0001_w == E0001_mb)
			{
				E0001_ext = 0;
			}
			else
			{
				E0001_ext += 1;
			}

			break;

		default: 
			throw gcnew System::ArgumentOutOfRangeException("E0001_lastState", "Must have a value between 1 or 3");
		}

		if (Math::Abs(E0001_fc) < Math::Abs(E0001_fb))
		{
			if (E0001_c != E0001_a)
			{
				E0001_d = E0001_a;
				E0001_fd = E0001_fa;
			}

			E0001_a = E0001_b;
			E0001_fa = E0001_fb;
			xlo = E0001_c;
			E0001_b = xlo;
			E0001_fb = E0001_fc;
			E0001_c = E0001_a;
			E0001_fc = E0001_fa;
		}

		E0001_tol = (0.5*Math::Max(E0001_abstol,E0001_reltol*Math::Abs((xlo))));
		E0001_m = (E0001_c+E0001_b)*.5;
		E0001_mb = E0001_m-E0001_b;

		if(!(Math::Abs(E0001_mb) > E0001_tol))
		{
			xhi = E0001_c;
			E0001_qrzero = E0001_fc >= 0.0 && E0001_fb <= 0.0 || E0001_fc < 0.0 && E0001_fb >= 0.0;
			if (E0001_qrzero)
				status = 0;
			else
				status = -1;
			return;
		}

		if (E0001_ext > 3)
		{
			E0001_w = E0001_mb;
			E0001_d = E0001_a;
			E0001_fd = E0001_fa;
			E0001_a = E0001_b;
			E0001_fa = E0001_fb;
			E0001_b += E0001_w;
			xlo = E0001_b;
			x = xlo;
			//
			//     GET-FUNCTION-VALUE
			//
			E0001_lastState = 3;
			//
			//     TO GET-FUNCTION-VALUE
			//
			status = 1;
			return;
		}

		E0001_tol = fifdsign(E0001_tol,E0001_mb);
		E0001_p = (E0001_b-E0001_a)*E0001_fb;

		if (E0001_first)
		{
			E0001_q = E0001_fa-E0001_fb;
			E0001_first = 0;
		}
		else
		{
			E0001_fdb = (E0001_fd - E0001_fb) / (E0001_d - E0001_b);
			E0001_fda = (E0001_fd - E0001_fa) / (E0001_d - E0001_a);
			E0001_p = E0001_fda * E0001_p;
			E0001_q = E0001_fdb * E0001_fa - E0001_fda * E0001_fb;
		}

		if (E0001_p < 0.0)
		{
			E0001_p = -E0001_p;
			E0001_q = -E0001_q;
		}

		if (E0001_ext == 3)
			E0001_p *= 2.0;
		if (E0001_p == 0.0 || E0001_p <= E0001_q * E0001_tol)
		{
			E0001_w = E0001_tol;
		}
		else if (E0001_p < E0001_mb * E0001_q)
		{
			E0001_w = E0001_p/E0001_q;
		}
		else
		{
			E0001_w = E0001_mb;
		}

		E0001_d = E0001_a;
		E0001_fd = E0001_fa;
		E0001_a = E0001_b;
		E0001_fa = E0001_fb;
		E0001_b += E0001_w;
		xlo = E0001_b;
		x = xlo;
		//
		//     GET-FUNCTION-VALUE
		//
		E0001_lastState = 3;
		//
		//     TO GET-FUNCTION-VALUE
		//
		status = 1;
		return;

	case 1:
		E0001_xxlo = zxlo;
		E0001_xxhi = zxhi;
		E0001_abstol = zabstl;
		E0001_reltol = zreltl;
		return;
	}

}

//****************************************************************************80***
//
//  Purpose: 
//
//    ERF_VALUES returns some values of the ERF or "error" function.
//
//  Definition:
//
//    ERF(X) = ( 2 / Math::Sqrt ( PI ) * integral ( 0 <= T <= X ) Math::Exp ( - T%2 ) dT
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//  Parameters:
//
//    Input/output, int %N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double %X, the argument of the function.
//
//    Output, double %FX, the value of the function.
//
//****************************************************************************80***
void StatClass::erf_values ( int %n_data, double %x, double %fx )
{
	const int N_MAX = 21;

	array<double> ^fx_vec = { 
		0.0000000000, 0.1124629160, 0.2227025892, 0.3286267595, 
		0.4283923550, 0.5204998778, 0.6038560908, 0.6778011938, 
		0.7421009647, 0.7969082124, 0.8427007929, 0.8802050696, 
		0.9103139782, 0.9340079449, 0.9522851198, 0.9661051465, 
		0.9763483833, 0.9837904586, 0.9890905016, 0.9927904292, 
		0.9953222650 
	};

	array<double> ^x_vec = { 
		0.0, 0.1, 0.2, 0.3, 
		0.4, 0.5, 0.6, 0.7, 
		0.8, 0.9, 1.0, 1.1, 
		1.2, 1.3, 1.4, 1.5, 
		1.6, 1.7, 1.8, 1.9, 
		2.0 
	};

	if ( n_data < 0 )
	{
		n_data = 0;
	}

	n_data = n_data + 1;

	if ( N_MAX < n_data )
	{
		n_data = 0;
		x = 0.0;
		fx = 0.0;
	}
	else
	{
		x = x_vec[n_data-1];
		fx = fx_vec[n_data-1];
	}
	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    ERROR_F evaluates the error function ERF.
//
//  Parameters:
//
//    Input, double %X, the argument.
//
//    Output, double ERROR_F, the value of the error function at X.
//
//****************************************************************************80
double StatClass::error_f ( double %x )
{
	double c = .564189583547756;
	array<double> ^a = {
		.771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
		.479137145607681e-01,.128379167095513e+00
	};
	array<double> ^b = {
		.301048631703895e-02,.538971687740286e-01,.375795757275549e+00
	};
	array<double> ^p = {
		-1.36864857382717e-07,5.64195517478974e-01,7.21175825088309e+00,
		4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
		4.51918953711873e+02,3.00459261020162e+02
	};
	array<double> ^q = {
		1.00000000000000e+00,1.27827273196294e+01,7.70001529352295e+01,
		2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
		7.90950925327898e+02,3.00459260956983e+02
	};
	array<double> ^r = {
		2.10144126479064e+00,2.62370141675169e+01,2.13688200555087e+01,
		4.65807828718470e+00,2.82094791773523e-01
	};
	array<double> ^s = {
		9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
		1.80124575948747e+01
	};

	double erf1,ax,bot,t,top,x2;

	ax = Math::Abs(x);
	if (ax <= 0.5)
	{
		t = x * x;
		top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0;
		bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0;
		erf1 = x * (top/bot);
		return erf1;
	}

	if (ax <= 4.0)
	{
		top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax + p[5]) * ax + p[6]) * ax + p[7];
		bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax + q[5]) * ax + q[6]) * ax + q[7];
		erf1 = 0.5 + (0.5-Math::Exp(-(x * x)) * top/bot);
		if(x < 0.0)
			erf1 = -erf1;
		return erf1;
	}

	if (ax < 5.8)
	{
		x2 = x * x;
		t = 1.0 / x2;
		top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
		bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0;
		erf1 = (c - top / (x2 * bot)) /ax;
		erf1 = 0.5 + (0.5 - Math::Exp(-x2) * erf1);
		if (x < 0.0)
			erf1 = -erf1;
		return erf1;
	}

	erf1 = fifdsign(1.0,x);
	return erf1;
}

//****************************************************************************80
//
//  Purpose:
// 
//    ERROR_FC evaluates the complementary error function ERFC.
//
//  Modified:
//
//    09 December 1999
//
//  Parameters:
//
//    Input, int %IND, chooses the scaling.
//    If IND is nonzero, then the value returned has been multiplied by
//    EXP(X*X).
//
//    Input, double %X, the argument of the function.
//
//    Output, double ERROR_FC, the value of the complementary 
//    error function.
//
//****************************************************************************80
double StatClass::error_fc ( int %ind, double %x )
{
	double c = .564189583547756e0;
	array<double> ^a = {
		.771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
		.479137145607681e-01,.128379167095513e+00
	};
	array<double> ^b = {
		.301048631703895e-02,.538971687740286e-01,.375795757275549e+00
	};
	array<double> ^p = {
		-1.36864857382717e-07,5.64195517478974e-01,7.21175825088309e+00,
		4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
		4.51918953711873e+02,3.00459261020162e+02
	};
	array<double> ^q = {
		1.00000000000000e+00,1.27827273196294e+01,7.70001529352295e+01,
		2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
		7.90950925327898e+02,3.00459260956983e+02
	};
	array<double> ^r = {
		2.10144126479064e+00,2.62370141675169e+01,2.13688200555087e+01,
		4.65807828718470e+00,2.82094791773523e-01
	};
	array<double> ^s = {
		9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
		1.80124575948747e+01
	};
	int K1 = 1;
	double erfc1,ax,bot,e,t,top,w;

	//
	//                     ABS(X) .LE. 0.5
	//
	ax = Math::Abs(x);
	if (ax <= 0.5e0)
	{
		t = x*x;
		top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1.0;
		bot = ((b[0]*t+b[1])*t+b[2])*t+1.0;
		erfc1 = 0.5e0+(0.5e0-x*(top/bot));
		if (ind != 0)
			erfc1 = Math::Exp(t)*erfc1;
		return erfc1;
	}

	//
	//                  0.5 .LT. ABS(X) .LE. 4
	//
	if (ax <= 4.0)
	{
		top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+p[7];
		bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax+q[7];
		erfc1 = top/bot;
	}
	else if(x <= -5.6e0)
	{

		//
		//                      ABS(X) .GT. 4
		//
		//
		//             LIMIT VALUE FOR LARGE NEGATIVE X
		//
		erfc1 = 2.0e0;
		if(ind != 0) erfc1 = 2.0e0*Math::Exp(x*x);
		return erfc1;
	}
	else if (ind != 0)
	{
			t = Math::Pow(1.0/ x,2.0);
			top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
			bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1.0;
			erfc1 = (c-t*top/bot)/ax;
	}
	else if(x > 100.0e0)
	{
		erfc1 = 0.0;
		return erfc1;
	}
	else if(x*x > -exparg(K1))
	{
		erfc1 = 0.0;
		return erfc1;
	}
	else
	{
		t = Math::Pow(1.0/ x,2.0);
		top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
		bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1.0;
		erfc1 = (c-t*top/bot)/ax;
	}
	//
	//                      FINAL ASSEMBLY
	//
	if (ind != 0)
	{
		if (x < 0.0)
			erfc1 = 2.0 * Math::Exp(x * x) - erfc1;
		return erfc1;
	}


	w = x*x;
	t = w;
	e = w-t;
	erfc1 = (0.5 +(0.5 - e)) * Math::Exp(-t) * erfc1;
	if (x < 0.0)
		erfc1 = 2.0e0-erfc1;
	return erfc1;
}

//****************************************************************************80
//
//  Purpose:
// 
//    ESUM evaluates Math::Exp ( MU + X ).
//
//  Parameters:
//
//    Input, int %MU, part of the argument.
//
//    Input, double %X, part of the argument.
//
//    Output, double ESUM, the value of Math::Exp ( MU + X ).
//
//****************************************************************************80
double StatClass::esum ( int %mu, double %x )
{
	double esum, w;

	if (x <= 0.0)
	{
		if (mu >= 0)
		{
			w = mu;
			esum = Math::Exp(w)*Math::Exp(x);
			return esum;
		}

		w = (double)mu+x;
		
		if (w <= 0.0)
		{
			w = mu;
			esum = Math::Exp(w)*Math::Exp(x);
			return esum;
		}
		
		esum = Math::Exp(w);
		return esum;
	}
	else
	{
		if (mu <= 0)
		{
			w = mu;
			esum = Math::Exp(w)*Math::Exp(x);
			return esum;
		}

		w = (double)mu+x;

		if (w >= 0.0)
		{
			w = mu;
			esum = Math::Exp(w)*Math::Exp(x);
			return esum;
		}

		esum = Math::Exp(w);
		return esum;
	}
}


//****************************************************************************80
// 
//  Purpose:
// 
//    EVAL_POL evaluates a polynomial at X.
//
//  Discussion:
//
//    EVAL_POL = A(0) + A(1)*X + ... + A(N)*X**N
//
//  Modified:
//
//    15 December 1999
//
//  Parameters:
//
//    Input, double precision A(0:N), coefficients of the polynomial.
//
//    Input, int %N, length of A.
//
//    Input, double %X, the point at which the polynomial 
//    is to be evaluated.
//
//    Output, double EVAL_POL, the value of the polynomial at X.
//
//****************************************************************************80
double StatClass::eval_pol ( array<double> ^a, int %n, double %x )
{
	double devlpl,term;
	int i;

	term = a[n-1];
	for ( i = n-1-1; i >= 0; i-- )
	{
		term = a[i]+term*x;
	}

	devlpl = term;
	return devlpl;
}

//****************************************************************************80
//
//  Purpose:
// 
//    EXPARG returns the largest or smallest legal argument for EXP.
//
//  Discussion:
//
//    Only an approximate limit for the argument of EXP is desired.
//
//  Modified:
//
//    09 December 1999
//
//  Parameters:
//
//    Input, int %L, indicates which limit is desired.
//    If L = 0, then the largest positive argument for EXP is desired.
//    Otherwise, the largest negative argument for EXP for which the
//    result is nonzero is desired.
//
//    Output, double EXPARG, the desired value.
//
//****************************************************************************80
double StatClass::exparg ( int %l )
{
	double exparg = 0.0;
	double const lnb = .69314718055995;
	double m = -1022.0;
	if (!l)
	{
		m = 1024.0;
	}
	exparg = 0.99999 * ( m * lnb);
	return exparg;
}

//****************************************************************************80***
//
//  Purpose: 
//
//    F_CDF_VALUES returns some values of the F CDF test function.
//
//  Discussion:
//
//    The value of F_CDF ( DFN, DFD, X ) can be evaluated in Mathematica by 
//    commands like:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[FRatioDistribution[ DFN, DFD ], X ]
//
//  Modified:
//
//    11 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input/output, int %N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int %A, int %B, the parameters of the function.
//
//    Output, double %X, the argument of the function.
//
//    Output, double %FX, the value of the function.
//
//****************************************************************************80***
void StatClass::f_cdf_values ( int %n_data, int %a, int %b, double %x, double %fx )
{
	const int N_MAX = 20;

	array<int> ^a_vec = { 
		1, 1, 5, 1, 
		2, 4, 1, 6, 
		8, 1, 3, 6, 
		1, 1, 1, 1,
		2, 3, 4, 5 
	};

	array<int> ^b_vec = { 
		1,  5,  1,  5, 
		10, 20,  5,  6, 
		16,  5, 10, 12, 
		5,  5,  5,  5,
		5,  5,  5,  5 
	};

	array<const double> ^fx_vec = { 
		0.500000, 0.499971, 0.499603, 0.749699, 
		0.750466, 0.751416, 0.899987, 0.899713, 
		0.900285, 0.950025, 0.950057, 0.950193, 
		0.975013, 0.990002, 0.994998, 0.999000, 
		0.568799, 0.535145, 0.514343, 0.500000 
	};

	array<const double> ^x_vec = { 
		1.00,  0.528, 1.89,  1.69, 
		1.60,  1.47,  4.06,  3.05, 
		2.09,  6.61,  3.71,  3.00, 
		10.01, 16.26, 22.78, 47.18, 
		1.00,  1.00,  1.00,  1.00 
	};

	if ( n_data < 0 )
	{
		n_data = 0;
	}

	n_data = n_data + 1;

	if ( N_MAX < n_data )
	{
		n_data = 0;
		a = 0;
		b = 0;
		x = 0.0;
		fx = 0.0;
	}
	else
	{
		a = a_vec[n_data-1];
		b = b_vec[n_data-1];
		x = x_vec[n_data-1];
		fx = fx_vec[n_data-1];
	}
	return;
}

//****************************************************************************80**
//
//  Purpose:
//
//    F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
//
//  Discussion:
//
//    The value of NONCENTRAL_F_CDF ( DFN, DFD, LAMDA, X ) can be evaluated
//    in Mathematica by commands like:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[NoncentralFRatioDistribution[ DFN, DFD, LAMBDA ], X ]
//
//  Modified:
//
//    12 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input/output, int %N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int %A, int %B, double %LAMBDA, the
//    parameters of the function.
//
//    Output, double %X, the argument of the function.
//
//    Output, double %FX, the value of the function.
//
//****************************************************************************80**
void StatClass::f_noncentral_cdf_values ( int %n_data, int %a, int %b, double %lambda, double %x, double %fx )
{
	const int N_MAX = 22;

	array<int> ^a_vec = 
	{
		1,  1,  1,  1, 
		1,  1,  1,  1, 
		1,  1,  2,  2, 
		3,  3,  4,  4, 
		5,  5,  6,  6, 
		8, 16 
	};

	array<int> ^b_vec = 
	{
		1,  5,  5,  5, 
		5,  5,  5,  5, 
		5,  5,  5, 10, 
		5,  5,  5,  5, 
		1,  5,  6, 12, 
		16,  8
	};

	array<const double> ^fx_vec = 
	{
		0.500000, 0.636783, 0.584092, 0.323443, 
		0.450119, 0.607888, 0.705928, 0.772178, 
		0.819105, 0.317035, 0.432722, 0.450270, 
		0.426188, 0.337744, 0.422911, 0.692767, 
		0.363217, 0.421005, 0.426667, 0.446402, 
		0.844589, 0.816368 
	};

	array<const double> ^lambda_vec = 
	{
		0.00,  0.000, 0.25,  1.00, 
		1.00,  1.00,  1.00,  1.00, 
		1.00,  2.00,  1.00,  1.00, 
		1.00,  2.00,  1.00,  1.00, 
		0.00,  1.00,  1.00,  1.00, 
		1.00,  1.00 
	};

	array<const double> ^x_vec = 
	{
		1.00,  1.00, 1.00,  0.50, 
		1.00,  2.00, 3.00,  4.00, 
		5.00,  1.00, 1.00,  1.00, 
		1.00,  1.00, 1.00,  2.00, 
		1.00,  1.00, 1.00,  1.00, 
		2.00,  2.00 
	};

	if ( n_data < 0 )
	{
		n_data = 0;
	}

	n_data = n_data + 1;

	if ( N_MAX < n_data )
	{
		n_data = 0;
		a = 0;
		b = 0;
		lambda = 0.0;
		x = 0.0;
		fx = 0.0;
	}
	else
	{
		a = a_vec[n_data-1];
		b = b_vec[n_data-1];
		lambda = lambda_vec[n_data-1];
		x = x_vec[n_data-1];
		fx = fx_vec[n_data-1];
	}

	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    FIFDSIGN transfers the sign of the variable "sign" to the variable "mag"
//
//  Parameters:
//
//  mag     -     magnitude
//  sign    -     sign to be transfered 
//
//****************************************************************************80
double StatClass::fifdsign ( double mag, double sign )
{
	return Math::Sign(sign) * Math::Abs(mag);
}

//****************************************************************************80
//
//  Purpose:
// 
//    FIFIDINT truncates a double number to a long integer
//
//  Parameters:
//
//  a - number to be truncated 
//
//****************************************************************************80
long StatClass::fifidint ( double a )
{
	if ( a < 1.0 ) 
	{
		return (long) 0;
	}
	else
	{ 
		return ( long ) a;
	}
}

//****************************************************************************80
//
//  Purpose:
// 
//    FIFMOD returns the modulo of a and b
//
//  Parameters:
//
//  a - numerator
//  b - denominator 
//
//****************************************************************************80
long StatClass::fifmod ( long a, long b )
{
	return ( a % b );
}

//****************************************************************************80
// 
//  Purpose:
// 
//    FPSER evaluates IX(A,B)(X) for very small B.
//
//  Discussion:
//
//    This routine is appropriate for use when 
//
//      B < min ( EPS, EPS * A ) 
//
//    and 
//
//      X <= 0.5.
//
//  Parameters:
//
//    Input, double %A, *B, parameters of the function.
//
//    Input, double %X, the point at which the function is to
//    be evaluated.
//
//    Input, double %EPS, a tolerance.
//
//    Output, double FPSER, the value of IX(A,B)(X).
//
//****************************************************************************80
double StatClass::fpser ( double %a, double %b, double %x, double %eps )
{
	int K1 = 1;
	double fpser,an,c,s,t,tol;

	fpser = 1.0;
	if (a > 1.e-3 * eps)
	{
		fpser = 0.0;
		t = a * Math::Log(x);
		if (t < exparg(K1))
			return fpser;
		fpser = Math::Exp(t);
	}
	//
	//                NOTE THAT 1/B(A,B) = B
	//
	fpser = b/ a*fpser;
	tol = eps/ a;
	an = a+1.0;
	t = x;
	s = t/an;

	do
	{
		an += 1.0;
		t = x*t;
		c = t/an;
		s += c;
	} while (Math::Abs(c) > tol);

	fpser *= (1.0 + a * s);
	return fpser;
}

//****************************************************************************80
//
//  Purpose:
// 
//    GAM1 computes 1 / GAMMA(A+1) - 1 for -0.5D+00 <= A <= 1.5
//
//  Parameters:
//
//    Input, double %A, forms the argument of the Gamma function.
//
//    Output, double GAM1, the value of 1 / GAMMA ( A + 1 ) - 1.
//
//****************************************************************************80
double StatClass::gam1 ( double %a )
{
	double s1 = .273076135303957e+00;
	double s2 = .559398236957378e-01;

	array<const double> ^p = {
		.577215664901533e+00,-.409078193005776e+00,-.230975380857675e+00,
		.597275330452234e-01,.766968181649490e-02,-.514889771323592e-02,
		.589597428611429e-03
	};

	array<const double> ^q = {
		.100000000000000e+01,.427569613095214e+00,.158451672430138e+00,
		.261132021441447e-01,.423244297896961e-02
	};

	array<const double> ^r = {
		-.422784335098468e+00,-.771330383816272e+00,-.244757765222226e+00,
		.118378989872749e+00,.930357293360349e-03,-.118290993445146e-01,
		.223047661158249e-02,.266505979058923e-03,-.132674909766242e-03
	};
	double gam1,bot,d,t,top,w,T1;

	if ((a < -0.5) || (a > 1.5))
		throw gcnew ArgumentException("a must be between -0.5 and 1.5", "a");

	t = a;
	d = a-0.5e0;
	if(d > 0.0)
		t = d-0.5e0;
	T1 = t;
	if (T1 < 0)
	{
		top = (((((((r[8]*t+r[7])*t+r[6])*t+r[5])*t+r[4])*t+r[3])*t+r[2])*t+r[1])*t+r[0];
		bot = (s2*t+s1)*t+1.0;
		w = top/bot;
		if (d <= 0.0)
		{
			gam1 = a*(w+0.5e0+0.5e0);
			return gam1;
		}
		else
		{
			gam1 = t*w/ a;
			return gam1;
		}
	}
	else if(T1 == 0)
	{
		gam1 = 0.0;
		return gam1;
	}
	else
	{
		top = (((((p[6]*t+p[5])*t+p[4])*t+p[3])*t+p[2])*t+p[1])*t+p[0];
		bot = (((q[4]*t+q[3])*t+q[2])*t+q[1])*t+1.0;
		w = top/bot;
		if (d <= 0.0)
		{
			gam1 = a*w;
			return gam1;
		}
		else
		{
			gam1 = t/ a*(w-0.5e0-0.5e0);
			return gam1;
		}
	}
}

//****************************************************************************80
//
//  Purpose:
// 
//    GAMMA_INC evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
//
//  Discussion:
//
//    This is certified spaghetti code.
//
//  Author:
//
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//
//  Parameters:
//
//    Input, double %A, *X, the arguments of the incomplete
//    gamma ratio.  A and X must be nonnegative.  A and X cannot
//    both be zero.
//
//    Output, double %ANS, *QANS.  On normal output,
//    ANS = P(A,X) and QANS = Q(A,X).  However, ANS is set to 2 if
//    A or X is negative, or both are 0, or when the answer is
//    computationally indeterminate because A is extremely large
//    and X is very close to A.
//
//    Input, int %IND, indicates the accuracy request:
//    0, as much accuracy as possible.
//    1, to within 1 unit of the 6-th significant digit, 
//    otherwise, to within 1 unit of the 3rd significant digit.
//
//****************************************************************************80
void StatClass::gamma_inc ( double %a, double %x, double %ans, double %qans, int %ind )
{
	double alog10 = 2.30258509299405e0;
	double d10 = -.185185185185185e-02;
	double d20 = .413359788359788e-02;
	double d30 = .649434156378601e-03;
	double d40 = -.861888290916712e-03;
	double d50 = -.336798553366358e-03;
	double d60 = .531307936463992e-03;
	double d70 = .344367606892378e-03;
	double rt2pin = .398942280401433e0;
	double rtpi = 1.77245385090552e0;
	double third = .333333333333333e0;
	array<const double> ^acc0 = {
		5.e-15,5.e-7,5.e-4
	};
	array<const double> ^big = {
		20.0e0,14.0e0,10.0e0
	};
	array<const double> ^d0 = {
		.833333333333333e-01,-.148148148148148e-01,.115740740740741e-02,
		.352733686067019e-03,-.178755144032922e-03,.391926317852244e-04,
		-.218544851067999e-05,-.185406221071516e-05,.829671134095309e-06,
		-.176659527368261e-06,.670785354340150e-08,.102618097842403e-07,
		-.438203601845335e-08
	};
	array<const double> ^d1 = {
		-.347222222222222e-02,.264550264550265e-02,-.990226337448560e-03,
		.205761316872428e-03,-.401877572016461e-06,-.180985503344900e-04,
		.764916091608111e-05,-.161209008945634e-05,.464712780280743e-08,
		.137863344691572e-06,-.575254560351770e-07,.119516285997781e-07
	};
	array<const double> ^d2 = {
		-.268132716049383e-02,.771604938271605e-03,.200938786008230e-05,
		-.107366532263652e-03,.529234488291201e-04,-.127606351886187e-04,
		.342357873409614e-07,.137219573090629e-05,-.629899213838006e-06,
		.142806142060642e-06
	};
	array<const double> ^d3 = {
		.229472093621399e-03,-.469189494395256e-03,.267720632062839e-03,
		-.756180167188398e-04,-.239650511386730e-06,.110826541153473e-04,
		-.567495282699160e-05,.142309007324359e-05
	};
	array<const double> ^d4 = {
		.784039221720067e-03,-.299072480303190e-03,-.146384525788434e-05,
		.664149821546512e-04,-.396836504717943e-04,.113757269706784e-04
	};
	array<const double> ^d5 = {
		-.697281375836586e-04,.277275324495939e-03,-.199325705161888e-03,
		.679778047793721e-04
	};
	array<const double> ^d6 = {
		-.592166437353694e-03,.270878209671804e-03
	};
	array<const double> ^e00 = {
		.25e-3,.25e-1,.14e0
	};
	array<const double> ^x00 = {
		31.0e0,17.0e0,9.7e0
	};
	int K1 = 1;
	int K2 = 0;
	double a2n,a2nm1,acc,am0,amn,an,an0,apn,b2n,b2nm1,c,c0,c1,c2,c3,c4,c5,c6,
		cma,e,e0,g,h,j,l,r,rta,rtx,s,sum,t,t1,tol,twoa,u,w,x0,y,z;
	int i,iop,m,max,n;
	double T3;
	array<double> ^wk = gcnew array<double>(20);
	int T4,T5;
	double T6,T7;

	//
	//  E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
	//  NUMBER FOR WHICH 1.0 + E .GT. 1.0 .
	//
	e = 2.22044604925031E-16;
	if(a < 0.0 || x < 0.0)
	{
		//
		//  ERROR RETURN
		//
		ans = 2.0e0;
		return;
	}

	if(a == 0.0 && x == 0.0)
	{
		//
		//  ERROR RETURN
		//
		ans = 2.0e0;
		return;
	}
	if(a*x == 0.0)
	{
		if (x <= a)
		{
			ans = 0.0;
			qans = 1.0;
			return;
		}
		ans = 1.0;
		qans = 0.0;
		return;
	}
	iop = ind+1;
	if(iop != 1 && iop != 2) iop = 3;
	acc = Math::Max(acc0[iop-1],e);
	e0 = e00[iop-1];
	x0 = x00[iop-1];

	//
	//  SELECT THE APPROPRIATE ALGORITHM
	//
	if (a < 1.0)
	{
		if (a == 0.5e0)
		{
			if(x >= 0.25e0)
			{
				T7 = Math::Sqrt(x);
				qans = error_fc ( K2, T7 );
				ans = 0.5e0+(0.5e0-qans);
				return;
			}
			T6 = Math::Sqrt(x);
			ans = error_f ( T6 );
			qans = 0.5e0+(0.5e0-ans);
			return;
		}

		if(x < 1.1e0)
		{
			//
			//  TAYLOR SERIES FOR P(A,X)/X**A
			//
			an = 3.0e0;
			c = x;
			sum = x/(a+3.0e0);
			tol = 3.0e0*acc/(a+1.0);

			do
			{
				an += 1.0;
				c = -(c*(x/an));
				t = c/(a+an);
				sum += t;
			} while (Math::Abs(t) > tol);

			j = a*x*((sum/6.0e0-0.5e0/(a+2.0e0))*x+1.0/(a+1.0));
			z = a*Math::Log(x);
			h = gam1(a);
			g = 1.0+h;
			if (x >= 0.25e0)
			{
				if(a < x/2.59e0)
				{
					l = rexp(z);
					w = 0.5e0+(0.5e0+l);
					qans = (w*j-l)*g-h;
					if (qans < 0.0)
					{
						ans = 1.0;
						qans = 0.0;
						return;
					}
					ans = 0.5e0+(0.5e0-qans);
					return;
				}
			}
			else
			{
				if(z > -.13394e0)
				{
					l = rexp(z);
					w = 0.5e0+(0.5e0+l);
					qans = (w*j-l)*g-h;
					if (qans < 0.0)
					{
						ans = 1.0;
						qans = 0.0;
						return;
					}
					ans = 0.5e0+(0.5e0-qans);
					return;
				}
			}

			w = Math::Exp(z);
			ans = w*g*(0.5e0+(0.5e0-j));
			qans = 0.5e0+(0.5e0-ans);
			return;
		}

		t1 = a*Math::Log(x)-x;
		u = a*Math::Exp(t1);
		if (u == 0.0)
		{
			ans = 1.0;
			qans = 0.0;
			return;
		}
		r = u*(1.0+gam1(a));
		//
		//  CONTINUED FRACTION EXPANSION
		//
		tol = Math::Max(5.0e0*e,acc);
		a2nm1 = a2n = 1.0;
		b2nm1 = x;
		b2n = x+(1.0-a);
		c = 1.0;

		do
		{
			a2nm1 = x*a2n+c*a2nm1;
			b2nm1 = x*b2n+c*b2nm1;
			am0 = a2nm1/b2nm1;
			c += 1.0;
			cma = c-a;
			a2n = a2nm1+cma*a2n;
			b2n = b2nm1+cma*b2n;
			an0 = a2n/b2n;
		} while ((Math::Abs(an0-am0) >= tol*an0));

		qans = r*an0;
		ans = 0.5e0+(0.5e0-qans);
		return;
	}

	if (a < big[iop-1])
	{
		if (a <= x && x < x0)
		{
			twoa = a+a;
			m = fifidint(twoa);
			if (twoa == (double)m)
			{
				i = m/2;
				if(a == (double)i)
				{
					// 
					//  FINITE SUMS FOR Q WHEN A .GE. 1 AND 2*A IS AN INTEGER
					//
					sum = Math::Exp(-x);
					t = sum;
					n = 1;
					c = 0.0;
					while (n != i)
					{
						n += 1;
						c += 1.0;
						t = x*t/c;
						sum += t;
					}
					ans = 1.0 - sum;
					qans = sum;
					return;
				}
				rtx = Math::Sqrt(x);
				sum = error_fc ( K2, rtx );
				t = Math::Exp(-x)/(rtpi*rtx);
				n = 0;
				c = -0.5e0;
			}
		}

		t1 = a*Math::Log(x)-x;
		r = Math::Exp(t1)/ gamma_x(a);

	}

	else
	{
		l = x/ a;
		if(l == 0.0)
		{
			ans = 0.0;
			qans = 1.0;
			return;
		}
		s = 0.5e0+(0.5e0-l);
		z = rlog(l);
		if (z >= 700.0e0/ a)
		{
			if(Math::Abs(s) <= 2.0e0*e)
			{
				//
				//  ERROR RETURN
				//
				ans = 2.0e0;
				return;
			}

			if(x <= a)
			{
				ans = 0.0;
				qans = 1.0;
				return;
			}
			ans = 1.0;
			qans = 0.0;
			return;
		}
		y = a*z;
		rta = Math::Sqrt(a);
		if (Math::Abs(s) <= e0/rta)
		{
			//
			//  TEMME EXPANSION FOR L = 1
			//
			if(a*e*e > 3.28e-3)
			{
				//
				//  ERROR RETURN
				//
				ans = 2.0e0;
				return;
			}
			c = 0.5e0+(0.5e0-y);
			w = (0.5e0-Math::Sqrt(y)*(0.5e0+(0.5e0-y/3.0e0))/rtpi)/c;
			u = 1.0/ a;
			z = Math::Sqrt(z+z);
			if(l < 1.0) z = -z;
			T5 = iop-2;
			if(T5 < 0)
			{
				c0 = ((((((d0[6]*z+d0[5])*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-
					third;
				c1 = (((((d1[5]*z+d1[4])*z+d1[3])*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
				c2 = ((((d2[4]*z+d2[3])*z+d2[2])*z+d2[1])*z+d2[0])*z+d20;
				c3 = (((d3[3]*z+d3[2])*z+d3[1])*z+d3[0])*z+d30;
				c4 = (d4[1]*z+d4[0])*z+d40;
				c5 = (d5[1]*z+d5[0])*z+d50;
				c6 = d6[0]*z+d60;
				t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u+c0;

				if (l >= 1.0)
				{
					qans = c*(w+rt2pin*t/rta);
					ans = 0.5e0+(0.5e0-qans);
					return;
				}
				else
				{
					ans = c*(w-rt2pin*t/rta);
					qans = 0.5e0+(0.5e0-ans);
					return;
				}
			}
			else if(T5 == 0)
			{

				c0 = (d0[1]*z+d0[0])*z-third;
				c1 = d1[0]*z+d10;
				t = (d20*u+c1)*u+c0;

				if (l >= 1.0)
				{
					qans = c*(w+rt2pin*t/rta);
					ans = 0.5e0+(0.5e0-qans);
					return;
				}
				else
				{
					ans = c*(w-rt2pin*t/rta);
					qans = 0.5e0+(0.5e0-ans);
					return;
				}
			}

			t = d0[0]*z-third;

			if (l >= 1.0)
			{
				qans = c*(w+rt2pin*t/rta);
				ans = 0.5e0+(0.5e0-qans);
				return;
			}
			else
			{
				ans = c*(w-rt2pin*t/rta);
				qans = 0.5e0+(0.5e0-ans);
				return;
			}
		}
		if (Math::Abs(s) <= 0.4e0)
		{
			//
			//  GENERAL TEMME EXPANSION
			//
			if(Math::Abs(s) <= 2.0e0*e && a*e*e > 3.28e-3)
			{
				//
				//  ERROR RETURN
				//
				ans = 2.0e0;
				return;
			}

			c = Math::Exp(-y);
			T3 = Math::Sqrt(y);
			w = 0.5e0 * error_fc ( K1, T3 );
			u = 1.0/ a;
			z = Math::Sqrt(z+z);
			if(l < 1.0) z = -z;
			T4 = iop-2;
			if (T4 < 0)
			{
				if(Math::Abs(s) <= 1.e-3)
				{
					c0 = ((((((d0[6]*z+d0[5])*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-
						third;
					c1 = (((((d1[5]*z+d1[4])*z+d1[3])*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
					c2 = ((((d2[4]*z+d2[3])*z+d2[2])*z+d2[1])*z+d2[0])*z+d20;
					c3 = (((d3[3]*z+d3[2])*z+d3[1])*z+d3[0])*z+d30;
					c4 = (d4[1]*z+d4[0])*z+d40;
					c5 = (d5[1]*z+d5[0])*z+d50;
					c6 = d6[0]*z+d60;
					t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u+c0;

					if (l >= 1.0)
					{
						qans = c*(w+rt2pin*t/rta);
						ans = 0.5e0+(0.5e0-qans);
						return;
					}
					else
					{
						ans = c*(w-rt2pin*t/rta);
						qans = 0.5e0+(0.5e0-ans);
						return;
					}
				}
				c0 = ((((((((((((d0[12]*z+d0[11])*z+d0[10])*z+d0[9])*z+d0[8])*z+d0[7])*z+d0[
					6])*z+d0[5])*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-third;
					c1 = (((((((((((d1[11]*z+d1[10])*z+d1[9])*z+d1[8])*z+d1[7])*z+d1[6])*z+d1[5]
					)*z+d1[4])*z+d1[3])*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
					c2 = (((((((((d2[9]*z+d2[8])*z+d2[7])*z+d2[6])*z+d2[5])*z+d2[4])*z+d2[3])*z+
						d2[2])*z+d2[1])*z+d2[0])*z+d20;
					c3 = (((((((d3[7]*z+d3[6])*z+d3[5])*z+d3[4])*z+d3[3])*z+d3[2])*z+d3[1])*z+
						d3[0])*z+d30;
					c4 = (((((d4[5]*z+d4[4])*z+d4[3])*z+d4[2])*z+d4[1])*z+d4[0])*z+d40;
					c5 = (((d5[3]*z+d5[2])*z+d5[1])*z+d5[0])*z+d50;
					c6 = (d6[1]*z+d6[0])*z+d60;
					t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u+c0;
					if (l >= 1.0)
					{
						qans = c*(w+rt2pin*t/rta);
						ans = 0.5e0+(0.5e0-qans);
						return;
					}
					else
					{
						ans = c*(w-rt2pin*t/rta);
						qans = 0.5e0+(0.5e0-ans);
						return;
					}
			}
			else if(T4 == 0)
			{
				c0 = (((((d0[5]*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-third;
				c1 = (((d1[3]*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
				c2 = d2[0]*z+d20;
				t = (c2*u+c1)*u+c0;
				if (l >= 1.0)
				{
					qans = c*(w+rt2pin*t/rta);
					ans = 0.5e0+(0.5e0-qans);
					return;
				}
				else
				{
					ans = c*(w-rt2pin*t/rta);
					qans = 0.5e0+(0.5e0-ans);
					return;
				}
			}
			else
			{
				t = ((d0[2]*z+d0[1])*z+d0[0])*z-third;
				if (l >= 1.0)
				{
					qans = c*(w+rt2pin*t/rta);
					ans = 0.5e0+(0.5e0-qans);
					return;
				}
				else
				{
					ans = c*(w-rt2pin*t/rta);
					qans = 0.5e0+(0.5e0-ans);
					return;
				}
			}
		}
		t = Math::Pow(1.0/ a,2.0);
		t1 = (((0.75e0*t-1.0)*t+3.5e0)*t-105.0e0)/(a*1260.0e0);
		t1 -= y;
		r = rt2pin*rta*Math::Exp(t1);
	}

	if (r == 0.0)
	{
		if(x <= a)
		{
			ans = 0.0;
			qans = 1.0;
			return;
		}
		ans = 1.0;
		qans = 0.0;
		return;
	}

	if (x > Math::Max(a,alog10))
	{
		if(x < x0)
		{
			//
			//  CONTINUED FRACTION EXPANSION
			//
			tol = Math::Max(5.0e0*e,acc);
			a2nm1 = a2n = 1.0;
			b2nm1 = x;
			b2n = x+(1.0-a);
			c = 1.0;

			do
			{
				a2nm1 = x*a2n+c*a2nm1;
				b2nm1 = x*b2n+c*b2nm1;
				am0 = a2nm1/b2nm1;
				c += 1.0;
				cma = c-a;
				a2n = a2nm1+cma*a2n;
				b2n = b2nm1+cma*b2n;
				an0 = a2n/b2n;
			} while ((Math::Abs(an0-am0) >= tol*an0));

			qans = r*an0;
			ans = 0.5e0+(0.5e0-qans);
			return;
		}

		//
		//  ASYMPTOTIC EXPANSION
		//
		amn = a-1.0;
		t = amn/ x;
		wk[0] = t;

		for ( n = 2; n <= 20; n++ )
		{
			amn -= 1.0;
			t *= (amn/ x);
			if (Math::Abs(t) <= 1.e-3) break;
			wk[n-1] = t;
		}

		if (Math::Abs(t) > 1.e-3)
			n = 20;

		sum = t;

		while (Math::Abs(t) > acc)
		{
			amn -= 1.0;
			t *= (amn/ x);
			sum += t;
		}

		max = n-1;
		for ( m = 1; m <= max; m++ )
		{
			n -= 1;
			sum += wk[n-1];
		}
		qans = r/ x*(1.0+sum);
		ans = 0.5e0+(0.5e0-qans);
		return;
	}

	//
	//  TAYLOR SERIES FOR P/R
	//
	apn = a+1.0;
	t = x/apn;
	wk[0] = t;
	for ( n = 2; n <= 20; n++ )
	{
		apn += 1.0;
		t *= (x/apn);
		if (t <= 1.e-3) break;
		wk[n-1] = t;
	}
	if (t > 1.e-3)
		n = 20;

	sum = t;
	tol = 0.5e0*acc;

	do
	{
		apn += 1.0;
		t *= (x/apn);
		sum += t;
	} while (t > tol);
	max = n-1;
	for ( m = 1; m <= max; m++ )
	{
		n -= 1;
		sum += wk[n-1];
	}
	ans = r/ a*(1.0+sum);
	qans = 0.5e0+(0.5e0-ans);
	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    GAMMA_INC_INV computes the inverse incomplete gamma ratio function.
//
//  Discussion:
//
//    The routine is given positive A, and nonnegative P and Q where P + Q = 1.
//    The value X is computed with the property that P(A,X) = P and Q(A,X) = Q.  
//    Schroder iteration is employed.  The routine attempts to compute X
//    to 10 significant digits if this is possible for the particular computer 
//    arithmetic being used.
//
//  Author:
//
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//
//  Parameters:
//
//    Input, double %A, the parameter in the incomplete gamma
//    ratio.  A must be positive.
//
//    Output, double %X, the computed point for which the
//    incomplete gamma functions have the values P and Q.
//
//    Input, double %X0, an optional initial approximation
//    for the solution X.  If the user does not want to supply an
//    initial approximation, then X0 should be set to 0, or a negative
//    value.
//
//    Input, double %P, *Q, the values of the incomplete gamma
//    functions, for which the corresponding argument is desired.
//
//    Output, int %IERR, error flag.
//    0, the solution was obtained. Iteration was not used.
//    0 < K, The solution was obtained. IERR iterations were performed.
//    -2, A <= 0
//    -3, No solution was obtained. The ratio Q/A is too large.
//    -4, P + Q /= 1
//    -6, 20 iterations were performed. The most recent value obtained 
//        for X is given.  This cannot occur if X0 <= 0.
//    -7, Iteration failed. No value is given for X.
//        This may occur when X is approximately 0.
//    -8, A value for X has been obtained, but the routine is not certain
//        of its accuracy.  Iteration cannot be performed in this
//        case. If X0 <= 0, this can occur only when P or Q is 
//        approximately 0. If X0 is positive then this can occur when A is
//        exceedingly close to X and A is extremely large (say A .GE. 1.E20).
//
//****************************************************************************80
void StatClass::gamma_inc_inv ( double %a, double %x, double %x0, double %p, double %q, int %ierr )
{
	double a0 = 3.31125922108741e0;
	double a1 = 11.6616720288968e0;
	double a2 = 4.28342155967104e0;
	double a3 = .213623493715853e0;
	double b1 = 6.61053765625462e0;
	double b2 = 6.40691597760039e0;
	double b3 = 1.27364489782223e0;
	double b4 = .036117081018842e0;
	double c = .577215664901533e0;
	double ln10 = 2.302585e0;
	double tol = 1.e-12;
	array<const double> ^amin = {
		500.0e0,100.0e0
	};
	array<const double> ^bmin = {
		1.e-28,1.e-13
	};
	array<const double> ^dmin = {
		1.e-06,1.e-04
	};
	array<const double> ^emin = {
		2.e-03,6.e-03
	};
	array<const double> ^eps0 = {
		1.e-10,1.e-08
	};
	int K8 = 0;
	double am1,amax,ap1,ap2,ap3,apn,b = 0,c1,c2,c3,c4,c5,d,e,e2,eps,g,h,pn,qg,qn,
		r,rta,s,s2,sum,t,u,w,xmax,xmin,xn,y,z;
	int iop;
	double T4,T5,T6,T7,T9;

	//
	//  E, XMIN, AND XMAX ARE MACHINE DEPENDENT CONSTANTS.
	//            E IS THE SMALLEST NUMBER FOR WHICH 1.0 + E .GT. 1.0.
	//            XMIN IS THE SMALLEST POSITIVE NUMBER AND XMAX IS THE
	//            LARGEST POSITIVE NUMBER.
	//
	e = 2.22044604925031E-16;
	xmin = 4.94065645841247e-324;
	xmax = double::MaxValue;
	x = 0.0;
	if(a <= 0.0)
	{
		ierr = -2;
		return;
	}
	t = p + q - 1.0;
	if (Math::Abs(t) > e)
	{
		ierr = -4;
		return;
	}
	ierr = 0;
	if (p == 0.0) return;
	if (q == 0.0)
	{
		//
		//                       SPECIAL CASES
		//
		x = xmax;
		return;
	}
	if (a == 1.0)
	{
		if (q >= 0.9e0)
		{
			T9 = -p;
			x = -alnrel(T9);
		}
		else
		{
			x = -Math::Log(q);
		}
		return;
	}
	e2 = 2.0e0*e;
	amax = 0.4e-10/(e*e);
	iop = 1;
	if(e > 1.e-10) iop = 2;
	eps = eps0[iop-1];
	xn = x0;
	if (x0 > 0.0)
	{
		if (p > 0.5e0) goto ShroderIterQ;
		goto ShroderIterP;
	}
	//
	//        SELECTION OF THE INITIAL APPROXIMATION XN OF X
	//                       WHEN A .LT. 1
	//
	if (a > 1.0)
	{
		//
		//        SELECTION OF THE INITIAL APPROXIMATION XN OF X
		//                       WHEN A .GT. 1
		//
		if (q > 0.5e0)
		{
			w = Math::Log(p);
		}
		else
		{
			w = Math::Log(q);
		}

		t = Math::Sqrt(-(2.0e0*w));
		s = t-(((a3*t+a2)*t+a1)*t+a0)/((((b4*t+b3)*t+b2)*t+b1)*t+1.0);
		if(q > 0.5e0) s = -s;
		rta = Math::Sqrt(a);
		s2 = s*s;
		xn = a+s*rta+(s2-1.0)/3.0e0+s*(s2-7.0e0)/(36.0e0*rta)-((3.0e0*s2+7.0e0)*
			s2-16.0e0)/(810.0e0*a)+s*((9.0e0*s2+256.0e0)*s2-433.0e0)/(38880.0e0*a*
			rta);
		xn = Math::Max(xn,0.0);
		if (a >= amin[iop-1])
		{
			x = xn;
			d = 0.5e0+(0.5e0-x/ a);
			if(Math::Abs(d) <= dmin[iop-1]) return;
		}

		if (p > 0.5e0)
		{
			if(xn < 3.0e0*a) goto ShroderIterQ;
			y = -(w+ gamma_log ( a ) );
			d = Math::Max(2.0e0,a*(a-1.0));
			if (y < ln10*d)
			{
				s = 1.0-a;
				z = Math::Log(y);
				c1 = -(s*z);
				c2 = -(s*(1.0+c1));
				c3 = s*((0.5e0*c1+(2.0e0-a))*c1+(2.5e0-1.5e0*a));
				c4 = -(s*(((c1/3.0e0+(2.5e0-1.5e0*a))*c1+((a-6.0e0)*a+7.0e0))*c1+(
					(11.0e0*a-46.0)*a+47.0e0)/6.0e0));
				c5 = -(s*((((-(c1/4.0e0)+(11.0e0*a-17.0e0)/6.0e0)*c1+((-(3.0e0*a)+13.0e0)*
					a-13.0e0))*c1+0.5e0*(((2.0e0*a-25.0e0)*a+72.0e0)*a-61.0e0))*c1+((
					(25.0e0*a-195.0e0)*a+477.0e0)*a-379.0e0)/12.0e0));
				xn = (((c5/y+c4)/y+c3)/y+c2)/y+c1+y;
				if(a > 1.0) goto ShroderIterQ;
				if(b > bmin[iop-1]) goto ShroderIterQ;
				x = xn;
				return;
			}

			t = a-1.0;
			T6 = -(t/(xn+1.0));
			xn = y+t*Math::Log(xn)-alnrel(T6);
			T7 = -(t/(xn+1.0));
			xn = y+t*Math::Log(xn)-alnrel(T7);
			goto ShroderIterQ;
		}

		ap1 = a+1.0;
		if(xn > 0.70e0*ap1) goto ShroderIterP;
		w += gamma_log ( ap1 );
		
		if (xn > 0.15e0*ap1)
		{

			apn = ap1;
			t = xn/apn;
			sum = 1.0+t;

			do
			{
				apn += 1.0;
				t *= (xn/apn);
				sum += t;
			} while (t > 1.e-4);
			t = w-Math::Log(sum);
			xn = Math::Exp((xn+t)/ a);
			xn *= (1.0-(a*Math::Log(xn)-xn-t)/(a-xn));
			goto ShroderIterP;
		}

		ap2 = a+2.0e0;
		ap3 = a+3.0e0;
		x = Math::Exp((w+x)/ a);
		x = Math::Exp((w+x-Math::Log(1.0+x/ap1*(1.0+x/ap2)))/ a);
		x = Math::Exp((w+x-Math::Log(1.0+x/ap1*(1.0+x/ap2)))/ a);
		x = Math::Exp((w+x-Math::Log(1.0+x/ap1*(1.0+x/ap2*(1.0+x/ap3))))/ a);
		xn = x;
		if (xn > 1.e-2*ap1)
		{

			apn = ap1;
			t = xn/apn;
			sum = 1.0+t;

			do
			{
				apn += 1.0;
				t *= (xn/apn);
				sum += t;
			} while (t > 1.e-4);
			t = w-Math::Log(sum);
			xn = Math::Exp((xn+t)/ a);
			xn *= (1.0-(a*Math::Log(xn)-xn-t)/(a-xn));
			goto ShroderIterP;
		}
		if(xn <= emin[iop-1]*ap1) return;
		goto ShroderIterP;
	}

	T4 = a+1.0;
	g = gamma_x(T4);
	qg = q*g;
	if (qg == 0.0)
	{
		x = xmax;
		ierr = -8;
		return;
	}
	b = qg/ a;

	if (qg > 0.6e0*a)
	{
		if (b*q <= 1.e-8)
		{
			if (p > 0.9e0)
			{
				T5 = -q;
				xn = Math::Exp((alnrel(T5)+ gamma_ln1 ( a ) ) / a );
			}
			else
			{
				xn = Math::Exp(Math::Log(p*g)/ a);
			}

			if (xn == 0.0)
			{
				ierr = -3;
				return;
			}
			t = 0.5e0+(0.5e0-xn/(a+1.0));
			xn /= t;
			if(p > 0.5e0) goto ShroderIterQ;
			goto ShroderIterP;
		}
		xn = Math::Exp(-(q/ a+c));

		if (xn == 0.0)
		{
			ierr = -3;
			return;
		}
		t = 0.5e0+(0.5e0-xn/(a+1.0));
		xn /= t;
		if(p > 0.5e0) goto ShroderIterQ;
		goto ShroderIterP;
	}

	if (a < 0.30e0 && b >= 0.35e0)
	{
		t = Math::Exp(-(b+c));
		u = t*Math::Exp(t);
		xn = t*Math::Exp(u);

		if(p > 0.5e0) goto ShroderIterQ;
		goto ShroderIterP;
	}

	if (b >= 0.45e0)
	{
		if (b*q <= 1.e-8)
		{
			if (p > 0.9e0)
			{
				T5 = -q;
				xn = Math::Exp((alnrel(T5)+ gamma_ln1 ( a ) ) / a );
			}
			else
			{
				xn = Math::Exp(Math::Log(p*g)/ a);
			}

			if (xn == 0.0)
			{
				ierr = -3;
				return;
			}
			t = 0.5e0+(0.5e0-xn/(a+1.0));
			xn /= t;
			if(p > 0.5e0) goto ShroderIterQ;
			goto ShroderIterP;
		}
		xn = Math::Exp(-(q/ a+c));

		if (xn == 0.0)
		{
			ierr = -3;
			return;
		}
		t = 0.5e0+(0.5e0-xn/(a+1.0));
		xn /= t;
		if(p > 0.5e0) goto ShroderIterQ;
		goto ShroderIterP;
	}

	if(b == 0.0)
	{
		x = xmax;
		ierr = -8;
		return;
	}

	y = -Math::Log(b);
	s = 0.5e0+(0.5e0-a);
	z = Math::Log(y);
	t = y-s*z;

	if (b >= 0.15e0)
	{
		xn = y-s*Math::Log(t)-Math::Log(1.0+s/(t+1.0));
		goto ShroderIterQ;
	}

	if (b > 0.01e0)
	{
		u = ((t+2.0e0*(3.0e0-a))*t+(2.0e0-a)*(3.0e0-a))/((t+(5.0e0-a))*t+2.0e0);
		xn = y-s*Math::Log(t)-Math::Log(u);
		goto ShroderIterQ;
	}

ShroderIterP:
	if (p <= 1.e10*xmin)
	{
		x = xn;
		ierr = -8;
		return;
	}
	am1 = a-0.5e0-0.5e0;

	do
	{
		if (a > amax)
		{
			d = 0.5e0+(0.5e0-xn/ a);
			if (Math::Abs(d) <= e2)
			{
				x = xn;
				ierr = -8;
				return;
			}
		}
		if (ierr >= 20)
		{
			ierr = -6;
			return;
		}
		ierr += 1;
		gamma_inc ( a, xn, pn, qn, K8 );
		if (pn == 0.0 || qn == 0.0)
		{
			x = xn;
			ierr = -8;
			return;
		}
		r = rcomp(a,xn);
		if (r == 0.0)
		{
			x = xn;
			ierr = -8;
			return;
		}
		t = (pn-p)/r;
		w = 0.5e0*(am1-xn);
		if (Math::Abs(t) <= 0.1e0 && Math::Abs(w*t) <= 0.1e0)
		{
			h = t*(1.0+w*t);
			x = xn*(1.0-h);
			if (x <= 0.0)
			{
				ierr = -7;
				return;
			}
			if(Math::Abs(w) >= 1.0 && Math::Abs(w)*t*t <= eps) return;
			d = Math::Abs(h);
			xn = x;
			if (d > tol) continue;
			if(d <= eps) return;
			if(Math::Abs(p-pn) <= tol*p) return;
		}
		else
		{
			x = xn*(1.0-t);
			if (x <= 0.0)
			{
				ierr = -7;
				return;
			}
			d = Math::Abs(t);
			xn = x;
			if (d > tol) continue;
			if (d <= eps) return;
			if(Math::Abs(p-pn) <= tol*p) return;
		}
	} while (true);

ShroderIterQ:
	//
	//                 SCHRODER ITERATION USING Q
	//
	if (q <= 1.e10*xmin)
	{
		x = xn;
		ierr = -8;
		return;
	}
	am1 = a-0.5e0-0.5e0;
	do
	{
		if (a > amax)
		{
			d = 0.5e0+(0.5e0-xn/ a);
			if (Math::Abs(d) <= e2)
			{
				x = xn;
				ierr = -8;
				return;
			}
		}
		if (ierr >= 20)
		{
			ierr = -6;
			return;
		}
		ierr += 1;
		gamma_inc ( a, xn, pn, qn, K8 );
		if (pn == 0.0 || qn == 0.0)
		{
			x = xn;
			ierr = -8;
			return;
		}
		r = rcomp(a,xn);
		if (r == 0.0)
		{
			x = xn;
			ierr = -8;
			return;
		}
		t = (q-qn)/r;
		w = 0.5e0*(am1-xn);
		if (Math::Abs(t) <= 0.1e0 && Math::Abs(w*t) <= 0.1e0)
		{
			h = t*(1.0+w*t);
			x = xn*(1.0-h);
			if (x <= 0.0)
			{
				ierr = -7;
				return;
			}
			if(Math::Abs(w) >= 1.0 && Math::Abs(w)*t*t <= eps) return;
			d = Math::Abs(h);
		}
		x = xn*(1.0-t);
		if (x <= 0.0)
		{
			ierr = -7;
			return;
		}
		d = Math::Abs(t);
		xn = x;
		
		if (d > tol) continue;
		if(d <= eps) break;
	} while (Math::Abs(q-qn) > tol*q);

}

//****************************************************************************80***
//
//  Purpose: 
//
//    GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
//
//  Discussion:
//
//    The (normalized) incomplete Gamma function P(A,X) is defined as:
//
//      PN(A,X) = 1/GAMMA(A) * Integral ( 0 <= T <= X ) T**(A-1) * Math::Exp(-T) dT.
//
//    With this definition, for all A and X,
//
//      0 <= PN(A,X) <= 1
//
//    and
//
//      PN(A,INFINITY) = 1.0
//
//    Mathematica can compute this value as
//
//      1 - GammaRegularized[A,X]
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//  Parameters:
//
//    Input/output, int %N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double %A, the parameter of the function.
//
//    Output, double %X, the argument of the function.
//
//    Output, double %FX, the value of the function.
//
//****************************************************************************80***
void StatClass::gamma_inc_values ( int %n_data, double %a, double %x, double %fx )
{
	const int N_MAX = 20;

	array<const double> ^a_vec =
	{ 
		0.1,  0.1,  0.1,  0.5, 
		0.5,  0.5,  1.0,  1.0, 
		1.0,  1.1,  1.1,  1.1, 
		2.0,  2.0,  2.0,  6.0, 
		6.0, 11.0, 26.0, 41.0
	};
	array<const double> ^fx_vec = 
	{ 
		0.7420263, 0.9119753, 0.9898955, 0.2931279, 
		0.7656418, 0.9921661, 0.0951626, 0.6321206, 
		0.9932621, 0.0757471, 0.6076457, 0.9933425, 
		0.0091054, 0.4130643, 0.9931450, 0.0387318, 
		0.9825937, 0.9404267, 0.4863866, 0.7359709 
	};
	array<const double> ^x_vec = 
	{ 
		3.1622777E-02, 3.1622777E-01, 1.5811388, 7.0710678E-02, 
		7.0710678E-01, 3.5355339, 0.1000000, 1.0000000, 
		5.0000000, 1.0488088E-01, 1.0488088, 5.2440442, 
		1.4142136E-01, 1.4142136, 7.0710678, 2.4494897, 
		1.2247449E+01, 1.6583124E+01, 2.5495098E+01, 4.4821870E+01 
	};

	if ( n_data < 0 )
	{
		n_data = 0;
	}

	n_data = n_data + 1;

	if ( N_MAX < n_data )
	{
		n_data = 0;
		a = 0.0;
		x = 0.0;
		fx = 0.0;
	}
	else
	{
		a = a_vec[n_data-1];
		x = x_vec[n_data-1];
		fx = fx_vec[n_data-1];
	}
	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    GAMMA_LN1 evaluates ln ( Gamma ( 1 + A ) ), for -0.2 <= A <= 1.25.
//
//  Parameters:
//
//    Input, double %A, defines the argument of the function.
//
//    Output, double GAMMA_LN1, the value of ln ( Gamma ( 1 + A ) ).
//
//****************************************************************************80
double StatClass::gamma_ln1 ( double %a )
{
	double p0 = .577215664901533e+00;
	double p1 = .844203922187225e+00;
	double p2 = -.168860593646662e+00;
	double p3 = -.780427615533591e+00;
	double p4 = -.402055799310489e+00;
	double p5 = -.673562214325671e-01;
	double p6 = -.271935708322958e-02;
	double q1 = .288743195473681e+01;
	double q2 = .312755088914843e+01;
	double q3 = .156875193295039e+01;
	double q4 = .361951990101499e+00;
	double q5 = .325038868253937e-01;
	double q6 = .667465618796164e-03;
	double r0 = .422784335098467e+00;
	double r1 = .848044614534529e+00;
	double r2 = .565221050691933e+00;
	double r3 = .156513060486551e+00;
	double r4 = .170502484022650e-01;
	double r5 = .497958207639485e-03;
	double s1 = .124313399877507e+01;
	double s2 = .548042109832463e+00;
	double s3 = .101552187439830e+00;
	double s4 = .713309612391000e-02;
	double s5 = .116165475989616e-03;
	double gamln1,w,x;

	if (a >= 0.6e0)
	{
		x = a-0.5e0-0.5e0;
		w = (((((r5*x+r4)*x+r3)*x+r2)*x+r1)*x+r0)/(((((s5*x+s4)*x+s3)*x+s2)*x+s1)*x
			+1.0);
		gamln1 = x*w;
	}
	else
	{
		w = ((((((p6*a+p5)*a+p4)*a+p3)*a+p2)*a+p1)*a+p0)/((((((q6*a+q5)*a+
			q4)*a+q3)*a+q2)*a+q1)*a+1.0);
		gamln1 = -(a*w);
	}
	return gamln1;
}

//****************************************************************************80
//
//  Purpose:
// 
//    GAMMA_LOG evaluates ln ( Gamma ( A ) ) for positive A.
//
//  Author:
//
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//
//  Reference:
//
//    Armido DiDinato and Alfred Morris,
//    Algorithm 708: 
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//
//  Parameters:
//
//    Input, double %A, the argument of the function.
//    A should be positive.
//
//    Output, double GAMMA_LOG, the value of ln ( Gamma ( A ) ).
//
//****************************************************************************80
double StatClass::gamma_log ( double %a )
{
	double c0 = .833333333333333e-01;
	double c1 = -.277777777760991e-02;
	double c2 = .793650666825390e-03;
	double c3 = -.595202931351870e-03;
	double c4 = .837308034031215e-03;
	double c5 = -.165322962780713e-02;
	double d = .418938533204673;
	double gamln,t,w;
	int i,n;
	double T1;

	if (a <= 0.8)
	{
		gamln = gamma_ln1 ( a ) - Math::Log ( a );
		return gamln;
	}

	if (a <= 2.25)
	{
		t = a-0.5e0-0.5e0;
		gamln = gamma_ln1 ( t );
		return gamln;
	}

	if (a < 10.0)
	{
		n = ( int ) ( a - 1.25e0 );
		t = a;
		w = 1.0;
		for ( i = 1; i <= n; i++ )
		{
			t -= 1.0;
			w = t*w;
		}
		T1 = t-1.0;
		gamln = gamma_ln1 ( T1 ) + Math::Log ( w );
		return gamln;
	}

	t = Math::Pow(1.0/a, 2.0);
	w = (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0) / a;
	gamln = d + w + (a-0.5e0) * (Math::Log(a)-1.0);
	return gamln;
}

//****************************************************************************80
//
//  Purpose: 
//
//    GAMMA_RAT1 evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
//
//  Parameters:
//
//    Input, double %A, *X, the parameters of the functions.
//    It is assumed that A <= 1.
//
//    Input, double %R, the value Math::Exp(-X) * X**A / Gamma(A).
//
//    Output, double %P, *Q, the values of P(A,X) and Q(A,X).
//
//    Input, double %EPS, the tolerance.
//
//****************************************************************************80
void StatClass::gamma_rat1 ( double %a, double %x, double %r, double %p, double %q,	double %eps )
{
	double a2n,a2nm1,am0,an,an0,b2n,b2nm1,c,cma,g,h,j,l,sum,t,tol,w,z,T1,T3;

	//
	//                SPECIAL CASES
	//
	if  (0.0 == a * x)
	{
		if (x <= a)
		{
			p = 0.0;
			q = 1.0;
		}
		else
		{
			p = 1.0;
			q = 0.0;
		}
		return;
	}

	if (a == 0.5e0)
	{
		if (x >= 0.25e0)
		{
			T3 = Math::Sqrt(x);
			int K3 = 0;
			q = error_fc ( K3, T3 );
			p = 1.0 - q;
		}
		else
		{
			T1 = Math::Sqrt(x);
			p = error_f ( T1 );
			q = 1.0 -p;
		}
		return;
	}

	if (x < 1.1e0)
	{
		//
		//             TAYLOR SERIES FOR P(A,X)/X**A
		//
		an = 3.0e0;
		c = x;
		sum = x/(a+3.0e0);
		tol = 0.1e0*eps/(a+1.0);
		do
		{
			an += 1.0;
			c = -(c*(x/an));
			t = c/(a+an);
			sum += t;
		} while (Math::Abs(t) > tol);

		j = a*x*((sum/6.0e0-0.5e0/(a+2.0e0))*x+1.0/(a+1.0));
		z = a*Math::Log(x);
		h = gam1(a);
		g = 1.0+h;
		
		if (((x < 0.25e0) && (z > -.13394e0)) || (a < x/2.59))
		{
			l = rexp(z);
			w = 0.5e0+(0.5e0+l);
			q = (w*j-l)*g-h;
			if (q < 0.0)
			{
				p = 1.0;
				q = 0.0;
				return;
			}
			p = 0.5e0+(0.5e0-q);
			return;
		}

		w = Math::Exp(z);
		p = w*g*(0.5e0+(0.5e0-j));
		q = 0.5e0+(0.5e0-p);
		return;
	}

	//
	//              CONTINUED FRACTION EXPANSION
	//
	a2nm1 = a2n = 1.0;
	b2nm1 = x;
	b2n = x+(1.0-a);
	c = 1.0;

	do
	{
		a2nm1 = x*a2n+c*a2nm1;
		b2nm1 = x*b2n+c*b2nm1;
		am0 = a2nm1/b2nm1;
		c += 1.0;
		cma = c-a;
		a2n = a2nm1+cma*a2n;
		b2n = b2nm1+cma*b2n;
		an0 = a2n/b2n;
	}
	while (Math::Abs(an0-am0) >= eps*an0);

	q = r*an0;
	p = 0.5e0+(0.5e0-q);
	return;
}

//****************************************************************************80***
//
//  Purpose: 
//
//    GAMMA_VALUES returns some values of the Gamma function.
//
//  Definition:
//
//    GAMMA(Z) = Integral ( 0 <= T < Infinity) T**(Z-1) EXP(-T) dT
//
//  Recursion:
//
//    GAMMA(X+1) = X*GAMMA(X)
//
//  Restrictions:
//
//    0 < X ( a software restriction).
//
//  Special values:
//
//    GAMMA(0.5) = Math::Sqrt(PI)
//
//    For N a positive integer, GAMMA(N+1) = N!, the standard factorial.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//  Parameters:
//
//    Input/output, int %N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double %X, the argument of the function.
//
//    Output, double %FX, the value of the function.
//
//****************************************************************************80***
void StatClass::gamma_values ( int %n_data, double %x, double %fx )
{
	const int N_MAX = 18;

	array<const double> ^fx_vec =
	{ 
		4.590845,     2.218160,     1.489192,     1.164230, 
		1.0000000000, 0.9513507699, 0.9181687424, 0.8974706963, 
		0.8872638175, 0.8862269255, 0.8935153493, 0.9086387329, 
		0.9313837710, 0.9617658319, 1.0000000000, 3.6288000E+05, 
		1.2164510E+17,    8.8417620E+30 
	};
	array<const double> ^x_vec = 
	{ 
		0.2,  0.4,  0.6,  0.8, 
		1.0,  1.1,  1.2,  1.3, 
		1.4,  1.5,  1.6,  1.7, 
		1.8,  1.9,  2.0, 10.0, 
		20.0, 30.0 
	};

	if ( n_data < 0 )
	{
		n_data = 0;
	}

	n_data = n_data + 1;

	if ( N_MAX < n_data )
	{
		n_data = 0;
		x = 0.0;
		fx = 0.0;
	}
	else
	{
		x = x_vec[n_data-1];
		fx = fx_vec[n_data-1];
	}
	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    GAMMA_X evaluates the gamma function.
//
//  Discussion:
//
//    This routine was renamed from "GAMMA" to avoid a conflict with the
//    C/C++ math library routine.
//
//  Author:
//
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//
//  Parameters:
//
//    Input, double %A, the argument of the Gamma function.
//
//    Output, double GAMMA_X, the value of the Gamma function.
//
//****************************************************************************80
double StatClass::gamma_x ( double %a )
{
	double d = .41893853320467274178;
	double pi = 3.1415926535898;
	double r1 = .820756370353826e-03;
	double r2 = -.595156336428591e-03;
	double r3 = .793650663183693e-03;
	double r4 = -.277777777770481e-02;
	double r5 = .833333333333333e-01;

	array<const double> ^p =
	{
		0.539637273585445e-03,
		0.261939260042690e-02,
		0.204493667594920e-01,
		0.730981088720487e-01,
		0.279648642639792,
		0.553413866010467,
		1.0
	};

	array<const double> ^q =
	{
		-0.832979206704073e-03,
		0.470059485860584e-02,
		0.225211131035340e-01,
		-0.170458969313360,
		-0.567902761974940e-01,
		0.113062953091122e+01,
		1.0
	};

	double Xgamm,bot,g,s = 0,t,top,w,x,z;
	int i,j,m,n;

	Xgamm = 0.0;
	x = a;
	if (Math::Abs(a) < 15.0)
	{
		//
		//            EVALUATION OF GAMMA(A) FOR ABS(A) .LT. 15
		//
		t = 1.0;
		m = fifidint(a)-1;

		//
		//     LET T BE THE PRODUCT OF A-J WHEN A .GE. 2
		//

		if (0 == m)
		{
			x -= 1.0;
		}
		else if (m > 0)
		{
			for ( j = 1; j <= m; j++ )
			{
				x -= 1.0;
				t = x*t;
			}

			x -= 1.0;
		}
		else if (m < 0)
		{
			//
			//     LET T BE THE PRODUCT OF A+J WHEN A .LT. 1
			//
			t = a;
			if (a <= 0.0)
			{
				m = -m-1;

				if (m != 0)
				{
					for ( j = 1; j <= m; j++ )
					{
						x += 1.0;
						t = x*t;
					}
				}

				x += 1.0;
				t = x * t;

				if (t == 0.0) return Xgamm;
			}
			//
			//     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
			//     CODE MAY BE OMITTED IF DESIRED.
			//
			if (Math::Abs(t) < 1.e-30)
			{
				if (Math::Abs(t) * double::MaxValue <= 1.0001) return Xgamm;
				Xgamm = 1.0/t;
				return Xgamm;
			}
		}

		//
		//     COMPUTE GAMMA(1 + X) FOR  0 .LE. X .LT. 1
		//
		top = p[0];
		bot = q[0];
		for ( i = 1; i < 7; i++ )
		{
			top = p[i]+x*top;
			bot = q[i]+x*bot;
		}
		Xgamm = top/bot;

		//
		//     TERMINATION
		//
		if (a >= 1.0)
		{
			Xgamm *= t;
			return Xgamm;
		}
		else
		{
			Xgamm /= t;
			return Xgamm;
		}
	}

	//
	//  EVALUATION OF GAMMA(A) FOR ABS(A) .GE. 15
	//
	if (Math::Abs(a) >= 1.0e03)
		return Xgamm;
	
	if (a <= 0.0)
	{
		x = -a;
		n = ( int ) x;
		t = x-(double)n;

		if(t > 0.9)
			t = 1.0-t;

		s = Math::Sin(pi*t)/pi;
		
		if(fifmod(n,2) == 0)
			s = -s;

		if (s == 0.0)
			return Xgamm;
	}

	//
	//     COMPUTE THE MODIFIED ASYMPTOTIC SUM
	//
	t = 1.0/(x*x);
	g = ((((r1*t+r2)*t+r3)*t+r4)*t+r5)/x;

	//
	//  FINAL ASSEMBLY
	//
	z = x;
	g = d + g + (z - 0.5)*(Math::Log(x) - 1.0);
	w = g;
	t = g - w;
	
	int K3 = 0;
	if (w > 0.99999 * exparg(K3))
		return Xgamm;

	Xgamm = Math::Exp(w)*(1.0+t);
	
	if (a < 0.0)
		Xgamm = 1.0/(Xgamm*s)/x;

	return Xgamm;
}

//****************************************************************************80
//
//  Purpose:
// 
//    GSUMLN evaluates the function ln(Gamma(A + B)).
//
//  Discussion:
//
//    GSUMLN is used for 1 <= A <= 2 and 1 <= B <= 2
//
//  Parameters:
//
//    Input, double %A, *B, values whose sum is the argument of
//    the Gamma function.
//
//    Output, double GSUMLN, the value of ln(Gamma(A+B)).
//
//****************************************************************************80
double StatClass::gsumln ( double %a, double %b )
{
	double gsumln, x;

	x = a + b - 2.0;

	if (x <= 0.25)
	{
		gsumln = gamma_ln1 ( 1.0 + x );
		return gsumln;
	}

	if (x <= 1.25)
	{
		gsumln = gamma_ln1 ( x ) + alnrel ( x );
		return gsumln;
	}

	gsumln = gamma_ln1 ( x - 1.0 ) + Math::Log ( x * ( 1.0 + x ) );
	return gsumln;
}

//****************************************************************************80
//
//  Purpose:
//  
//    IPMPAR returns integer machine constants. 
//
//  Discussion:
//
//    Input arguments 1 through 3 are queries about integer arithmetic.
//    We assume integers are represented in the N-digit, base-A form
//
//      sign * ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
//
//    where 0 <= X(0:N-1) < A.
//
//    Then:
//
//      IPMPAR(1) = A, the base of integer arithmetic;
//      IPMPAR(2) = N, the number of base A digits;
//      IPMPAR(3) = A**N - 1, the largest magnitude.
//
//    It is assumed that the single and double precision floating
//    point arithmetics have the same base, say B, and that the
//    nonzero numbers are represented in the form
//
//      sign * (B**E) * (X(1)/B + ... + X(M)/B**M)
//
//    where X(1:M) is one of { 0, 1,..., B-1 }, and 1 <= X(1) and
//    EMIN <= E <= EMAX.
//
//    Input argument 4 is a query about the base of real arithmetic:
//
//      IPMPAR(4) = B, the base of single and double precision arithmetic.
//
//    Input arguments 5 through 7 are queries about single precision
//    floating point arithmetic:
//
//     IPMPAR(5) = M, the number of base B digits for single precision.
//     IPMPAR(6) = EMIN, the smallest exponent E for single precision.
//     IPMPAR(7) = EMAX, the largest exponent E for single precision.
//
//    Input arguments 8 through 10 are queries about double precision
//    floating point arithmetic:
//
//     IPMPAR(8) = M, the number of base B digits for double precision.
//     IPMPAR(9) = EMIN, the smallest exponent E for double precision.
//     IPMPAR(10) = EMAX, the largest exponent E for double precision.
//
//  Reference:
//
//    Phyllis Fox, Andrew Hall, and Norman Schryer,
//    Algorithm 528,
//    Framework for a Portable FORTRAN Subroutine Library,
//    ACM Transactions on Mathematical Software,
//    Volume 4, 1978, pages 176-188.
//
//  Parameters:
//
//    Input, int %I, the index of the desired constant.
//
//    Output, int IPMPAR, the value of the desired constant.
//
//****************************************************************************80
//int StatClass::ipmpar ( const int i )
//{
//	array<int> ^imach = gcnew array<int>(11);
//
//	imach[1] = 2;
//	imach[2] = 31;
//	imach[3] = 2147483647;
//	imach[4] = 2;
//	imach[5] = 24;
//	imach[6] = -125;
//	imach[7] = 128;
//	imach[8] = 53;
//	imach[9] = -1021;
//	imach[10] = 1024;
//
//	return imach[i];
//}


//****************************************************************************80**
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_CDF_VALUES returns values of the negative binomial CDF.
//
//  Discussion:
//
//    Assume that a coin has a probability P of coming up heads on
//    any one trial.  Suppose that we plan to flip the coin until we
//    achieve a total of S heads.  If we let F represent the number of
//    tails that occur in this process, then the value of F satisfies
//    a negative binomial PDF:
//
//      PDF(F,S,P) = Choose ( F from F+S-1 ) * P**S * (1-P)**F
//
//    The negative binomial CDF is the probability that there are F or
//    fewer failures upon the attainment of the S-th success.  Thus,
//
//      CDF(F,S,P) = sum ( 0 <= G <= F ) PDF(G,S,P)
//
//  Modified:
//
//    07 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    F C Powell,
//    Statistical Tables for Sociology, Biology and Physical Sciences,
//    Cambridge University Press, 1982.
//
//  Parameters:
//
//    Input/output, int %N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int %F, the maximum number of failures.
//
//    Output, int %S, the number of successes.
//
//    Output, double %P, the probability of a success on one trial.
//
//    Output, double %CDF, the probability of at most F failures before the 
//    S-th success.
//
//****************************************************************************80**
void StatClass::negative_binomial_cdf_values ( int %n_data, int %f, int %s, double %p, double %cdf )
{
	const int N_MAX = 27;

	array<const double> ^cdf_vec = 
	{
		0.6367, 0.3633, 0.1445, 
		0.5000, 0.2266, 0.0625, 
		0.3438, 0.1094, 0.0156, 
		0.1792, 0.0410, 0.0041, 
		0.0705, 0.0109, 0.0007, 
		0.9862, 0.9150, 0.7472, 
		0.8499, 0.5497, 0.2662, 
		0.6513, 0.2639, 0.0702, 
		1.0000, 0.0199, 0.0001 
	};
	array<int> ^f_vec = 
	{
		4,  3,  2, 
		3,  2,  1, 
		2,  1,  0, 
		2,  1,  0, 
		2,  1,  0, 
		11, 10,  9, 
		17, 16, 15, 
		9,  8,  7, 
		2,  1,  0 
	};
	array<const double> ^p_vec = 
	{
		0.50, 0.50, 0.50, 
		0.50, 0.50, 0.50, 
		0.50, 0.50, 0.50, 
		0.40, 0.40, 0.40, 
		0.30, 0.30, 0.30, 
		0.30, 0.30, 0.30, 
		0.10, 0.10, 0.10, 
		0.10, 0.10, 0.10, 
		0.01, 0.01, 0.01 
	};
	array<const int> ^s_vec = 
	{
		4, 5, 6, 
		4, 5, 6, 
		4, 5, 6, 
		4, 5, 6, 
		4, 5, 6, 
		1, 2, 3, 
		1, 2, 3, 
		1, 2, 3, 
		0, 1, 2 
	};

	if ( n_data < 0 )
	{
		n_data = 0;
	}

	n_data = n_data + 1;

	if ( N_MAX < n_data )
	{
		n_data = 0;
		f = 0;
		s = 0;
		p = 0.0;
		cdf = 0.0;
	}
	else
	{
		f = f_vec[n_data-1];
		s = s_vec[n_data-1];
		p = p_vec[n_data-1];
		cdf = cdf_vec[n_data-1];
	}

	return;
}

//****************************************************************************80***
//
//  Purpose: 
//
//    NORMAL_CDF_VALUES returns some values of the Normal CDF.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//  Parameters:
//
//    Input/output, int %N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double %X, the argument of the function.
//
//    Output double %FX, the value of the function.
//
//****************************************************************************80***
void StatClass::normal_cdf_values ( int %n_data, double %x, double %fx )
{
	const int N_MAX = 13;

	array<const double> ^fx_vec = 
	{ 
		0.500000000000000, 0.539827837277029, 0.579259709439103, 
		0.617911422188953, 0.655421741610324, 0.691462461274013, 
		0.725746882249927, 0.758036347776927, 0.788144601416604, 
		0.815939874653241, 0.841344746068543, 0.933192798731142, 
		0.977249868051821 
	};
	array<const double> ^x_vec = 
	{ 
		0.00, 0.10, 0.20, 
		0.30, 0.40, 0.50, 
		0.60, 0.70, 0.80, 
		0.90, 1.00, 1.50, 
		2.00 
	};

	if ( n_data < 0 )
	{
		n_data = 0;
	}

	n_data = n_data + 1;

	if ( N_MAX < n_data )
	{
		n_data = 0;
		x = 0.0;
		fx = 0.0;
	}
	else
	{
		x = x_vec[n_data-1];
		fx = fx_vec[n_data-1];
	}

	return;
}


//****************************************************************************80***
//
//  Purpose: 
//
//    POISSON_CDF_VALUES returns some values of the Poisson CDF.
//
//  Discussion:
//
//    CDF(X)(A) is the probability of at most X successes in unit time,
//    given that the expected mean number of successes is A.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 653-658.
//
//  Parameters:
//
//    Input/output, int %N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double %A, the parameter of the function.
//
//    Output, int %X, the argument of the function.
//
//    Output, double %FX, the value of the function.
//
//****************************************************************************80***
void StatClass::poisson_cdf_values ( int %n_data, double %a, int %x, double %fx )
{
	const int N_MAX = 21;

	array<const double> ^a_vec = 
	{ 
		0.02, 0.10, 0.10, 0.50, 
		0.50, 0.50, 1.00, 1.00, 
		1.00, 1.00, 2.00, 2.00, 
		2.00, 2.00, 5.00, 5.00, 
		5.00, 5.00, 5.00, 5.00, 
		5.00 
	};
	array<const double> ^fx_vec = 
	{ 
		0.980, 0.905, 0.995, 0.607, 
		0.910, 0.986, 0.368, 0.736, 
		0.920, 0.981, 0.135, 0.406, 
		0.677, 0.857, 0.007, 0.040, 
		0.125, 0.265, 0.441, 0.616, 
		0.762 
	};
	array<int> ^x_vec = 
	{ 
		0, 0, 1, 0, 
		1, 2, 0, 1, 
		2, 3, 0, 1, 
		2, 3, 0, 1, 
		2, 3, 4, 5, 
		6 
	};

	if ( n_data < 0 )
	{
		n_data = 0;
	}

	n_data = n_data + 1;

	if ( N_MAX < n_data )
	{
		n_data = 0;
		a = 0.0;
		x = 0;
		fx = 0.0;
	}
	else
	{
		a = a_vec[n_data-1];
		x = x_vec[n_data-1];
		fx = fx_vec[n_data-1];
	}
	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    PSI evaluates the psi or digamma function, d/dx ln(gamma(x)).
//
//  Discussion:
//
//    The main computation involves evaluation of rational Chebyshev
//    approximations.  PSI was written at Argonne National Laboratory 
//    for FUNPACK, and subsequently modified by A. H. Morris of NSWC.
//
//  Reference:
//
//    William Cody, Strecok and Thacher,
//    Chebyshev Approximations for the Psi Function,
//    Mathematics of Computation,
//    Volume 27, 1973, pages 123-127.
//
//  Parameters:
//
//    Input, double %XX, the argument of the psi function.
//
//    Output, double PSI, the value of the psi function.  PSI 
//    is assigned the value 0 when the psi function is undefined.
//
//****************************************************************************80
double StatClass::psi ( double %xx )
{
	double dx0 =  1.46163214496836;
	double piov4 = .785398163397448;
	array<const double> ^p1 =
	{
		.895385022981970e-02,.477762828042627e+01,.142441585084029e+03,
		.118645200713425e+04,.363351846806499e+04,.413810161269013e+04,
		.130560269827897e+04
	};
	array<const double> ^p2 = 
	{
		-.212940445131011e+01,-.701677227766759e+01,-.448616543918019e+01,
		-.648157123766197
	};
	array<const double> ^q1 = 
	{
		.448452573429826e+02,.520752771467162e+03,.221000799247830e+04,
		.364127349079381e+04,.190831076596300e+04,.691091682714533e-05
	};
	array<const double> ^q2 = 
	{
		.322703493791143e+02,.892920700481861e+02,.546117738103215e+02,
		.777788548522962e+01
	};
	double psi,aug,den,sgn,upper,w,x,xmax1,xmx0,xsmall,z;
	int i,m,n,nq;
	//
	//     MACHINE DEPENDENT CONSTANTS ...
	//        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
	//                 WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED
	//                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
	//                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
	//                 PSI MAY BE REPRESENTED AS ALOG(X).
	//        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X)
	//                 MAY BE REPRESENTED BY 1/X.
	//
	// xmax1 = ipmpar(K1);
	xmax1 = 2147483647;
	xmax1 = Math::Min(xmax1, 1.0 / double::Epsilon);
	xsmall = 1.e-9;
	x = xx;
	aug = 0.0;

	if (x == 0.0) return 0.0;

	if (x < 0.5)
	{
		//
		//     X .LT. 0.5,  USE REFLECTION FORMULA
		//     PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
		//
		if(Math::Abs(x) <= xsmall)
		{
			//
			//     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE
			//     FOR  PI*COTAN(PI*X)
			//
			aug = -(1.0/x);
		}
		else
		{
			//
			//     REDUCTION OF ARGUMENT FOR COTAN
			//
			w = -x;
			sgn = piov4;
			if (w <= 0.0)
			{
				w = -w;
				sgn = -sgn;
			}

			//
			//     MAKE AN ERROR EXIT IF X .LE. -XMAX1
			//
			if (w >= xmax1)
				return 0.0;
			nq = fifidint(w);
			w -= (double)nq;
			nq = fifidint(w*4.0);
			w = 4.0e0*(w-(double)nq*.25);
			//
			//     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
			//     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
			//     QUADRANT AND DETERMINE SIGN
			//
			n = nq/2;
			if(n+n != nq) w = 1.0-w;
			z = piov4*w;
			m = n/2;
			if(m+m != n) sgn = -sgn;
			//
			//     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X)
			//
			n = (nq+1)/2;
			m = n/2;
			m += m;
			if(m == n)
			{
				//
				//     CHECK FOR SINGULARITY
				//
				if (z == 0.0)
					return 0.0;
				//
				//     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
				//     SIN/COS AS A SUBSTITUTE FOR TAN
				//
				aug = sgn*(Math::Cos(z)/Math::Sin(z)*4.0e0);
			}
			else
			{
				aug = sgn*(Math::Sin(z)/Math::Cos(z)*4.0e0);
			}
		}
		x = 1.0-x;
	}

	if (x <= 3.0)
	{
		//
		//     0.5 .LE. X .LE. 3.0
		//
		den = x;
		upper = p1[0]*x;
		for ( i = 1; i <= 5; i++ )
		{
			den = (den+q1[i-1])*x;
			upper = (upper+p1[i+1-1])*x;
		}
		den = (upper+p1[6])/(den+q1[5]);
		xmx0 = x-dx0;
		psi = den*xmx0+aug;
		return psi;
	}

	//
	//     IF X .GE. XMAX1, PSI = LN(X)
	//
	if (x < xmax1)
	{
		//
		//     3.0 .LT. X .LT. XMAX1
		//
		w = 1.0/(x*x);
		den = w;
		upper = p2[0]*w;
		for ( i = 1; i <= 3; i++ )
		{
			den = (den+q2[i-1])*w;
			upper = (upper+p2[i+1-1])*w;
		}
		aug = upper/(den+q2[3])-0.5e0/x+aug;
	}

	psi = aug+Math::Log(x);
	return psi;
}

//****************************************************************************80***
//
//  Purpose: 
//
//    PSI_VALUES returns some values of the Psi or Digamma function.
//
//  Discussion:
//
//    PSI(X) = d LN ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
//
//    PSI(1) = - Euler's constant.
//
//    PSI(X+1) = PSI(X) + 1 / X.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//  Parameters:
//
//    Input/output, int %N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double %X, the argument of the function.
//
//    Output, double %FX, the value of the function.
//
//****************************************************************************80***
void StatClass::psi_values ( int %n_data, double %x, double %fx )
{
	const int N_MAX = 11;

	array<double> ^fx_vec =
	{ 
		-0.5772156649, -0.4237549404, -0.2890398966, 
		-0.1691908889, -0.0613845446, 0.0364899740, 
		0.1260474528,  0.2085478749,  0.2849914333, 
		0.3561841612,  0.4227843351
	};

	array<double> ^x_vec =
	{ 
		1.0,  1.1,  1.2,  
		1.3,  1.4,  1.5,  
		1.6,  1.7,  1.8,  
		1.9,  2.0
	};

	if ( n_data < 0 )
	{
		n_data = 0;
	}

	n_data = n_data + 1;

	if ( N_MAX < n_data )
	{
		n_data = 0;
		x = 0.0;
		fx = 0.0;
	}
	else
	{
		x = x_vec[n_data-1];
		fx = fx_vec[n_data-1];
	}
	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    RCOMP evaluates Math::Exp(-X) * X**A / Gamma(A).
//
//  Parameters:
//
//    Input, double %A, *X, arguments of the quantity to be computed.
//
//    Output, double RCOMP, the value of Math::Exp(-X) * X**A / Gamma(A).
//
//  Local parameters:
//
//    RT2PIN = 1/SQRT(2*PI)
//
//****************************************************************************80
double StatClass::rcomp ( double %a, double %x )
{
	double rt2pin = .398942280401433;
	double rcomp,t,t1,u;
	rcomp = 0.0;
	if (a < 20.0)
	{
		t = a*Math::Log(x)-x;
		if (a < 1.0)
		{
			rcomp = a * Math::Exp(t) * (1.0 + gam1(a));
			return rcomp;
		}
		else
		{
			rcomp = Math::Exp(t )/ gamma_x(a);
			return rcomp;
		}
	}
	else
	{
		u = x / a;
		if (u == 0.0) return rcomp;
		t = Math::Pow(1.0 / a , 2.0);
		t1 = (((0.75 * t - 1.0) * t + 3.5) * t - 105.0) / (a * 1260.0);
		t1 -= (a * rlog(u));
		rcomp = rt2pin * Math::Sqrt(a) * Math::Exp(t1);
		return rcomp;
	}
}

//****************************************************************************80
//
//  Purpose:
// 
//    REXP evaluates the function EXP(X) - 1.
//
//  Modified:
//
//    09 December 1999
//
//  Parameters:
//
//    Input, double %X, the argument of the function.
//
//    Output, double REXP, the value of EXP(X)-1.
//
//****************************************************************************80
double StatClass::rexp ( double %x )
{
	double p1 = .914041914819518e-09;
	double p2 = .238082361044469e-01;
	double q1 = -.499999999085958e+00;
	double q2 = .107141568980644e+00;
	double q3 = -.119041179760821e-01;
	double q4 = .595130811860248e-03;
	double rexp,w;

	if (Math::Abs(x) > 0.15)
	{
		w = Math::Exp(x);
		if(x > 0.0)
			rexp = w*(1.0 - 1.0/w);
		else
			rexp = w - 1.0;
		return rexp;
	}

	rexp = x*(((p2*x+p1)*x+1.0)/((((q4*x+q3)*x+q2)*x+q1)*x+1.0));
	return rexp;

}

//****************************************************************************80
//
//  Purpose:
// 
//    RLOG computes  X - 1 - LN(X).
//
//  Modified:
//
//    09 December 1999
//
//  Parameters:
//
//    Input, double %X, the argument of the function.
//
//    Output, double RLOG, the value of the function.
//
//****************************************************************************80
double StatClass::rlog ( double %x )
{
	double a = .566749439387324e-01;
	double b = .456512608815524e-01;
	double p0 = .333333333333333e+00;
	double p1 = -.224696413112536e+00;
	double p2 = .620886815375787e-02;
	double q1 = -.127408923933623e+01;
	double q2 = .354508718369557e+00;
	double rlog,r,t,u,w,w1;

	if ((x > 0.61) && (x < 1.57))
	{
		//
		//              ARGUMENT REDUCTION
		//
		if (x < 0.82e0)
		{
			u = x-0.7e0;
			u /= 0.7e0;
			w1 = a-u*0.3e0;
		}
		else if(x > 1.18e0)
		{
			u = 0.75e0 * x - 1.0;
			w1 = b+u/3.0;
		}
		else
		{
			u = x - 1.0;
			w1 = 0.0;
		}

		//
		//               SERIES EXPANSION
		//
		r = u/(u+2.0e0);
		t = r*r;
		w = ((p2*t+p1)*t+p0)/((q2*t+q1)*t+1.0);
		rlog = 2.0e0*t*(1.0/(1.0-r)-r*w)+w1;
	}
	else
	{
		if (x <= 0.0)
			throw gcnew ArgumentException(L"x must be greater than 0.0", L"x");

		r = x - 1.0;
		rlog = r-Math::Log(x);
	}
	return rlog;
}

//****************************************************************************80
//
//  Purpose:
// 
//    RLOG1 evaluates the function X - ln ( 1 + X ).
//
//  Parameters:
//
//    Input, double %X, the argument.
//
//    Output, double RLOG1, the value of X - ln ( 1 + X ).
//
//****************************************************************************80
double StatClass::rlog1 ( double %x )
{
	double a = .566749439387324e-01;
	double b = .456512608815524e-01;
	double p0 = .333333333333333;
	double p1 = -.224696413112536;
	double p2 = .620886815375787e-02;
	double q1 = -.127408923933623e+01;
	double q2 = .354508718369557;
	double rlog1,h,r,t,w,w1;

	if (x >= -0.39 && x <= 0.57)
	{
		//
		//              ARGUMENT REDUCTION
		//
		if(x < -0.18)
		{
			h = x + 0.3;
			h /= 0.7;
			w1 = a - h * 0.3;
		}
		else if (x > 0.18e0)
		{
			h = 0.75 * x - 0.25;
			w1 = b + h / 3.0;
		}
		else
		{
			h = x;
			w1 = 0.0;
		}
		//
		//               SERIES EXPANSION
		//
		r = h / (h + 2.0);
		t = r * r;
		w = ((p2 * t + p1) * t + p0)/((q2 * t + q1) * t + 1.0);
		rlog1 = 2.0 * t * (1.0 / (1.0 - r) - r * w) + w1;
		return rlog1;
	}
	else
	{
		w = x + 1.0;
		if (w <= 0.0)
			throw gcnew ArgumentException("x must be > -1", "x");
		rlog1 = x - Math::Log(w);
		return rlog1;
	}
}

//****************************************************************************80***
//
//  Purpose: 
//
//    STUDENT_CDF_VALUES returns some values of the Student CDF.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//  Parameters:
//
//    Input/output, int %N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int %A, the parameter of the function.
//
//    Output, double %X, the argument of the function.
//
//    Output, double %FX, the value of the function.
//
//****************************************************************************80***
void StatClass::student_cdf_values ( int %n_data, int %a, double %x, double %fx )
{
	const int N_MAX = 13;

	array<int> ^a_vec =	{	1,	2,	3,	4,	5,	2,	5,	2,	5,	2,	3,	4, 	5 };

	array<double> ^fx_vec =	{ 0.60,	0.60, 0.60,	0.60, 0.60,	0.75, 0.75,	0.95, 0.95,	0.99, 0.99,	0.99, 0.99 };

	array<double> ^x_vec =	{ 0.325, 0.289,	0.277, 0.271, 0.267, 0.816, 0.727, 2.920, 2.015, 6.965,	4.541, 3.747, 3.365	};

	if ( n_data < 0 )
		n_data = 1;
	else
		++n_data;

	if ( N_MAX < n_data )
	{
		n_data = 0;
		a = 0;
		x = 0.0;
		fx = 0.0;
	}
	else
	{
		a = a_vec[n_data-1];
		x = x_vec[n_data-1];
		fx = fx_vec[n_data-1];
	}

	return;
}

//****************************************************************************80
//
//  Purpose:
// 
//    STVALN provides starting values for the inverse of the normal distribution.
//
//  Discussion:
//
//    The routine returns X such that 
//      P = CUMNOR(X),  
//    that is, 
//      P = Integral from -infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU.
//
//  Reference:
//
//    Kennedy and Gentle,
//    Statistical Computing, 
//    Marcel Dekker, NY, 1980, page 95,
//    QA276.4  K46
//
//  Parameters:
//
//    Input, double %P, the probability whose normal deviate 
//    is sought.
//
//    Output, double STVALN, the normal deviate whose probability
//    is P.
//
//****************************************************************************80
double StatClass::stvaln ( double %p )
{
	array<double> ^xden = 
	{
		0.993484626060e-1,
		0.588581570495e0,
		0.531103462366e0,
		0.103537752850e0,
		0.38560700634e-2
	};
	array<double> ^xnum = 
	{
		-0.322232431088e0,
		-1.000000000000e0,
		-0.342242088547e0,
		-0.204231210245e-1,
		-0.453642210148e-4
	};

	int K1 = 5;
	double stvaln = 0;
	double sign = 0;
	double y = 0;
	double z = 0;

	if ((p < 0.0 ) || (p > 1.0))
		throw gcnew ArgumentException(L"stvaln: argument must be between 0 and 1.", L"p");

	if(p <= 0.5e0)
	{
		sign = -1.0;
		z = p;
	}
	else
	{
		sign = 1.0;
		z = 1.0-p;
	}

	y = Math::Sqrt(-(2.0e0*Math::Log(z)));
	stvaln = y + eval_pol( xnum, K1, y ) / eval_pol( xden, K1, y );
	stvaln = sign * stvaln;
	return stvaln;
}
//**************************************************************************80
