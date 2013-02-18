// dcdflib.h

#pragma once

using namespace System;

namespace dcdflib {

	public ref class StatClass
	{
	private:

		/* 
		Methods E0000 and E0001 are used to find the zero of a
		function f(x). The methods use reverse communication, which
		means they return to the caller everytime they require a new
		f(x), x pair. The functions utilized static variables to
		maintain state between calls.
		By moving the static declarations to private instance variables,
		the methods should be thread safe when utilized in a .Net environment.
		*/

		// These were static declarations in StatClass::E0001
		// I moved them to instance variables to improve thread
		// safety (hopefully)
		double E0001_a;
		double E0001_abstol;
		double E0001_b;
		double E0001_c;
		double E0001_d;
		double E0001_fa;
		double E0001_fb;
		double E0001_fc;
		double E0001_fd;
		double E0001_fda;
		double E0001_fdb;
		double E0001_m;
		double E0001_mb;
		double E0001_p;
		double E0001_q;
		double E0001_reltol;
		double E0001_tol;
		double E0001_w;
		double E0001_xxhi;
		double E0001_xxlo;
		int E0001_ext;
		int E0001_lastState;
		unsigned long E0001_first;
		unsigned long E0001_qrzero;

		// These were static declarations in StatClass::E0000
		double E0000_absstp;
		double E0000_abstol;
		double E0000_big;
		double E0000_fbig;
		double E0000_fsmall;
		double E0000_relstp;
		double E0000_reltol;
		double E0000_small;
		double E0000_step;
		double E0000_stpmul;
		double E0000_xhi;
		double E0000_xlb;
		double E0000_xlo;
		double E0000_xsave;
		double E0000_xub;
		double E0000_yy;
		int E0000_lastState;
		unsigned long E0000_qbdd;
		unsigned long E0000_qcond;
		unsigned long E0000_qdum1;
		unsigned long E0000_qdum2;
		unsigned long E0000_qincr;
		unsigned long E0000_qlim;
		unsigned long E0000_qok;
		unsigned long E0000_qup;
		double tol;
		double atol;
		double zero;
		double inf;
		double one;
		double tent4;
		double oneless;

	public:
		///<summary>
		/// Constructor
		///</summary>
		StatClass();
		
		///<summary>
		/// algdiv computes ln ( Gamma ( B ) / Gamma ( A + B ) ) when 8 .le. B
		///</summary>
		///<remarks>
		///    In this algorithm, DEL(X) is the function defined by
		///    ln ( Gamma(X) ) = ( X - 0.5 ) * ln ( X ) - X + 0.5 * ln ( 2 * PI ) + DEL(X)
		///</remarks>
		///<param name="a">denominator term</param>
		///<param name="b">numerator term</param>
		///<returns>The value of ln(Gamma(a)/Gamma(a+b))</returns>
		double algdiv ( double %a, double %b );

		///<summary>
		/// alnrel evaluates the function ln ( 1 + A )
		///</summary>
		///<remarks>
		///  Reference:
		///
		///    Armido DiDinato, Alfred Morris,
		///    Algorithm 708: 
		///    Significant Digit Computation of the Incomplete Beta Function Ratios,
		///    ACM Transactions on Mathematical Software,
		///    Volume 18, 1993, pages 360-373.
		///</remarks>
		///<param name="a">Must be .gt. -1.0</param>
		///<returns>The value of ln ( 1 + A )</returns>
		double alnrel ( double %a );

		///<summary>
		///APSER computes the incomplete beta ratio I(SUB(1-X))(B,A)
		///</summary>
		///<remarks>
		///    APSER is used only for cases where
		///      A &lt;= min ( EPS, EPS * B ), 
		///      B * X &lt;= 1, and 
		///      X &lt;= 0.5
		///</remarks>
		///<param name="a">Parmeter of incomplete beta ratio</param>		
		///<param name="b">Parmeter of incomplete beta ratio</param>		
		///<param name="x">Parmeter of incomplete beta ratio</param>		
		///<param name="eps">Calculation tolerance</param>		
		///<returns>The computed value of the incomplete beta ratio</returns>
		double apser ( double %a, double %b, double %x, double %eps );

		///<summary>BCORR evaluates DEL(A0) + DEL(B0) - DEL(A0 + B0)</summary>
		///<remarks>
		///    The function DEL(A) is a remainder term that is used in the expression:
		///
		///      ln ( Gamma ( A ) ) = ( A - 0.5 ) * ln ( A ) 
		///        - A + 0.5 * ln ( 2 * PI ) + DEL ( A ),
		///
		///    or, in other words, DEL ( A ) is defined as:
		///
		///      DEL ( A ) = ln ( Gamma ( A ) ) - ( A - 0.5 ) * ln ( A ) 
		///        + A + 0.5 * ln ( 2 * PI ).
		///</remarks>
		///<param name="a0">Argument</param>		
		///<param name="b0">Argument</param>		
		///<returns>The value of the function</returns>
		double bcorr ( double %a0, double %b0 );

		///<summary>BETA evaluates the beta function</summary>
		///<param name="a">Argument or the beta function</param>		
		///<param name="b">Argument or the beta function</param>		
		///<returns>The value of the beta function</returns>
		double beta ( double a, double b );

		///<summary>BETA_ASYM computes an asymptotic expansion for IX(A,B), for large A and B</summary>
		///<remarks>
		///The input parameters A and B should be nonnegative
		///It is assumed that both A and B are greater than or equal to 15
		///It is assumed that 0 &lt;= lambda
		///</remarks>
		///<param name="a">&gt;0, preferrably &gt;15</param>
		///<param name="b">&gt;0, preferrably &gt;15</param>
		///<param name="eps">The calculation tolerance</param>
		///<param name="lambda">The value of ( A + B ) * Y - B</param>
		///<returns>Asymptotic expansion of IX(A,B)</returns>
		double beta_asym ( double %a, double %b, double %lambda, double %eps );

		///<summary>BETA_FRAC evaluates a continued fraction expansion for IX(A,B)</summary>
		///<param name="a">&gt;0, assumed also to be &gt;15</param>
		///<param name="b">&gt;0, assumed also to be &gt;15</param>
		///<param name="x">0 &lt;= x &lt;= 1</param>
		///<param name="y">y = 1 - x</param>
		///<param name="lambda">The value of ( A + B ) * Y - B</param>
		///<param name="eps">The calculation tolerance</param>
		double beta_frac ( double %a, double %b, double %x, double %y, double %lambda, double %eps );
		
		///<summary>BETA_GRAT evaluates an asymptotic expansion for IX(A,B)</summary>
		///<param name="a">15 &lt;= a</param>
		///<param name="b">b &lt;= 1</param>
		///<param name="x">0 &lt;= x &lt;= 1</param>
		///<param name="y">y = 1 - x</param>
		///<param name="w">a quantity to which the result of the computation is to be added on output</param>
		///<param name="eps">The calculation tolerance</param>
		///<param name="ierr">0 no error, 1 error detected</param>
		///<remarks>It is assumed that 15 &lt;= A and B &lt;= 1, and that B is less than A</remarks>
		void beta_grat ( double %a, double %b, double %x, double %y, double %w,	double %eps, int %ierr );
		
		///<summary>BETA_INC evaluates the incomplete beta function IX(A,B)</summary>
		///<remarks>    Alfred H Morris, Jr,
		///    Naval Surface Weapons Center,
		///    Dahlgren, Virginia.
		///</remarks>
		///<param name="a">0 &lt; a</param>
		///<param name="b">0 &lt; b</param>
		///<param name="x">0 &lt;= x &lt;= 1</param>
		///<param name="y">y = 1 - x</param>
		///<param name="w">the value of IX(A,B)</param>
		///<param name="w1">the value of 1 - IX(A,B)</param>
		///<param name="ierr">error flag
		///<list type="bullet">
		/// <item>
		/// <description>0, no error was detected</description>
		/// </item>
		/// <item>
		/// <description>1, A or B is negative</description>
		/// </item>
		/// <item>
		/// <description>2, A = B = 0</description>
		/// </item>
		/// <item>
		/// <description>3, X &lt; 0 or 1 &lt; X</description>
		/// </item>
		/// <item>
		/// <description>4, Y &lt; 0 or 1 &lt; Y</description>
		/// </item>
		/// <item>
		/// <description>5, X + Y /= 1</description>
		/// </item>
		/// <item>
		/// <description>6, X = A = 0</description>
		/// </item>
		/// <item>
		/// <description>7, Y = B = 0</description>
		/// </item>
		///</list>
		///</param>
		void beta_inc ( double %a, double %b, double %x, double %y, double %w, double %w1, int %ierr );
		
		///<summary>BETA_INC_VALUES returns some test values of the incomplete Beta function</summary>
		///<remarks>
		//    The incomplete Beta function may be written
		///
		///      BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
		///                      / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
		///
		///    Thus,
		///
		///      BETA_INC(A,B,0.0) = 0.0
		///      BETA_INC(A,B,1.0) = 1.0
		///
		///    Note that in Mathematica, the expressions:
		///
		///      BETA[A,B]   = Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
		///      BETA[X,A,B] = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
		///
		///    and thus, to evaluate the incomplete Beta function requires:
		///
		///      BETA_INC(A,B,X) = BETA[X,A,B] / BETA[A,B]
		///
		///  Reference:
		///
		///    Milton Abramowitz and Irene Stegun,
		///    Handbook of Mathematical Functions,
		///    US Department of Commerce, 1964.
		///
		///    Karl Pearson,
		///    Tables of the Incomplete Beta Function,
		///    Cambridge University Press, 1968.
		///</remarks>
		///<param name="n_data">The user sets N_DATA to 0 before the
		///    first call.  On each call, the routine increments N_DATA by 1, and
		///    returns the corresponding data; when there is no more data, the
		///    output value of N_DATA will be 0 again.
		///</param>
		///<param name="a">a value for index n_data</param>
		///<param name="b">b value for index n_data</param>
		///<param name="x">x function argument</param>
		///<param name="fx">fx function argument</param>
		void beta_inc_values ( int %n_data, double %a, double %b, double %x, double %fx );
		
		///<summary>BETA_LOG evaluates the logarithm of the beta function</summary>
		///<remarks>
		///    Armido DiDinato and Alfred Morris,
		///    Algorithm 708: 
		///    Significant Digit Computation of the Incomplete Beta Function Ratios,
		///    ACM Transactions on Mathematical Software, 
		///    Volume 18, 1993, pages 360-373.
		///</remarks>
		///<param name="a0">0 &lt; a0</param>
		///<param name="b0">0 &lt; b0</param>
		///<returns>logarithm of the beta function</returns>
		double beta_log ( double %a0, double %b0 );
		
		///<summary>BETA_PSER uses a power series expansion to evaluate IX(A,B)(X)</summary>
		///<remarks>BETA_PSER is used when B &lt;= 1 or B*X &lt;= 0.7</remarks>
		///<param name="a">a value for integral limit</param>
		///<param name="b">b value for integral limit</param>
		///<param name="x">x function argument</param>
		///<param name="eps">calculation tolerance</param>
		///<returns>function value</returns>
		double beta_pser ( double %a, double %b, double %x, double %eps );

		///<summary>BETA_RCOMP evaluates X**A * Y**B / Beta(A,B)</summary>
		///<param name="a">a value for integral limit</param>
		///<param name="b">b value for integral limit</param>
		///<param name="x">component of numerator</param>
		///<param name="y">component of numerator</param>
		///<returns>the value of evaluates X**A * Y**B / Beta(A,B)</returns>
		double beta_rcomp ( double %a, double %b, double %x, double %y );

		///<summary>BETA_RCOMP1 evaluates Math::Exp(MU) * X**A * Y**B / Beta(A,B)</summary>
		///<param name="mu">mu integer</param>
		///<param name="a">a value for integral limit</param>
		///<param name="b">b value for integral limit</param>
		///<param name="x">component of numerator</param>
		///<param name="y">component of numerator</param>
		///<returns>the value of Math::Exp(MU) * X**A * Y**B / Beta(A,B)</returns>
		double beta_rcomp1 ( int %mu, double %a, double %b, double %x, double %y );

		///<summary>BETA_UP evaluates IX(A,B) - IX(A+N,B) where N is a positive integer</summary>
		///<param name="a">a value for integral limit</param>
		///<param name="b">b value for integral limit</param>
		///<param name="x">component of numerator</param>
		///<param name="y">component of numerator</param>
		///<param name="n">integer</param>
		///<param name="eps">tolerance</param>
		///<returns>the value of Math::Exp(MU) * X**A * Y**B / Beta(A,B)</returns>
		double beta_up ( double %a, double %b, double %x, double %y, int %n, double %eps );

		///<summary>BINOMIAL_CDF_VALUES returns some values of the binomial CDF</summary>
		///<remarks>
		///CDF(X)(A,B) is the probability of at most X successes in A trials,
		///given that the probability of success on a single trial is B.
		/// Reference:
		///
		///   Milton Abramowitz and Irene Stegun,
		///   Handbook of Mathematical Functions,
		///   US Department of Commerce, 1964.
		///
		///   Daniel Zwillinger,
		///   CRC Standard Mathematical Tables and Formulae,
		///   30th Edition, CRC Press, 1996, pages 651-652.
		///</remarks>
		///<param name="n_data">The user sets N_DATA to 0 before the
		///    first call.  On each call, the routine increments N_DATA by 1, and
		///    returns the corresponding data; when there is no more data, the
		///    output value of N_DATA will be 0 again</param>
		///<param name="a">a value for integral limit</param>
		///<param name="b">b value for integral limit</param>
		///<param name="x">component of numerator</param>
		///<param name="fx">function value</param>
		///<returns>the value of Math::Exp(MU) * X**A * Y**B / Beta(A,B)</returns>
		void binomial_cdf_values ( int %n_data, int %a, double %b, int %x, double %fx );

		///<summary>
		///CDFBET evaluates the CDF of the Beta Distribution
		///</summary>
		///<remarks>
		///    This routine calculates any one parameter of the beta distribution 
		///    given the others.
		///
		///    The value P of the cumulative distribution function is calculated 
		///    directly by code associated with the reference.
		///
		///    Computation of the other parameters involves a seach for a value that
		///    produces the desired value of P.  The search relies on the
		///    monotonicity of P with respect to the other parameters.
		///
		///    The beta density is proportional to t%(A-1) * (1-t)%(B-1).
		///
		///  Reference:
		///
		///    Armido DiDinato and Alfred Morris,
		///    Algorithm 708: 
		///    Significant Digit Computation of the Incomplete Beta Function Ratios,
		///    ACM Transactions on Mathematical Software, 
		///    Volume 18, 1993, pages 360-373.
		///</remarks>
		///<param name="which">which indicates which of the next four argument
		///    values is to be calculated from the others.
		///<list type="table">
		///<listheader>
		///<term>value</term>
		///<description>Determines the calculation path</description>
		///</listheader>
		///<item>
		///<term>1</term>
		///<description>Calculate P and Q from X, Y, A and B</description>
		///</item>
		///<item>
		///<term>2</term>
		///<description>Calculate X and Y from P, Q, A and B</description>
		///</item>
		///<item>
		///<term>3</term>
		///<description>3, Calculate A from P, Q, X, Y and B</description>
		///</item>
		///<item>
		///<term>4</term>
		///<description>Calculate B from P, Q, X, Y and A</description>
		///</item>
		///</list>
		///</param>
		///<param name="p">The integral from 0 to X of the chi-square distribution: 0.0 to 1.0</param>
		///<param name="q">1 - p: 0.0 to 1.0</param>
		///<param name="x">the upper limit of integration of the beta density: 0.0 to 1.0</param>
		///<param name="y">1 - x: 0.0 to 1.0</param>
		///<param name="a">the first parameter of the beta density: 0.0 to + infinity</param>
		///<param name="b">the second parameter of the beta density: 0.0 to + infinity</param>
		///<param name="status">reports the status of the computation.
		///     0, if the calculation completed correctly;
		///    -I, if the input parameter number I is out of range;
		///    +1, if the answer appears to be lower than lowest search bound;
		///    +2, if the answer appears to be higher than greatest search bound;
		///    +3, if P + Q /= 1;
		///    +4, if X + Y /= 1.
		///</param>
		///<param name="bound">is only defined if STATUS is nonzero.
		///    If STATUS is negative, then this is the value exceeded by parameter I.
		///    if STATUS is 1 or 2, this is the search bound that was exceeded.
		///</param>
		void cdfbet ( int %which, double %p, double %q, double %x, double %y, double %a, double %b, int %status, double %bound );
		
		///<summary>CDFBIN evaluates the CDF of the Binomial distribution</summary>
		///<remarks>This routine calculates any one parameter of the binomial distribution 
		///    given the others.
		///
		///    The value P of the cumulative distribution function is calculated 
		///    directly.
		///
		///    Computation of the other parameters involves a seach for a value that
		///    produces the desired value of P.  The search relies on the
		///    monotonicity of P with respect to the other parameters.
		///
		///    P is the probablility of S or fewer successes in XN binomial trials,
		///    each trial having an individual probability of success of PR.  
		///
		///    Milton Abramowitz and Irene Stegun,
		///    Handbook of Mathematical Functions 
		///    1966, Formula 26.5.24.
		///</remarks>
		///<param name="which">indicates which of argument values is to 
		///    be calculated from the others.
		///    1: Calculate P and Q from S, XN, PR and OMPR;
		///    2: Calculate S from P, Q, XN, PR and OMPR;
		///    3: Calculate XN from P, Q, S, PR and OMPR;
		///    4: Calculate PR and OMPR from P, Q, S and XN.
		///</param>
		///<param name="p">the cumulation, from 0 to S,
		///    of the binomial distribution.  If P is an input value, it should 
		///    lie in the range [0,1].
		///</param>
		///<param name="q">equal to 1-P</param>
		///<param name="s">the number of successes observed.
		///    Whether this is an input or output value, it should lie in the
		///    range [0,XN].
		///</param>
		///<param name="xn">the number of binomial trials.
		///    If this is an input value it should lie in the range: (0, +infinity).
		///    If it is an output value it will be searched for in the
		///    range [1.0D-300, 1.0D+300].
		///</param>
		///<param name="pr">the probability of success in each 
		///    binomial trial.  Whether this is an input or output value, it should
		///    lie in the range: [0,1].
		///</param>
		///<param name="ompr">equal to 1-PR.  Whether this is an
		///    input or output value, it should lie in the range [0,1].  Also, it should
		///    be the case that PR + OMPR = 1.
		///</param>
		///<param name="status">reports the status of the computation.
		///     0, if the calculation completed correctly;
		///    -I, if the input parameter number I is out of range;
		///    +1, if the answer appears to be lower than lowest search bound;
		///    +2, if the answer appears to be higher than greatest search bound;
		///    +3, if P + Q /= 1;
		///    +4, if PR + OMPR /= 1.
		///</param>
		///<param name="bound">is only defined if STATUS is nonzero.
		///    If STATUS is negative, then this is the value exceeded by parameter I.
		///    if STATUS is 1 or 2, this is the search bound that was exceeded.
		///</param>
		void cdfbin ( int %which, double %p, double %q, double %s, double %xn, double %pr, double %ompr, int %status, double %bound );
		
		///<summary> 
		///CDFCHI evaluates the CDF of the chi square distribution.
		///</summary>
		///<remarks>
		///    This routine calculates any one parameter of the chi square distribution 
		///    given the others.
		///
		///    The value P of the cumulative distribution function is calculated 
		///    directly.
		///
		///    Computation of the other parameters involves a seach for a value that
		///    produces the desired value of P.  The search relies on the
		///    monotonicity of P with respect to the other parameters.
		///
		///    The CDF of the chi square distribution can be evaluated 
		///    within Mathematica by commands such as:
		///
		///      Needs["Statistics`ContinuousDistributions`"]
		///      CDF [ ChiSquareDistribution [ DF ], X ]
		///
		///    Milton Abramowitz and Irene Stegun,
		///    Handbook of Mathematical Functions 
		///    1966, Formula 26.4.19.
		///
		///    Stephen Wolfram,
		///    The Mathematica Book,
		///    Fourth Edition,
		///    Wolfram Media / Cambridge University Press, 1999.
		///</remarks>
		///<param name="which">indicates which argument is to be calculated
		///    from the others.
		///    1: Calculate P and Q from X and DF;
		///    2: Calculate X from P, Q and DF;
		///    3: Calculate DF from P, Q and X.
		///</param>
		///<param name="p">the integral from 0 to X of 
		///    the chi-square distribution.  If this is an input value, it should
		///    lie in the range [0,1].
		///</param>
		///<param name="q">equal to 1 - p</param>
		///<param name="x">the upper limit of integration 
		///    of the chi-square distribution.  If this is an input 
		///    value, it should lie in the range: [0, +infinity).  If it is an output
		///    value, it will be searched for in the range: [0,1.0D+300].
		///</param>
		///<param name="df">the degrees of freedom of the
		///    chi-square distribution.  If this is an input value, it should lie
		///    in the range: (0, +infinity).  If it is an output value, it will be
		///    searched for in the range: [ 1.0D-300, 1.0D+300].
		///</param>
		///<param name="status">reports the status of the computation.
		///     0, if the calculation completed correctly;
		///    -I, if the input parameter number I is out of range;
		///    +1, if the answer appears to be lower than lowest search bound;
		///   +2, if the answer appears to be higher than greatest search bound;
		///    +3, if P + Q /= 1;
		///    +10, an error was returned from CUMGAM.
		///</param>
		///<param name="bound">is only defined if STATUS is nonzero.
		///    If STATUS is negative, then this is the value exceeded by parameter I.
		///    if STATUS is 1 or 2, this is the search bound that was exceeded.
		///</param>
		void cdfchi ( int %which, double %p, double %q, double %x, double %df, int %status, double %bound );
		
		///<summary>CDFCHN evaluates the CDF of the Noncentral Chi-Square</summary>
		///<remarks>
		///    This routine calculates any one parameter of the noncentral chi-square
		///    distribution given values for the others.
		///
		///    The value P of the cumulative distribution function is calculated 
		///    directly.
		///
		///    Computation of the other parameters involves a seach for a value that
		///    produces the desired value of P.  The search relies on the
		///    monotonicity of P with respect to the other parameters.
		///
		///    The computation time required for this routine is proportional
		///    to the noncentrality parameter (PNONC).  Very large values of
		///    this parameter can consume immense computer resources.  This is
		///    why the search range is bounded by 10,000.
		///
		///    The CDF of the noncentral chi square distribution can be evaluated 
		///    within Mathematica by commands such as:
		///
		///      Needs["Statistics`ContinuousDistributions`"]
		///      CDF[ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
		///
		///    Milton Abramowitz and Irene Stegun,
		///    Handbook of Mathematical Functions 
		///    1966, Formula 26.5.25.
		///
		///    Stephen Wolfram,
		///    The Mathematica Book,
		///    Fourth Edition,
		///    Wolfram Media / Cambridge University Press, 1999.
		///</remarks>
		///<param name="which">indicates which argument is to be calculated
		///    from the others.
		///    1: Calculate P and Q from X, DF and PNONC;
		///    2: Calculate X from P, DF and PNONC;
		///    3: Calculate DF from P, X and PNONC;
		///    4: Calculate PNONC from P, X and DF.
		///</param>
		///<param name="p">the integral from 0 to X of 
		///    the noncentral chi-square distribution.  If this is an input
		///    value, it should lie in the range: [0, 1.0-1.0D-16).
		///</param>
		///<param name="q">is generally not used by this 
		///    subroutine and is only included for similarity with other routines.
		///    However, if P is to be computed, then a value will also be computed
		///    for Q.
		///</param>
		///<param name="x">the upper limit of integration of the 
		///    noncentral chi-square distribution.  If this is an input value, it
		///    should lie in the range: [0, +infinity).  If it is an output value,
		///    it will be sought in the range: [0,1.0D+300].
		///</param>
		///<param name="df">the number of degrees of freedom 
		///    of the noncentral chi-square distribution.  If this is an input value,
		///    it should lie in the range: (0, +infinity).  If it is an output value,
		///   it will be searched for in the range: [ 1.0D-300, 1.0D+300].
		///</param>
		///<param name="pnonc">the noncentrality parameter of 
		///    the noncentral chi-square distribution.  If this is an input value, it
		///    should lie in the range: [0, +infinity).  If it is an output value,
		///    it will be searched for in the range: [0,1.0D+4]
		///</param>
		///<param name="status">reports on the calculation.
		///    0, if calculation completed correctly;
		///    -I, if input parameter number I is out of range;
		///    1, if the answer appears to be lower than the lowest search bound;
		///    2, if the answer appears to be higher than the greatest search bound.
		///</param>
		///<param name="bound">is only defined if STATUS is nonzero.
		///    If STATUS is negative, then this is the value exceeded by parameter I.
		///    if STATUS is 1 or 2, this is the search bound that was exceeded.
		///</param>		
		void cdfchn ( int %which, double %p, double %q, double %x, double %df, double %pnonc, int %status, double %bound );
		
		///<summary>CDFF evaluates the CDF of the F distribution</summary>
		///<remarks>
		///    This routine calculates any one parameter of the F distribution 
		///    given the others.
		///
		///    The value P of the cumulative distribution function is calculated 
		///    directly.
		///
		///    Computation of the other parameters involves a seach for a value that
		///    produces the desired value of P.  The search relies on the
		///    monotonicity of P with respect to the other parameters.
		///
		///    The value of the cumulative F distribution is not necessarily
		///    monotone in either degree of freedom.  There thus may be two
		///    values that provide a given CDF value.  This routine assumes
		///    monotonicity and will find an arbitrary one of the two values.
		///
		///  Modified:
		///
		///    14 April 2007
		///
		///  Reference:
		///
		///    Milton Abramowitz, Irene Stegun,
		///    Handbook of Mathematical Functions 
		///    1966, Formula 26.6.2.
		///</remarks>
		///<param name="which">indicates which argument is to be calculated
		///    from the others.
		///    1: Calculate P and Q from F, DFN and DFD;
		///    2: Calculate F from P, Q, DFN and DFD;
		///    3: Calculate DFN from P, Q, F and DFD;
		///    4: Calculate DFD from P, Q, F and DFN.
		///</param>
		///<param name="p">The integral from 0 to F of 
		///    the F-density.  If it is an input value, it should lie in the
		///    range [0,1].
		///</param>
		///<param name="q">Equal to 1-P.  If Q is an input
		///    value, it should lie in the range [0,1].  If Q is an output value,
		///    it will lie in the range [0,1].
		///</param>
		///<param name="f">The upper limit of integration 
		///    of the F-density.  If this is an input value, it should lie in the
		///    range [0, +infinity).  If it is an output value, it will be searched
		///    for in the range [0,1.0D+300].
		///</param>
		///<param name="dfn">The number of degrees of 
		///    freedom of the numerator sum of squares.  If this is an input value,
		///    it should lie in the range: (0, +infinity).  If it is an output value,
		///    it will be searched for in the range: [ 1.0D-300, 1.0D+300].
		///</param>
		///<param name="dfd">The number of degrees of freedom 
		///    of the denominator sum of squares.  If this is an input value, it should
		///    lie in the range: (0, +infinity).  If it is an output value, it will
		///    be searched for in the  range: [ 1.0D-300, 1.0D+300].
		///</param>
		///<param name="status">Reports the status of the computation.
		///     0, if the calculation completed correctly;
		///    -I, if the input parameter number I is out of range;
		///    +1, if the answer appears to be lower than lowest search bound;
		///    +2, if the answer appears to be higher than greatest search bound;
		///    +3, if P + Q /= 1.
		///</param>
		///<param name="bound">Is only defined if STATUS is nonzero.
		///    If STATUS is negative, then this is the value exceeded by parameter I.
		///    if STATUS is 1 or 2, this is the search bound that was exceeded.
		///</param>
		void cdff ( int %which, double %p, double %q, double %f, double %dfn, double %dfd, int %status, double %bound );

		void cdffnc ( int %which, double %p, double %q, double %f, double %dfn,
			double %dfd, double %phonc, int %status, double %bound );
		void cdfgam ( int %which, double %p, double %q, double %x, double %shape,
			double %scale, int %status, double %bound );
		void cdfnbn ( int %which, double %p, double %q, double %s, double %xn,
			double %pr, double %ompr, int %status, double %bound );
		void cdfnor ( int %which, double %p, double %q, double %x, double %mean,
			double %sd, int %status, double %bound );
		void cdfpoi ( int %which, double %p, double %q, double %s, double %xlam,
			int %status, double %bound );
		void cdft ( int %which, double %p, double %q, double %t, double %df,
			int %status, double %bound );
		void chi_noncentral_cdf_values ( int %n_data, double %x, double %lambda, 
			int %df, double %cdf );
		void chi_square_cdf_values ( int %n_data, int %a, double %x, double %fx );
		void cumbet ( double %x, double %y, double %a, double %b, double %cum, double %ccum );
		void cumbin ( double %s, double %xn, double %pr, double %ompr,
			double %cum, double %ccum );
		void cumchi ( double %x, double %df, double %cum, double %ccum );
		void cumchn ( double %x, double %df, double %pnonc, double %cum,
			double %ccum );
		void cumf ( double %f, double %dfn, double %dfd, double %cum, double %ccum );
		void cumfnc ( double %f, double %dfn, double %dfd, double %pnonc,
			double %cum, double %ccum );
		void cumgam ( double %x, double %a, double %cum, double %ccum );
		void cumnbn ( double %s, double %xn, double %pr, double %ompr,
			double %cum, double %ccum );
		void cumnor ( double %arg, double %result, double %ccum );
		void cumpoi ( double %s, double %xlam, double %cum, double %ccum );
		void cumt ( double %t, double %df, double %cum, double %ccum );
		double dbetrm ( double %a, double %b );
		double dexpm1 ( double %x );
		double dinvnr ( double %p, double %q );
		void dinvr ( int %status, double %x, double %fx,
			unsigned long %qleft, unsigned long %qhi );
		double dlanor ( double %x );
		void dstinv ( double %zsmall, double %zbig, double %zabsst,
			double %zrelst, double %zstpmu, double %zabsto, double %zrelto );
		double dstrem ( double %z );
		void dstzr ( double %zxlo, double %zxhi, double %zabstl, double %zreltl );
		double dt1 ( double %p, double %q, double %df );
		void dzror ( int %status, double %x, double %fx, double %xlo,
			double %xhi, unsigned long %qleft, unsigned long %qhi );
		void erf_values ( int %n_data, double %x, double %fx );
		double error_f ( double %x );
		double error_fc ( int %ind, double %x );
		double esum ( int %mu, double %x );
		double eval_pol ( array<double> ^a, int %n, double %x );
		double exparg ( int %l );
		void f_cdf_values ( int %n_data, int %a, int %b, double %x, double %fx );
		void f_noncentral_cdf_values ( int %n_data, int %a, int %b, double %lambda, 
			double %x, double %fx );
		double fifdsign ( double mag, double sign );
		long fifidint ( double a );
		long fifmod ( long a, long b );
		double fpser ( double %a, double %b, double %x, double %eps );
		double gam1 ( double %a );
		void gamma_inc ( double %a, double %x, double %ans, double %qans, int %ind );
		void gamma_inc_inv ( double %a, double %x, double %x0, double %p, double %q,
			int %ierr );
		void gamma_inc_values ( int %n_data, double %a, double %x, double %fx );
		double gamma_ln1 ( double %a );
		double gamma_log ( double %a );
		void gamma_rat1 ( double %a, double %x, double %r, double %p, double %q,
			double %eps );
		void gamma_values ( int %n_data, double %x, double %fx );
		double gamma_x ( double %a );
		double gsumln ( double %a, double %b );
		void negative_binomial_cdf_values ( int %n_data, int %f, int %s, double %p, 
			double %cdf );
		void normal_cdf_values ( int %n_data, double %x, double %fx );
		void poisson_cdf_values ( int %n_data, double %a, int %x, double %fx );
		double psi ( double %xx );
		void psi_values ( int %n_data, double %x, double %fx );
		double rcomp ( double %a, double %x );
		double rexp ( double %x );
		double rlog ( double %x );
		double rlog1 ( double %x );
		void student_cdf_values ( int %n_data, int %a, double %x, double %fx );
		double stvaln ( double %p );

protected:
		virtual void E0000 ( int IENTRY, int %status, double %x, double %fx,
			unsigned long %qleft, unsigned long %qhi, double %zabsst,
			double %zabsto, double %zbig, double %zrelst,
			double %zrelto, double %zsmall, double %zstpmu );
		virtual void E0001 ( int IENTRY, int %status, double %x, double %fx,
			double %xlo, double %xhi, unsigned long %qleft,
			unsigned long %qhi, double %zabstl, double %zreltl,
			double %zxhi, double %zxlo );

	};
}
