using System;
using System.Collections.Generic;
using System.Text;
using dcdflib;

namespace dcdflibTest
{
    class Program
    {

        //****************************************************************************80
        //
        //  Purpose:
        //
        //   TEST005 tests BETA_INC and BETA_INC_VALUES.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test005()
        {
            double a = 0;
            double b = 0;
            double ccdf_compute = 0;
            double ccdf_lookup = 0;
            double cdf_compute = 0;
            double cdf_lookup = 0;
            int ierror = 0;
            int n_data = 0;
            double x = 0;
            double y = 0;

            Console.WriteLine();
            Console.WriteLine("TEST005");
            Console.WriteLine("  BETA_INC computes the incomplete Beta ratio.");
            Console.WriteLine("  BETA_INC_CDF_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("    X         Y         A         B         CDF           CDF");
            Console.WriteLine("                                           (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            StatClass calc = new StatClass();

            for (; ; )
            {
                calc.beta_inc_values(ref n_data, ref a, ref b, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                y = 1.0 - x;

                calc.beta_inc(ref a, ref b, ref x, ref y, ref cdf_compute, ref ccdf_compute, ref ierror);

                Object[] o = { x, y, a, b, cdf_lookup, cdf_compute };
                Console.WriteLine("   {0,-10:0.0#######}{1,-10:0.0#######}{2,-10:0.0#######}{3,-10:0.0#######}{4,-14:0.0############}{5,-14:0.0############}", o);
            }

            Console.WriteLine();
            Console.WriteLine("    X         Y         A         B         CDF           CDF");
            Console.WriteLine("                                           (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            for (; ; )
            {
                calc.beta_inc_values(ref n_data, ref a, ref b, ref x, ref cdf_lookup);

                if (n_data == 0)
                {
                    break;
                }

                ccdf_lookup = 1.0 - cdf_lookup;

                y = 1.0 - x;

                calc.beta_inc(ref a, ref b, ref x, ref y, ref cdf_compute, ref ccdf_compute, ref ierror);

                Object[] o = { x, y, a, b, ccdf_lookup, ccdf_compute };
                Console.WriteLine("   {0,-10:0.0#######}{1,-10:0.0#######}{2,-10:0.0#######}{3,-10:0.0#######}{4,-14:0.0############}{5,-14:0.0############}", o);
            }

            Console.WriteLine();
            Console.WriteLine("    X         Y         A         B         CDF         iErr");
            Console.WriteLine("                                          (Computed)");
            Console.WriteLine();

            // I added this test for a specific execution path. Output agress with native version.
            a = 1.567;
            b = 2.201;
            x = 0.8911;
            y = 1.0 - x;

            calc.beta_inc(ref a, ref b, ref x, ref y, ref cdf_compute, ref ccdf_compute, ref ierror);

            Object[] mo = { x, y, a, b, ccdf_compute, ierror };
            Console.WriteLine("   {0,-10:0.0#######}{1,-10:0.0#######}{2,-10:0.0#######}{3,-10:0.0#######}{4,-14:0.0############}   {5}", mo);

            a = 0.0;
            b = 0.5;
            x = 0;
            y = 1;

            calc.beta_inc(ref a, ref b, ref x, ref y, ref cdf_compute, ref ccdf_compute, ref ierror);

            mo[0] = x;
            mo[1] = y;
            mo[2] = a;
            mo[3] = b;
            mo[4] = ccdf_compute;
            mo[5] = ierror;

            Console.WriteLine("   {0,-10:0.0#######}{1,-10:0.0#######}{2,-10:0.0#######}{3,-10:0.0#######}{4,-14:0.0############}   {5}", mo);

            a = 0.0;
            b = 0.5;
            x = 1;
            y = 0;

            calc.beta_inc(ref a, ref b, ref x, ref y, ref cdf_compute, ref ccdf_compute, ref ierror);

            mo[0] = x;
            mo[1] = y;
            mo[2] = a;
            mo[3] = b;
            mo[4] = ccdf_compute;
            mo[5] = ierror;

            Console.WriteLine("   {0,-10:0.0#######}{1,-10:0.0#######}{2,-10:0.0#######}{3,-10:0.0#######}{4,-14:0.0############}   {5}", mo);

            a = 0.5;
            b = 0.0;
            x = 0;
            y = 1;

            calc.beta_inc(ref a, ref b, ref x, ref y, ref cdf_compute, ref ccdf_compute, ref ierror);

            mo[0] = x;
            mo[1] = y;
            mo[2] = a;
            mo[3] = b;
            mo[4] = ccdf_compute;
            mo[5] = ierror;

            Console.WriteLine("   {0,-10:0.0#######}{1,-10:0.0#######}{2,-10:0.0#######}{3,-10:0.0#######}{4,-14:0.0############}   {5}", mo);

            a = 0.5;
            b = 0.0;
            x = 1;
            y = 0;

            calc.beta_inc(ref a, ref b, ref x, ref y, ref cdf_compute, ref ccdf_compute, ref ierror);

            mo[0] = x;
            mo[1] = y;
            mo[2] = a;
            mo[3] = b;
            mo[4] = ccdf_compute;
            mo[5] = ierror;

            Console.WriteLine("   {0,-10:0.0#######}{1,-10:0.0#######}{2,-10:0.0#######}{3,-10:0.0#######}{4,-14:0.0############}   {5}", mo);

            // Test beta_log
            a = 7.99;
            b = 7.99;
            x = calc.beta_log(ref a, ref b);

            mo[0] = "";
            mo[1] = "";
            mo[2] = a;
            mo[3] = b;
            mo[4] = x;
            mo[5] = 0;

            Console.WriteLine("   {0,-10:0.0#######}{1,-10:0.0#######}{2,-10:0.0#######}{3,-10:0.0#######}{4,-14:0.0############}   {5}", mo);
            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests CDFBET.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test01()
        {
            double a = 0;
            double b = 0;
            double bound = 0;
            double p = 0;
            double q = 0;
            int status = 0;
            int which = 0;
            double x = 0;
            double y = 0;

            Console.WriteLine();
            Console.WriteLine("TEST01");
            Console.WriteLine("  CDFBET computes one missing parameter from the");
            Console.WriteLine("    BETA CDF:");
            Console.WriteLine();
            Console.WriteLine("   BETA_CDF ( (P,Q), (X,Y), A, B )");
            Console.WriteLine();
            Console.WriteLine("      P           Q           X           Y           A           B");
            Console.WriteLine();

            for (which = 1; which <= 4; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    x = 0.25;
                    y = 1.0 - x;
                    a = 2.0;
                    b = 3.0;
                }
                else if (which == 2)
                {
                    p = 0.261719;
                    q = 1.0 - p;
                    x = -1.0;
                    y = -1.0;
                    a = 2.0;
                    b = 3.0;
                }
                else if (which == 3)
                {
                    p = 0.261719;
                    q = 1.0 - p;
                    x = 0.25;
                    y = 1.0 - x;
                    a = -1.0;
                    b = 3.0;
                }
                else if (which == 4)
                {
                    p = 0.261719;
                    q = 1.0 - p;
                    x = 0.25;
                    y = 1.0 - x;
                    a = 2.0;
                    b = -1.0;
                }

                StatClass calc = new StatClass();

                calc.cdfbet(ref which, ref p, ref q, ref x, ref y, ref a, ref b, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine();
                    Console.WriteLine("  CDFBET returned STATUS = {0}", status);
                    return;
                }

                Object[] o = { p, q, x, y, a, b };
                Console.WriteLine("  {0,-10:0.0#######}  {1,-10:0.0#######}    {2,-10:0.0#######}  {3,-10:0.0#######}  {4,-10:0.0#######}  {5,-10:0.0#######}", o);
            }

            return;

        }
        //****************************************************************************80


        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests CDFBIN.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test02()
        {
            double bound = 0;
            double ompr = 0;
            double p = 0;
            double pr = 0;
            double q = 0;
            double s = 0;
            int status = 0;
            int which = 0;
            double xn = 0;

            Console.WriteLine();
            Console.WriteLine("TEST02");
            Console.WriteLine("  CDFBIN computes one missing parameter from the");
            Console.WriteLine("    Binomial CDF:");
            Console.WriteLine();
            Console.WriteLine("   BINOMIAL_CDF ( (P,Q), S, XN, (PR,OMPR) )");
            Console.WriteLine();
            Console.WriteLine("      P           Q           S          XN               PR             OMPR");
            Console.WriteLine();

            StatClass calc = new StatClass();

            for (which = 1; which <= 4; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    s = 5.0;
                    xn = 8.0;
                    pr = 0.875;
                    ompr = 1.0 - pr;
                }
                else if (which == 2)
                {
                    p = 0.067347;
                    q = 1.0 - p;
                    s = -1.0;
                    xn = 8.0;
                    pr = 0.875;
                    ompr = 1.0 - pr;
                }
                else if (which == 3)
                {
                    p = 0.067347;
                    q = 1.0 - p;
                    s = 5.0;
                    xn = -1.0;
                    pr = 0.875;
                    ompr = 1.0 - pr;
                }
                else if (which == 4)
                {
                    p = 0.067347;
                    q = 1.0 - p;
                    s = 5.0;
                    xn = 8.0;
                    pr = -1.0;
                    ompr = -1.0;
                }

                calc.cdfbin(ref which, ref p, ref q, ref s, ref xn, ref pr, ref ompr, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine();
                    Console.WriteLine("  CDFBIN returned STATUS = {0}", status);
                    continue;
                }
                Object[] o = { p, q, s, xn, pr, ompr };
                Console.WriteLine(" {0,10:0.000000}  {1,10:0.000000}  {2,10:0.000000}  {3,10:0.000000}  {4,14:0.000000}  {5,14:0.000000}", o);
            }

            return;
        }
        //****************************************************************************80


        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests CDFCHI.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test03()
        {
            double bound = 0;
            double df = 0;
            double p = 0;
            double q = 0;
            int status = 0;
            int which = 0;
            double x = 0;

            Console.WriteLine();
            Console.WriteLine("TEST03");
            Console.WriteLine("  CDFCHI computes one missing parameter from the");
            Console.WriteLine("    Chi Square CDF:");
            Console.WriteLine();
            Console.WriteLine("   CHI_CDF ( (P,Q), X, DF )");
            Console.WriteLine();
            Console.WriteLine("      P           Q           X          DF");
            Console.WriteLine();

            StatClass calc = new StatClass();

            for (which = 1; which <= 3; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    x = 5.0;
                    df = 8.0;
                }
                else if (which == 2)
                {
                    p = 0.242424;
                    q = 1.0 - p;
                    x = -1.0;
                    df = 8.0;
                }
                else if (which == 3)
                {
                    p = 0.242424;
                    q = 1.0 - p;
                    x = 5.0;
                    df = -1.0;
                }

                calc.cdfchi(ref which, ref p, ref q, ref x, ref df, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine();
                    Console.WriteLine("  CDFCHI returned STATUS = {0}", status);
                    return;
                }

                Object[] o = { p, q, x, df };
                Console.WriteLine(" {0,10:0.000000}  {1,10:0.000000}  {2,10:0.000000}  {3,10:0.000000}", o);
            }
            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests CDFCHN.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test04()
        {
            double bound = 0;
            double df = 0;
            double p = 0;
            double pnonc = 0;
            double q = 0;
            int status = 0;
            int which = 0;
            double x = 0;

            Console.WriteLine();
            Console.WriteLine("TEST04");
            Console.WriteLine("  CDFCHN computes one missing parameter from the");
            Console.WriteLine("    Chi Square CDF:");
            Console.WriteLine();
            Console.WriteLine("   CHI_Noncentral_CDF ( (P,Q), X, DF, PNONC )");
            Console.WriteLine();
            Console.WriteLine("      P           Q           X          DF             PNONC");
            Console.WriteLine();

            StatClass calc = new StatClass();

            for (which = 1; which <= 4; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    x = 5.0;
                    df = 8.0;
                    pnonc = 0.5;
                }
                else if (which == 2)
                {
                    p = 0.211040;
                    q = 1.0 - p;
                    x = -1.0;
                    df = 8.0;
                    pnonc = 0.5;
                }
                else if (which == 3)
                {
                    p = 0.211040;
                    q = 1.0 - p;
                    x = 5.0;
                    df = -1.0;
                    pnonc = 0.5;
                }
                else if (which == 4)
                {
                    p = 0.211040;
                    q = 1.0 - p;
                    x = 5.0;
                    df = 8.0;
                    pnonc = -1.0;
                }

                calc.cdfchn(ref which, ref p, ref q, ref x, ref df, ref pnonc, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine();
                    Console.WriteLine("  CDFCHN returned STATUS = {0}", status);
                    return;
                }

                Object[] o = { p, q, x, df, pnonc };
                Console.WriteLine(" {0,10:0.000000}  {1,10:0.000000}  {2,10:0.000000}  {3,10:0.000000}  {4,14:0.000000}", o);
            }
            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests CDFF.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test05()
        {
            double bound = 0;
            double dfd = 0;
            double dfn = 0;
            double f = 0;
            double p = 0;
            double q = 0;
            int status = 0;
            int which = 0;

            Console.WriteLine();
            Console.WriteLine("TEST05");
            Console.WriteLine("  CDFF computes one missing parameter from the");
            Console.WriteLine("    F CDF:");
            Console.WriteLine();
            Console.WriteLine("   F_CDF ( (P,Q), F, DFN, DFD )");
            Console.WriteLine();
            Console.WriteLine("         P         Q           F         DFN         DFD");
            Console.WriteLine();

            StatClass calc = new StatClass();

            for (which = 1; which <= 4; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    f = 5.0;
                    dfn = 8.0;
                    dfd = 3.0;
                }
                else if (which == 2)
                {
                    p = 0.893510;
                    q = 1.0 - p;
                    f = -1.0;
                    dfn = 8.0;
                    dfd = 3.0;
                }
                else if (which == 3)
                {
                    p = 0.893510;
                    q = 1.0 - p;
                    f = 5.0;
                    dfn = -1.0;
                    dfd = 3.0;
                }
                else if (which == 4)
                {
                    p = 0.893510;
                    q = 1.0 - p;
                    f = 5.0;
                    dfn = 8.0;
                    dfd = -1.0;
                }

                calc.cdff(ref which, ref p, ref q, ref f, ref dfn, ref dfd, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine();
                    Console.WriteLine("  CDFF returned STATUS = {0}", status);
                    continue;
                }

                Object[] o = { p, q, f, dfn, dfd };
                Console.WriteLine(" {0,10:0.000000}  {1,10:0.000000}  {2,10:0.000000}  {3,10:0.000000}  {4,10:0.000000}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests CDFFNC.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test06()
        {
            double bound = 0;
            double dfd = 0;
            double dfn = 0;
            double f = 0;
            double p = 0;
            double pnonc = 0;
            double q = 0;
            int status = 0;
            int which = 0;

            Console.WriteLine();
            Console.WriteLine("TEST06");
            Console.WriteLine("  CDFFNC computes one missing parameter from the");
            Console.WriteLine("    noncentral F CDF:");
            Console.WriteLine();
            Console.WriteLine("   F_noncentral_CDF ( (P,Q), F, DFN, DFD, PNONC )");
            Console.WriteLine();
            Console.WriteLine("       P        Q           F         DFN         DFD        PNONC");
            Console.WriteLine();

            StatClass calc = new StatClass();

            for (which = 1; which <= 5; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    f = 5.0;
                    dfn = 8.0;
                    dfd = 3.0;
                    pnonc = 17.648016;
                }
                else if (which == 2)
                {
                    p = 0.60;
                    q = 1.0 - p;
                    f = -1.0;
                    dfn = 8.0;
                    dfd = 3.0;
                    pnonc = 17.648016;
                }
                else if (which == 3)
                {
                    p = 0.60;
                    q = 1.0 - p;
                    f = 5.0;
                    dfn = -1.0;
                    dfd = 3.0;
                    pnonc = 17.648016;
                }
                else if (which == 4)
                {
                    p = 0.60;
                    q = 1.0 - p;
                    f = 5.0;
                    dfn = 8.0;
                    dfd = -1.0;
                    pnonc = 17.648016;
                }
                else if (which == 5)
                {
                    p = 0.60;
                    q = 1.0 - p;
                    f = 5.0;
                    dfn = 8.0;
                    dfd = 3.0;
                    pnonc = -1.0;
                }

                calc.cdffnc(ref which, ref p, ref q, ref f, ref dfn, ref dfd, ref pnonc, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine();
                    Console.WriteLine("  CDFFNC returned STATUS = {0}", status);
                    continue;
                }

                Object[] o = { p, q, f, dfn, dfd, pnonc };
                Console.WriteLine(" {0,10:0.000000}  {1,10:0.000000}  {2,10:0.000000}  {3,10:0.000000}  {4,10:0.000000}  {5,10:0.000000}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests CDFGAM.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test07()
        {
            double bound = 0;
            double p = 0;
            double q = 0;
            double scale = 0;
            double shape = 0;
            int status = 0;
            int which = 0;
            double x = 0;

            Console.WriteLine();
            Console.WriteLine("TEST07");
            Console.WriteLine("  CDFGAM computes one missing parameter from the");
            Console.WriteLine("    Gamma CDF:");
            Console.WriteLine();
            Console.WriteLine("   Gamma_CDF ( (P,Q), X, SHAPE, SCALE )");
            Console.WriteLine();
            Console.WriteLine("      P           Q            X        SHAPE       SCALE");
            Console.WriteLine();

            StatClass calc = new StatClass();

            for (which = 1; which <= 4; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    x = 5.0;
                    shape = 8.0;
                    scale = 3.0;
                }
                else if (which == 2)
                {
                    p = 0.981998;
                    q = 1.0 - p;
                    x = -1.0;
                    shape = 8.0;
                    scale = 3.0;
                }
                else if (which == 3)
                {
                    p = 0.981998;
                    q = 1.0 - p;
                    x = 5.0;
                    shape = -1.0;
                    scale = 3.0;
                }
                else if (which == 4)
                {
                    p = 0.981998;
                    q = 1.0 - p;
                    x = 5.0;
                    shape = 8.0;
                    scale = -1.0;
                }

                calc.cdfgam(ref which, ref p, ref q, ref x, ref shape, ref scale, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine();
                    Console.WriteLine("  CDFGAM returned STATUS = {0}", status);
                    continue;
                }

                Object[] o = { p, q, x, shape, scale };
                Console.WriteLine(" {0,10:0.000000}  {1,10:0.000000}  {2,10:0.000000}  {3,10:0.000000}  {4,10:0.000000}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests CDFNBN.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test08()
        {
            double bound = 0;
            double f = 0;
            double ompr = 0;
            double p = 0;
            double pr = 0;
            double q = 0;
            double s = 0;
            int status = 0;
            int which = 0;

            Console.WriteLine();
            Console.WriteLine("TEST08");
            Console.WriteLine("  CDFNBN computes one missing parameter from the");
            Console.WriteLine("    Negative_Binomial CDF:");
            Console.WriteLine();
            Console.WriteLine("   Negative_BINOMIAL_CDF ( (P,Q), F, S, (PR,OMPR) )");
            Console.WriteLine();
            Console.WriteLine("         P         Q           F           S         PR        OMPR");
            Console.WriteLine();

            StatClass calc = new StatClass();

            for (which = 1; which <= 4; which++)
            {
                if (which == 1)
                {
                    p = -1.0;
                    q = -1.0;
                    f = 3.0;
                    s = 5.0;
                    pr = 0.875;
                    ompr = 1.0 - pr;
                }
                else if (which == 2)
                {
                    p = 0.988752;
                    q = 1.0 - p;
                    f = -1.0;
                    s = 5.0;
                    pr = 0.875;
                    ompr = 1.0 - pr;
                }
                else if (which == 3)
                {
                    p = 0.988752;
                    q = 1.0 - p;
                    f = 3.0;
                    s = -1.0;
                    pr = 0.875;
                    ompr = 1.0 - pr;
                }
                else if (which == 4)
                {
                    p = 0.988752;
                    q = 1.0 - p;
                    f = 3.0;
                    s = 5.0;
                    pr = -1.0;
                    ompr = -1.0;
                }

                calc.cdfnbn(ref which, ref p, ref q, ref f, ref s, ref pr, ref ompr, ref status, ref bound);

                if (status != 0)
                {
                    Console.WriteLine();
                    Console.WriteLine("  CDFNBN returned STATUS = {0}", status);
                    continue;
                }

                Object[] o = { p, q, f, s, pr, ompr };
                Console.WriteLine(" {0,10:0.000000}  {1,10:0.000000}  {2,10:0.000000}  {3,10:0.000000}  {4,10:0.000000}  {5,10:0.000000}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests CDFNOR.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test09()
        {
            double bound = 0.0;
            double mean = 0.0;
            double p = 0.0;
            double q = 0.0;
            double sd = 0.0;
            int status = 0;
            int which = 0;
            double x = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST09");
            Console.WriteLine("  CDFNOR computes one missing parameter from the");
            Console.WriteLine("    Normal CDF:");
            Console.WriteLine();
            Console.WriteLine("   Normal_CDF ( (P,Q), X, MEAN, SD )");
            Console.WriteLine();
            Console.WriteLine("         P         Q           X          MEAN      SD");
            Console.WriteLine();

            StatClass calc = new StatClass();
            
            for ( which = 1; which <= 4; which++ )
            {
                if ( which == 1 )
                {
                    p = -1.0;
                    q = -1.0;
                    x = 3.0;
                    mean = 5.0;
                    sd = 0.875;
                }
                else if ( which == 2 )
                {
                    p = 0.011135;
                    q = 1.0 - p;
                    x = -1.0;
                    mean = 5.0;
                    sd = 0.875;
                }
                else if ( which == 3 )
                {
                    p = 0.011135;
                    q = 1.0 - p;
                    x = 3.0;
                    mean = -1.0;
                    sd = 0.875;
                }
                else if ( which == 4 )
                {
                    p = 0.011135;
                    q = 1.0 - p;
                    x = 3.0;
                    mean = 5.0;
                    sd = -1.0;
                }

                calc.cdfnor ( ref which, ref p, ref q, ref x, ref mean, ref sd, ref status, ref bound );

                if ( status != 0 )
                {
                    Console.WriteLine();
                    Console.WriteLine("  CDFNOR returned STATUS = {0}", status);
                }
                 
                Object[] o = { p, q, x, mean, sd };
                Console.WriteLine(" {0,10:0.000000}  {1,10:0.000000}  {2,10:0.000000}  {3,10:0.000000}  {4,10:0.000000}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests CDFPOI.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test10()
        {
            double bound = 0.0;
            double p = 0.0;
            double q = 0.0;
            double s = 0.0;
            int status = 0;
            int which = 0;
            double xlam = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST10");
            Console.WriteLine("  CDFPOI computes one missing parameter from the");
            Console.WriteLine("    Poisson CDF:");
            Console.WriteLine();
            Console.WriteLine("   POISSON_CDF ( (P,Q), S, XLAM )");
            Console.WriteLine();
            Console.WriteLine("         P         Q           S          XLAM");
            Console.WriteLine();

            StatClass calc = new StatClass();
            
            for ( which = 1; which <= 3; which++ )
            {
                if ( which == 1 )
                {
                    p = -1.0;
                    q = -1.0;
                    s = 3.0;
                    xlam = 5.0;
                }
                else if ( which == 2 )
                {
                    p = 0.265026;
                    q = 1.0 - p;
                    s = -1.0;
                    xlam = 5.0;
                }
                else if ( which == 3 )
                {
                    p = 0.265026;
                    q = 1.0 - p;
                    s = 3.0;
                    xlam = -1.0;
                }

                calc.cdfpoi ( ref which, ref p, ref q, ref s, ref xlam, ref status, ref bound );

                if ( status != 0 )
                {
                    Console.WriteLine();
                    Console.WriteLine("  CDFPOI returned STATUS = {0}", status);
                    continue;
                }
                Object[] o = { p, q, s, xlam };
                Console.WriteLine(" {0,10:0.000000}  {1,10:0.000000}  {2,10:0.000000}  {3,10:0.000000}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests CDFT.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test11()
        {
            double bound = 0.0;
            double df = 0.0;
            double p = 0.0;
            double q = 0.0;
            int status = 0;
            double t = 0.0;
            int which = 0;

            Console.WriteLine();
            Console.WriteLine("TEST11");
            Console.WriteLine("  CDFT computes one missing parameter from the");
            Console.WriteLine("    T CDF:");
            Console.WriteLine();
            Console.WriteLine("   T_CDF ( (P,Q), T, DF )");
            Console.WriteLine();
            Console.WriteLine("         P         Q           T          DF");
            Console.WriteLine();

            StatClass calc = new StatClass();
            
            for ( which = 1; which <= 3; which++ )
            {
                if ( which == 1 )
                {
                    p = -1.0;
                    q = -1.0;
                    t = 3.0;
                    df = 5.0;
                }
                else if ( which == 2 )
                {
                    p = 0.984950;
                    q = 1.0 - p;
                    t = -1.0;
                    df = 5.0;
                }
                else if ( which == 3 )
                {
                    p = 0.984950;
                    q = 1.0 - p;
                    t = 3.0;
                    df = -1.0;
                }

                calc.cdft ( ref which, ref p, ref q, ref t, ref df, ref status, ref bound );

                if ( status != 0 )
                {
                    Console.WriteLine();
                    Console.WriteLine("  CDFT returned STATUS = {0}", status);
                }
                Object[] o = { p, q, t, df };
                Console.WriteLine(" {0,10:0.000000}  {1,10:0.000000}  {2,10:0.000000}  {3,10:0.000000}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tests CUMBET, BETA_INC_VALUES.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test12()
        {
            double a = 0.0;
            double b = 0.0;
            double ccdf_compute = 0.0;
            double ccdf_lookup = 0.0;
            double cdf_compute = 0.0;
            double cdf_lookup = 0.0;
            int n_data = 0;
            double x = 0.0;
            double y = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST12");
            Console.WriteLine("  CUMBET computes the Beta CDF");
            Console.WriteLine("    and the complementary CDF.");
            Console.WriteLine("  BETA_INC_CDF_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("         X         Y           A           B        CDF        CDF");
            Console.WriteLine("                                                  (Lookup)   (Computed)");
            Console.WriteLine();

            StatClass calc = new StatClass();
            n_data = 0;

            for ( ; ; )
            {
                calc.beta_inc_values ( ref n_data, ref a, ref b, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                y = 1.0 - x;

                calc.cumbet ( ref x, ref y, ref a, ref b, ref cdf_compute, ref ccdf_compute );

                Object[] o = { x, y, a, b, cdf_lookup, cdf_compute };
                Console.WriteLine(" {0,10:0.000000}  {1,10:0.000000}  {2,10:0.000000}  {3,10:0.000000}  {4,10:0.000000}  {5,10:0.000000}", o);
            }

            Console.WriteLine();
            Console.WriteLine("         X         Y           A           B        CDF        CDF");
            Console.WriteLine("                                                  (Lookup)   (Computed)");
            Console.WriteLine();

            n_data = 0;

            for ( ; ; )
            {
                calc.beta_inc_values ( ref n_data, ref a, ref b, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                ccdf_lookup = 1.0 - cdf_lookup;

                y = 1.0 - x;

                calc.cumbet ( ref x, ref y, ref a, ref b, ref cdf_compute, ref ccdf_compute );

                Object[] o = { x, y, a, b, ccdf_lookup, ccdf_compute };
                Console.WriteLine(" {0,10:0.000000}  {1,10:0.000000}  {2,10:0.000000}  {3,10:0.000000}  {4,10:0.000000}  {5,10:0.000000}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST13 tests CUMBIN, BINOMIAL_CDF_VALUES.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test13()
        {
            double ccdf_compute = 0.0;
            double ccdf_lookup = 0.0;
            double cdf_compute = 0.0;
            double cdf_lookup = 0.0;
            int n_data = 0;
            double ompr = 0.0;
            int s = 0;
            double s_double = 0.0;
            double pr = 0.0;
            int x = 0;
            double x_double = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST13");
            Console.WriteLine("  CUMBIN computes the Binomial CDF");
            Console.WriteLine("    and the complementary CDF.");
            Console.WriteLine("  BINOMIAL_CDF_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("   X   S    Pr       CDF           CDF");
            Console.WriteLine("                    (Lookup)      (Computed)");
            Console.WriteLine();

            StatClass calc = new StatClass();

            n_data = 0;

            for ( ; ; )
            {
                calc.binomial_cdf_values(ref n_data, ref x, ref pr, ref s, ref cdf_lookup);

                if ( n_data == 0 )
                {
                    break;
                }

                ompr = 1.0 - pr;

                s_double = ( double ) s;
                x_double = ( double ) x;

                calc.cumbin ( ref s_double, ref x_double, ref pr, ref ompr, ref cdf_compute, ref ccdf_compute );

                Object[] o = { s, x, pr, cdf_lookup, cdf_compute };
                Console.WriteLine("  {0,2}  {1,2}  {2,-8:0.000000}  {3,-12:0.0#########}  {4,-12:0.0#########}", o);
            }

            Console.WriteLine();
            Console.WriteLine("   X   S    Pr       1-CDF         CCDF");
            Console.WriteLine("                    (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            for ( ; ; )
            {
                calc.binomial_cdf_values ( ref n_data, ref x, ref pr, ref s, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                ccdf_lookup = 1.0 - cdf_lookup;

                ompr = 1.0 - pr;

                s_double = ( double ) s;
                x_double = ( double ) x;

                calc.cumbin ( ref s_double, ref x_double, ref pr, ref ompr, ref cdf_compute, ref ccdf_compute );

                Object[] o = { s, x, pr, ccdf_lookup, ccdf_compute };
                Console.WriteLine("  {0,2}  {1,2}   {2,-8:0.0#####} {3,-12:0.0#########}  {4,-12:0.0#########}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST14 tests CUMCHI, CHI_SQUARE_CDF_VALUES.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test14()
        {
            double ccdf_compute = 0.0;
            double ccdf_lookup = 0.0;
            double cdf_compute = 0.0;
            double cdf_lookup = 0.0;
            int df = 0;
            double df_double = 0.0;
            int n_data = 0;
            double x = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST14");
            Console.WriteLine("  CUMCHI computes the chi square CDF");
            Console.WriteLine("    and the complementary CDF.");
            Console.WriteLine("  CHI_SQUARE_CDF_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("    X       DF    CDF           CDF");
            Console.WriteLine("                 (Lookup)      (Computed)");
            Console.WriteLine();

            StatClass calc = new StatClass();

            n_data = 0;

            for ( ; ; )
            {
                calc.chi_square_cdf_values ( ref n_data, ref df, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                df_double = ( double ) df;

                calc.cumchi ( ref x, ref df_double, ref cdf_compute, ref ccdf_compute );

                Object[] o = { x, df, cdf_lookup, cdf_compute };
                Console.WriteLine("  {0,-8:0.0#####}  {1,2}  {2,-12:0.0#########}  {3,-12:0.0#########}", o);
            }

            Console.WriteLine();
            Console.WriteLine("    X       DF    1-CDF         CCDF");
            Console.WriteLine("                 (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            for ( ; ; )
            {
                calc.chi_square_cdf_values ( ref n_data, ref df, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                ccdf_lookup = 1.0 - cdf_lookup;

                df_double = ( double ) df;

                calc.cumchi ( ref x, ref df_double, ref cdf_compute, ref ccdf_compute );

                Object[] o = { x, df, ccdf_lookup, ccdf_compute };
                Console.WriteLine("  {0,-8:0.0#####}  {1,2}  {2,-12:0.0#########}  {3,-12:0.0#########}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST15 tests CUMCHN, CHI_NONCENTRAL_CDF_VALUES.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test15()
        {
            double ccdf_compute = 0.0;
            double ccdf_lookup = 0.0;
            double cdf_compute = 0.0;
            double cdf_lookup = 0.0;
            int df = 0;
            double df_double = 0.0;
            double lambda = 0.0;
            int n_data = 0;
            double x = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST15");
            Console.WriteLine("  CUMCHN computes the cumulative density");
            Console.WriteLine("    function for the noncentral chi-squared");
            Console.WriteLine("    distribution.");
            Console.WriteLine("  CHI_NONCENTRAL_CDF_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("    DF    Lambda    X         CDF           CDF");
            Console.WriteLine("                             (Lookup)      (Computed)");
            Console.WriteLine();

            StatClass calc = new StatClass();

            n_data = 0;

            for ( ; ; )
            {
                calc.chi_noncentral_cdf_values ( ref n_data, ref x, ref lambda, ref df, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                df_double = ( double ) df;

                calc.cumchn(ref x, ref df_double, ref lambda, ref cdf_compute, ref ccdf_compute);

                Object[] o = { df, lambda, x, cdf_lookup, cdf_compute };
                Console.WriteLine(" {0,6}   {1,-8:0.0#####} {2,-8:0.0#####}  {3,-12:0.0#########}  {4,-12:0.0#########}", o);
            }

            Console.WriteLine();
            Console.WriteLine("    DF    Lambda    X         1-CDF         CCDF");
            Console.WriteLine("                             (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            for ( ; ; )
            {
                calc.chi_noncentral_cdf_values ( ref n_data, ref x, ref lambda, ref df, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                ccdf_lookup = 1.0 - cdf_lookup;

                df_double = ( double ) df;

                calc.cumchn ( ref x, ref df_double, ref lambda, ref cdf_compute, ref ccdf_compute );

                Object[] o = { df, lambda, x, ccdf_lookup, ccdf_compute };
                Console.WriteLine(" {0,6}   {1,-8:0.0#####} {2,-8:0.0#####}  {3,-12:0.0#########}  {4,-12:0.0#########}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST16 tests CUMF, F_CDF_VALUES.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test16()
        {
            double ccdf_compute = 0.0;
            double ccdf_lookup = 0.0;
            double cdf_compute = 0.0;
            double cdf_lookup = 0.0;
            int dfd = 0;
            double dfd_double = 0.0;
            int dfn = 0;
            double dfn_double = 0.0;
            int n_data = 0;
            double x = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST16");
            Console.WriteLine("  CUMF computes the F CDF");
            Console.WriteLine("    and the complementary CDF.");
            Console.WriteLine("  F_CDF_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("    X      DFN DFD    CDF            CDF");
            Console.WriteLine("                     (Lookup)      (Computed)");
            Console.WriteLine();

            StatClass calc = new StatClass();

            n_data = 0;

            for ( ; ; )
            {
                calc.f_cdf_values ( ref n_data, ref dfn, ref dfd, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                dfn_double = ( double ) dfn;
                dfd_double = ( double ) dfd;

                calc.cumf ( ref x, ref dfn_double, ref dfd_double, ref cdf_compute, ref ccdf_compute );

                Object[] o = { x, dfn, dfd, cdf_lookup, cdf_compute };
                Console.WriteLine("  {0,-8:0.0#####} {1,2}  {2,2}     {3,-12:0.0#########}  {4,-12:0.0#########}", o);

            }

            Console.WriteLine();
            Console.WriteLine("    X      DFN DFD    1-CDF         CCDF");
            Console.WriteLine("                     (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            for ( ; ; )
            {
                calc.f_cdf_values ( ref n_data, ref dfn, ref dfd, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                ccdf_lookup = 1.0 - cdf_lookup;

                dfn_double = ( double ) dfn;
                dfd_double = ( double ) dfd;

                calc.cumf ( ref x, ref dfn_double, ref dfd_double, ref cdf_compute, ref ccdf_compute );

                Object[] o = { x, dfn, dfd, ccdf_lookup, ccdf_compute };
                Console.WriteLine("  {0,-8:0.0#####} {1,2}  {2,2}     {3,-12:0.0#########}  {4,-12:0.0#########}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST17 tests CUMFNC, F_NONCENTRAL_CDF_VALUES.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test17()
        {
            double ccdf_compute = 0.0;
            double ccdf_lookup = 0.0;
            double cdf_compute = 0.0;
            double cdf_lookup = 0.0;
            int dfd = 0;
            double dfd_double = 0.0;
            int dfn = 0;
            double dfn_double = 0.0;
            double lambda = 0.0;
            int n_data = 0;
            double x = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST17");
            Console.WriteLine("  CUMFNC computes the noncentral F CDF");
            Console.WriteLine("    and the complementary CDF.");
            Console.WriteLine("  F_NONCENTRAL_CDF_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("    X      DFN DFD    LAMBDA      CDF           CDF");
            Console.WriteLine("                                 (Lookup)      (Computed)");
            Console.WriteLine();

            StatClass calc = new StatClass();

            n_data = 0;

            for ( ; ; )
            {
                calc.f_noncentral_cdf_values ( ref n_data, ref dfn, ref dfd, ref lambda, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                dfn_double = ( double ) dfn;
                dfd_double = ( double ) dfd;

                calc.cumfnc ( ref x, ref dfn_double, ref dfd_double, ref lambda, ref cdf_compute, ref ccdf_compute );

                Object[] o = { x, dfn, dfd, lambda, cdf_lookup, cdf_compute };
                Console.WriteLine("  {0,-8:0.0#####} {1,2}  {2,2}      {3,-8:0.0#####} {4,-12:0.0#########}  {5,-12:0.0#########}", o);
            }

            Console.WriteLine();
            Console.WriteLine("    X      DFN DFD    LAMBDA      1-CDF         CCDF");
            Console.WriteLine("                                 (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            for ( ; ; )
            {
                calc.f_noncentral_cdf_values ( ref n_data, ref dfn, ref dfd, ref lambda, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                ccdf_lookup = 1.0 - cdf_lookup;

                dfn_double = ( double ) dfn;
                dfd_double = ( double ) dfd;

                calc.cumfnc ( ref x, ref dfn_double, ref dfd_double, ref lambda, ref cdf_compute, ref ccdf_compute );

                Object[] o = { x, dfn, dfd, lambda, ccdf_lookup, ccdf_compute };
                Console.WriteLine("  {0,-8:0.0#####} {1,2}  {2,2}      {3,-8:0.0#####} {4,-12:0.0#########}  {5,-12:0.0#########}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST18 tests CUMGAM, GAMMA_INC_VALUES.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test18()
        {
            double a = 0.0;
            double ccdf_compute = 0.0;
            double ccdf_lookup = 0.0;
            double cdf_compute = 0.0;
            double cdf_lookup = 0.0;
            int n_data = 0;
            double x = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST18");
            Console.WriteLine("  CUMGAM computes the Gamma CDF");
            Console.WriteLine("    and the complementary CDF.");
            Console.WriteLine("  GAMMA_INC_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("    A         X           CDF           CDF");
            Console.WriteLine("                        (Lookup)      (Computed)");
            Console.WriteLine();

            StatClass calc = new StatClass();

            n_data = 0;

            for ( ; ; )
            {
                calc.gamma_inc_values ( ref n_data, ref a, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                calc.cumgam ( ref x, ref a, ref cdf_compute, ref ccdf_compute );

                Object[] o = { a, x, cdf_lookup, cdf_compute };
                Console.WriteLine("  {0,-9:#0.0####} {1,-9:#0.0####}  {2,-12:0.0#########}  {3,-12:0.0#########}", o);
            }

            Console.WriteLine();
            Console.WriteLine("    A         X           CDF           CDF");
            Console.WriteLine("                        (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            for ( ; ; )
            {
                calc.gamma_inc_values ( ref n_data, ref a, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                ccdf_lookup = 1.0 - cdf_lookup;

                calc.cumgam ( ref x, ref a, ref cdf_compute, ref ccdf_compute );

                Object[] o = { a, x, ccdf_lookup, ccdf_compute };
                Console.WriteLine("  {0,-9:#0.0####} {1,-9:#0.0####}  {2,-12:0.0#########}  {3,-12:0.0#########}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST19 tests CUMNBN, NEGATIVE_BINOMIAL_CDF_VALUES.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test19()
        {
            double ccdf_compute = 0.0;
            double ccdf_lookup = 0.0;
            double cdf_compute = 0.0;
            double cdf_lookup = 0.0;
            int f = 0;
            double f_double = 0.0;
            int n_data = 0;
            double ompr = 0.0;
            int s = 0;
            double s_double = 0.0;
            double pr = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST19");
            Console.WriteLine("  CUMNBN computes the Negative Binomial CDF");
            Console.WriteLine("    and the complementary CDF.");
            Console.WriteLine("  NEGATIVE_BINOMIAL_CDF_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("   F   S    Pr        CDF           CDF");
            Console.WriteLine("                      (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            StatClass calc = new StatClass();

            for ( ; ; )
            {
                calc.negative_binomial_cdf_values ( ref n_data, ref f, ref s, ref pr, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                ompr = 1.0 - pr;

                f_double = ( double ) f;
                s_double = ( double ) s;

                calc.cumnbn ( ref f_double, ref s_double, ref pr, ref ompr, ref cdf_compute, ref ccdf_compute );

                Object[] o = { f, s, pr, cdf_lookup, cdf_compute };
                Console.WriteLine("  {0,2}  {1,2}   {2,-8:0.0#####}  {3,-12:0.0#########}  {4,-12:0.0#########}", o);
            }

            Console.WriteLine();
            Console.WriteLine("   F   S    Pr        1-CDF         CCDF");
            Console.WriteLine("                      (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            for ( ; ; )
            {
                calc.negative_binomial_cdf_values ( ref n_data, ref f, ref s, ref pr, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                ccdf_lookup = 1.0 - cdf_lookup;

                ompr = 1.0 - pr;

                f_double = ( double ) f;
                s_double = ( double ) s;

                calc.cumnbn ( ref f_double, ref s_double, ref pr, ref ompr, ref cdf_compute, ref ccdf_compute );

                Object[] o = { f, s, pr, ccdf_lookup, ccdf_compute };
                Console.WriteLine("  {0,2}  {1,2}   {2,-8:0.0#####}  {3,-12:0.0#########}  {4,-12:0.0#########}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20 tests CUMNOR, NORMAL_CDF_VALUES.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test20()
        {
            double ccdf_compute = 0.0;
            double ccdf_lookup = 0.0;
            double cdf_compute = 0.0;
            double cdf_lookup = 0.0;
            int n_data = 0;
            double x = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST20");
            Console.WriteLine("  CUMNOR computes the Normal CDF");
            Console.WriteLine("    and the complementary CDF.");
            Console.WriteLine("  NORMAL_CDF_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("    X         CDF           CDF");
            Console.WriteLine("              (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            StatClass calc = new StatClass();

            for ( ; ; )
            {
                calc.normal_cdf_values ( ref n_data, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                calc.cumnor ( ref x, ref cdf_compute, ref ccdf_compute );

                Object[] o = { x, cdf_lookup, cdf_compute };
                Console.WriteLine("  {0,-8:0.0#####}  {1,-12:0.0#########}  {2,-12:0.0#########}", o);
            }

            Console.WriteLine();
            Console.WriteLine("    X         1-CDF         CCDF");
            Console.WriteLine("              (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            for ( ; ; )
            {
                calc.normal_cdf_values ( ref n_data, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                ccdf_lookup = 1.0 - cdf_lookup;

                calc.cumnor ( ref x, ref cdf_compute, ref ccdf_compute );

                Object[] o = { x, ccdf_lookup, ccdf_compute };
                Console.WriteLine("  {0,-8:0.0#####}  {1,-12:0.0#########}  {2,-12:0.0#########}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST21 tests CUMPOI, POISSON_CDF_VALUES.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test21()
        {
            double ccdf_compute = 0.0;
            double ccdf_lookup = 0.0;
            double cdf_compute = 0.0;
            double cdf_lookup = 0.0;
            double lambda = 0.0;
            int n_data = 0;
            int x = 0;
            double x_double = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST21");
            Console.WriteLine("  CUMPOI computes the Poisson CDF");
            Console.WriteLine("    and the complementary CDF.");
            Console.WriteLine("  POISSON_CDF_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("     X    LAMBDA      CDF           CDF");
            Console.WriteLine("                     (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            StatClass calc = new StatClass();

            for ( ; ; )
            {
                calc.poisson_cdf_values ( ref n_data, ref lambda, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                x_double = ( double ) x;

                calc.cumpoi ( ref x_double, ref lambda, ref cdf_compute, ref ccdf_compute );

                Object[] o = { x, lambda, cdf_lookup, cdf_compute };
                Console.WriteLine("    {0,2}    {1,-8:0.0#####}  {2,-12:0.0#########}  {3,-12:0.0#########}", o);
            }

            Console.WriteLine();
            Console.WriteLine("     X    LAMBDA      1-CDF         CCDF");
            Console.WriteLine("                     (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            for ( ; ; )
            {
                calc.poisson_cdf_values ( ref n_data, ref lambda, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                x_double = ( double ) x;
                ccdf_lookup = 1.0 - cdf_lookup;

                calc.cumpoi ( ref x_double, ref lambda, ref cdf_compute, ref ccdf_compute );

                Object[] o = { x, lambda, ccdf_lookup, ccdf_compute };
                Console.WriteLine("    {0,2}    {1,-8:0.0#####}  {2,-12:0.0#########}  {3,-12:0.0#########}", o);

            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST22 tests CUMT, STUDENT_CDF_VALUES.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test22()
        {
            double ccdf_compute = 0.0;
            double ccdf_lookup = 0.0;
            double cdf_compute = 0.0;
            double cdf_lookup = 0.0;
            int df = 0;
            double df_double = 0.0;
            int n_data = 0;
            double x = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST22");
            Console.WriteLine("  CUMT computes the Student T CDF");
            Console.WriteLine("    and the complementary CDF.");
            Console.WriteLine("  STUDENT_CDF_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("    X       DF    CDF           CDF");
            Console.WriteLine("                 (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            StatClass calc = new StatClass();

            for ( ; ; )
            {
                calc.student_cdf_values ( ref n_data, ref df, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }
                df_double = ( double ) df;

                calc.cumt ( ref x, ref df_double, ref cdf_compute, ref ccdf_compute );

                Object[] o = { x, df, cdf_lookup, cdf_compute };
                Console.WriteLine("  {0,-8:0.0#####}  {1,2}    {2,-12:0.0#########} {3,-12:0.0#########}", o);
            }

            Console.WriteLine();
            Console.WriteLine("    X       DF    1-CDF         CCDF");
            Console.WriteLine("                 (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            for (; ; )
            {
                calc.student_cdf_values ( ref n_data, ref df, ref x, ref cdf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                ccdf_lookup = 1.0 - cdf_lookup;

                df_double = ( double ) df;

                calc.cumt ( ref x, ref df_double, ref cdf_compute, ref ccdf_compute );

                Object[] o = { x, df, ccdf_lookup, ccdf_compute };
                Console.WriteLine("  {0,-8:0.0#####}  {1,2}    {2,-12:0.0#########} {3,-12:0.0#########}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST23 tests BETA, GAMMA_X.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test23()
        {
            double a = 0.0;
            double b = 0.0;
            double apb = 0.0;
            double beta1 = 0.0;
            double beta2 = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST23");
            Console.WriteLine("  BETA evaluates the Beta function;");
            Console.WriteLine("  GAMMA_X evaluates the Gamma function.");

            a = 2.2;
            b = 3.7;
            apb = a + b;

            StatClass calc = new StatClass();

            beta1 = calc.beta ( a, b );
            beta2 = calc.gamma_x ( ref a ) * calc.gamma_x ( ref b ) / calc.gamma_x ( ref apb );

            Console.WriteLine();
            Console.WriteLine("  Argument A =                   {0,12:0.0#########}", a);
            Console.WriteLine("  Argument B =                   {0,12:0.0#########}", b);
            Console.WriteLine("  Beta(A,B) =                    {0,12:0.0#########}", beta1);
            Console.WriteLine("  (Expected value = 0.0454 )");
            Console.WriteLine();
            Console.WriteLine("  Gamma(A)*Gamma(B)/Gamma(A+B) = {0,12:0.0#########}", beta2);

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST24 tests ERROR_F, ERROR_FC, ERF_VALUES..
        //
        //  Modified:
        //
        //    17 November 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test24()
        {
            double erf_compute = 0.0;
            double erf_lookup = 0.0;
            double erfc_compute = 0.0;
            double erfc_lookup = 0.0;
            int ind = 0;
            int n_data = 0;
            double x = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST24");
            Console.WriteLine("  ERROR_F computes the error function ERF;");
            Console.WriteLine("  ERROR_FC the complementary error function ERFC.");
            Console.WriteLine("  ERF_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("    X         ERF           ERF");
            Console.WriteLine("              (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            StatClass calc = new StatClass();

            for ( ; ; )
            {
                calc.erf_values ( ref n_data, ref x, ref erf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                erf_compute = calc.error_f ( ref x );

                Object[] o = { x, erf_lookup, erf_compute };
                Console.WriteLine("   {0,-8:0.0#####} {1,-12:0.0#########} {2,-12:0.0#########}", o);
            }

            Console.WriteLine();
            Console.WriteLine("    X         ERFC          ERFC");
            Console.WriteLine("              (Lookup)      (Computed)");
            Console.WriteLine();

            ind = 0;
            n_data = 0;

            for ( ; ; )
            {
                calc.erf_values ( ref n_data, ref x, ref erf_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                erfc_lookup = 1.0 - erf_lookup;
                erfc_compute = calc.error_fc ( ref ind, ref x );

                Object[] o = { x, erfc_lookup, erfc_compute };
                Console.WriteLine("   {0,-8:0.0#####} {1,-12:0.0#########} {2,-12:0.0#########}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST25 tests XGAMM, GAMMA_VALUES.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test25()
        {
            double gamma_compute = 0.0;
            double gamma_lookup = 0.0;
            int n_data = 0;
            double x = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST25");
            Console.WriteLine("  XGAMM computes the Gamma function;");
            Console.WriteLine("  GAMMA_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("    X         GAMMA         GAMMA");
            Console.WriteLine("              (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            StatClass calc = new StatClass();

            for ( ; ; )
            {
                calc.gamma_values ( ref n_data, ref x, ref gamma_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                gamma_compute = calc.gamma_x(ref x);

                Object[] o = { x, gamma_lookup, gamma_compute };
                //Console.WriteLine("   {0,-8:0.0#####} {1,-12:0.0#########} {2,-12:0.0#########}", o);
                Console.WriteLine("   {0,-8:0.0#####} {1,-12:E} {2,-12:E}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST26 tests GAMMA_INC, GAMMA_INC_INV.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test26()
        {
            double a = 3.0;
            int i = 0;
            int ierror = 0;
            int ind = 1;
            double p = 0.0;
            double q = 0.0;
            int test_num = 10;
            double x = 0.0;
            double x0 = 0.0;
            double x2 = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST26");
            Console.WriteLine("  GAMMA_INC evaluates the incomplete Gamma ratio;");
            Console.WriteLine("  GAMMA_INC_INV inverts it.");
            Console.WriteLine();
            Console.WriteLine("  Parameters:");
            Console.WriteLine();
            Console.WriteLine("    A = {0,-12:E}", a );
            Console.WriteLine();
            Console.WriteLine("    X             P            Q          Inverse");
            Console.WriteLine();

            StatClass calc = new StatClass();

            for (i = 0; i <= test_num; i++)
            {
                x = ( double ) i / ( double ) test_num;

                calc.gamma_inc ( ref a, ref x, ref p, ref q, ref ind );

                calc.gamma_inc_inv ( ref a, ref x2, ref x0, ref p, ref q, ref ierror );

                Object[] o = { x, p, q, x2 };
                Console.WriteLine("  {0,4:0.0#}       {1,-12:0.0#########} {2,-12:0.0#########} {3,-12:0.0#########}", o);
            }

            return;
        }
        //****************************************************************************80

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST27 tests PSI, PSI_VALUES.
        //
        //  Modified:
        //
        //    14 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        static void test27()
        {
            double psi_compute = 0.0;
            double psi_lookup = 0.0;
            int n_data = 0;
            double x = 0.0;

            Console.WriteLine();
            Console.WriteLine("TEST27");
            Console.WriteLine("  PSI computes the Psi function;");
            Console.WriteLine("  PSI_VALUES looks up some values.");
            Console.WriteLine();
            Console.WriteLine("    X         PSI           PSI");
            Console.WriteLine("              (Lookup)      (Computed)");
            Console.WriteLine();

            n_data = 0;

            StatClass calc = new StatClass();

            for (; ; )
            {
                calc.psi_values ( ref n_data, ref x, ref psi_lookup );

                if ( n_data == 0 )
                {
                    break;
                }

                psi_compute = calc.psi ( ref x );

                Object[] o = { x, psi_lookup, psi_compute };
                Console.WriteLine("  {0,-8:0.0#####} {1,-12:0.0#########} {2,-12:0.0#########}", o);
            }

            return;
        }
        //****************************************************************************80


        static void Main(string[] args)
        {

            Console.WriteLine();
            Console.WriteLine("DCDFLIB_PRB");
            Console.WriteLine("  C++/CLI safe version");
            Console.WriteLine("  Test the routines in the DCDFLIB library.");

            test005();
            test01();
            test02();
            test03();
            test04();
            test05();
            test06();
            test07();
            test08();
            test09();

            test10();
            test11();
            test12();
            test13();
            test14();
            test15();
            test16();
            test17();
            test18();
            test19();

            test20();
            test21();
            test22();
            test23();
            test24();
            test25();
            test26();
            test27();

            Console.WriteLine();
            Console.WriteLine("DCDFLIB_PRB");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine();

            return;
        }
    }
}
