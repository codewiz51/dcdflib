
using namespace Microsoft::VisualStudio::TestTools::UnitTesting;
using namespace dcdflib;
namespace DCDFlibUTest {
    using namespace System;
    ref class StatClassTest;
    
    
    /// <summary>
///This is a test class for StatClassTest and is intended
///to contain all StatClassTest Unit Tests
///</summary>
	[TestClass]
	public ref class StatClassTest
	{

	private: Microsoft::VisualStudio::TestTools::UnitTesting::TestContext^  testContextInstance;
			 /// <summary>
			 ///Gets or sets the test context which provides
			 ///information about and functionality for the current test run.
			 ///</summary>
	public: property Microsoft::VisualStudio::TestTools::UnitTesting::TestContext^  TestContext
			{
				Microsoft::VisualStudio::TestTools::UnitTesting::TestContext^  get()
				{
					return testContextInstance;
				}
				System::Void set(Microsoft::VisualStudio::TestTools::UnitTesting::TestContext^  value)
				{
					testContextInstance = value;
				}
			}

#pragma region Additional test attributes
			// 
			//You can use the following additional attributes as you write your tests:
			//
			//Use ClassInitialize to run code before running the first test in the class
			//public: [ClassInitialize]
			//static System::Void MyClassInitialize(TestContext^  testContext)
			//{
			//}
			//
			//Use ClassCleanup to run code after all tests in a class have run
			//public: [ClassCleanup]
			//static System::Void MyClassCleanup()
			//{
			//}
			//
			//Use TestInitialize to run code before running each test
			//public: [TestInitialize]
			//System::Void MyTestInitialize()
			//{
			//}
			//
			//Use TestCleanup to run code after each test has run
			//public: [TestCleanup]
			//System::Void MyTestCleanup()
			//{
			//}
			//
#pragma endregion
			/// <summary>
			///A test for stvaln
			///</summary>
	public: [TestMethod]
			void stvalnTest()
			{
				StatClass^  target = gcnew StatClass();
				
				// Test p < 0.5
				double p = .005;
				double pExpected = .005;
				double expected = -2.57582929;
				double actual = target->stvaln(p);
				double diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(pExpected, p);
				Assert::IsTrue(diff < 1.0e-06);

				// Test p == 0.5
				p = .5;
				pExpected = .5;
				expected = -0.00000001;
				actual = target->stvaln(p);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(pExpected, p);
				Assert::IsTrue(diff < 1.0e-06);

				// Test p > 0.5
				p = .9;
				pExpected = .9;
				expected = 1.28155156;
				actual = target->stvaln(p);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(pExpected, p);
				Assert::IsTrue(diff < 1.0e-06);

				// Test ArgumentException
				bool exceptionThrown = false;
				try
				{
					// Setting p > 1.0 will throw an ArgumentException
					p = 1.5;
					actual = target->stvaln(p);
				}
				catch (ArgumentException^)
				{
					exceptionThrown = true;
				}
				Assert::IsTrue(exceptionThrown);

			}
			/// <summary>
			///A test for rlog1
			///</summary>
	public: [TestMethod]
			void rlog1Test()
			{
				StatClass^  target = gcnew StatClass();
				double x = -0.98;
				double xExpected = -0.98;
				double expected = 2.93202300;
				double actual = target->rlog1(x);
				double diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 1.0e-06);
				
				x = -0.27;
				xExpected = -0.27;
				expected = 0.04471074;
				actual = target->rlog1(x);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 1.0e-06);

				x = -0.09;
				xExpected = -0.09;
				expected = 0.00431068;
				actual = target->rlog1(x);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 1.0e-06);

				x = 0.09;
				xExpected = 0.09;
				expected = 0.00382230;
				actual = target->rlog1(x);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 9.0e-09);

				x = 0.40;
				xExpected = 0.40;
				expected = 0.06352776;
				actual = target->rlog1(x);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 9.0e-09);

				x = 2.10;
				xExpected = 2.10;
				expected = 0.96859789;
				actual = target->rlog1(x);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 1.0e-06);

				// Test ArgumentException
				bool exceptionThrown = false;
				try
				{
					// Setting p <= -1.0 will throw an ArgumentException
					x = -1.01;
					actual = target->rlog1(x);
				}
				catch (ArgumentException^)
				{
					exceptionThrown = true;
				}
				Assert::IsTrue(exceptionThrown);

			}
			/// <summary>
			///A test for rlog
			///</summary>
	public: [TestMethod]
			void rlogTest()
			{
				StatClass^  target = (gcnew StatClass()); // TODO: Initialize to an appropriate value
				double x = 0.05;
				double xExpected = 0.05;
				double expected = 2.04573227;
				double actual = target->rlog(x);
				double diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 1.0e-06);
				
				x = 0.62;
				xExpected = 0.62;
				expected = 0.09803580;
				actual = target->rlog(x);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 1.0e-06);

				x = 1.05;
				xExpected = 1.05;
				expected = 0.00120984;
				actual = target->rlog(x);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 1.0e-06);

				x = 1.19;
				xExpected = 1.19;
				expected = 0.01604669;
				actual = target->rlog(x);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 9.0e-09);

				x = 1.59;
				xExpected = 1.59;
				expected = 0.12626598;
				actual = target->rlog(x);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 1.0e-06);

				x = 3.00;
				xExpected = 3.00;
				expected = 0.90138771;
				actual = target->rlog(x);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 1.0e-06);

				// Test ArgumentException
				bool exceptionThrown = false;
				try
				{
					// Setting p <= -1.0 will throw an ArgumentException
					x = -1.01;
					actual = target->rlog(x);
				}
				catch (ArgumentException^)
				{
					exceptionThrown = true;
				}
				Assert::IsTrue(exceptionThrown);
			}
			/// <summary>
			///A test for rexp
			///</summary>
	public: [TestMethod]
			void rexpTest()
			{
				StatClass^  target = gcnew StatClass();
				double x = -1.1;
				double xExpected = -1.1;
				double expected = -0.66712892;
				double actual = target->rexp(x);
				double diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 1.0e-06);

				x = 0.001;
				xExpected = 0.001;
				expected = 0.00100050;
				actual = target->rexp(x);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 1.0e-06);

				x = 0.50;
				xExpected = 0.50;
				expected = 0.64872127;
				actual = target->rexp(x);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 1.0e-06);

				x = 2.00;
				xExpected = 2.00;
				expected = 6.38905610;
				actual = target->rexp(x);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(xExpected, x);
				Assert::IsTrue(diff < 1.0e-06);
			}
			/// <summary>
			///A test for gam1
			///</summary>
	public: [TestMethod]
			void gam1Test()
			{
				StatClass^  target = (gcnew StatClass());
				double a = -0.4;
				double aExpected = -0.4;
				double expected = -0.32849503;
				double actual = target->gam1(a);
				double diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(aExpected, a);
				Assert::IsTrue(diff < 1.0e-06);

				a = 0.0;
				aExpected = 0.0;
				expected = 0.00000000;
				actual = target->gam1(a);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(aExpected, a);
				Assert::IsTrue(diff < 1.0e-06);

				a = 1.e-3;
				aExpected = 1.e-3;
				expected = 0.00057656;
				actual = target->gam1(a);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(aExpected, a);
				Assert::IsTrue(diff < 1.0e-06);

				a = 1.24;
				aExpected = 1.24;
				expected = -0.11234681;
				actual = target->gam1(a);
				diff = Math::Abs(Math::Abs(expected) - Math::Abs(actual));
				Assert::AreEqual(aExpected, a);
				Assert::IsTrue(diff < 1.0e-06);

				// Test ArgumentException
				bool exceptionThrown = false;
				try
				{
					// Setting p <= -1.0 will throw an ArgumentException
					a = -2.01;
					actual = target->gam1(a);
				}
				catch (ArgumentException^)
				{
					exceptionThrown = true;
				}
				Assert::IsTrue(exceptionThrown);
			}
			/// <summary>
			///A test for test005
			///</summary>
	public: [TestMethod]
			void test005()
			{
				StatClass^  target = gcnew StatClass();
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

				n_data = 0;

				for (; ; )
				{
					target->beta_inc_values(n_data, a, b, x, cdf_lookup);

					if (n_data == 0)
					{
						break;
					}

					y = 1.0 - x;

					target->beta_inc(a, b, x, y, cdf_compute, ccdf_compute, ierror);
					double diff = Math::Abs(Math::Abs(cdf_lookup) - Math::Abs(cdf_compute));
					Assert::IsTrue(diff < 1.0e-04);

				}

				n_data = 0;

				for (; ; )
				{
					target->beta_inc_values(n_data, a, b, x, cdf_lookup);

					if (n_data == 0)
					{
						break;
					}

					ccdf_lookup = 1.0 - cdf_lookup;

					y = 1.0 - x;

					target->beta_inc(a, b, x, y, cdf_compute, ccdf_compute, ierror);

					double diff = Math::Abs(Math::Abs(ccdf_lookup) - Math::Abs(ccdf_compute));
					Assert::IsTrue(diff < 1.0e-04);
				}


			}
			/// <summary>
			///A test for test01
			///</summary>
	public: [TestMethod]
			void test01()
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
				
				StatClass^  target = gcnew StatClass();

				which = 1;
				p = -1.0;
				q = -1.0;
				x = 0.25;
				y = 1.0 - x;
				a = 2.0;
				b = 3.0;
				target->cdfbet(which, p, q, x, y, a, b, status, bound);
				double expectedP = 0.26171875;
				double expectedQ = 0.73828125;
				double diffP = Math::Abs(Math::Abs(expectedP) - Math::Abs(p));
				double diffQ = Math::Abs(Math::Abs(expectedQ) - Math::Abs(q));
				Assert::IsTrue(diffP < 1.0e-05);
				Assert::IsTrue(diffQ < 1.0e-05);

				which = 2;
				p = 0.261719;
				q = 1.0 - p;
				x = -1.0;
				y = -1.0;
				a = 2.0;
				b = 3.0;
				double expectedX = 0.25;
				double expectedY = 0.75;
				target->cdfbet(which, p, q, x, y, a, b, status, bound);
				double diffX = Math::Abs(Math::Abs(expectedX) - Math::Abs(x));
				double diffY = Math::Abs(Math::Abs(expectedY) - Math::Abs(y));
				Assert::IsTrue(diffX < 1.0e-05);
				Assert::IsTrue(diffY < 1.0e-05);

				which = 3;
				p = 0.261719;
				q = 1.0 - p;
				x = 0.25;
				y = 1.0 - x;
				a = -1.0;
				b = 3.0;
				target->cdfbet(which, p, q, x, y, a, b, status, bound);
				double expectedA = 2.0;
				double diffA = Math::Abs(Math::Abs(expectedA) - Math::Abs(a));
				Assert::IsTrue(diffA < 1.0e-05);

				which = 4;
				p = 0.261719;
				q = 1.0 - p;
				x = 0.25;
				y = 1.0 - x;
				a = 2.0;
				b = -1.0;
				target->cdfbet(which, p, q, x, y, a, b, status, bound);
				double expectedB = 3.0;
				double diffB = Math::Abs(Math::Abs(expectedB) - Math::Abs(b));
				Assert::IsTrue(diffB < 1.0e-05);

			}
			/// <summary>
			///A test for test02
			///</summary>
	public: [TestMethod]
			void test02()
			{
				StatClass^  target = (gcnew StatClass());
				double bound = 0;
				double ompr = 0;
				double p = 0;
				double pr = 0;
				double q = 0;
				double s = 0;
				int status = 0;
				int which = 0;
				double xn = 0;

                which = 1;
				p = -1.0;
                q = -1.0;
                s = 5.0;
                xn = 8.0;
                pr = 0.875;
                ompr = 1.0 - pr;
				target->cdfbin(which, p, q, s, xn, pr, ompr, status, bound);
				double expectedP = 0.067347;
				double expectedQ = 0.932653;
				double diffP = Math::Abs(Math::Abs(expectedP) - Math::Abs(p));
				double diffQ = Math::Abs(Math::Abs(expectedQ) - Math::Abs(q));
				Assert::IsTrue(diffP < 1.0e-05);
				Assert::IsTrue(diffQ < 1.0e-05);

                which = 2;
                p = 0.067347;
                q = 1.0 - p;
                s = -1.0;
                xn = 8.0;
                pr = 0.875;
                ompr = 1.0 - pr;
				target->cdfbin(which, p, q, s, xn, pr, ompr, status, bound);
				double expectedS = 5.0;
				double diffS = Math::Abs(Math::Abs(expectedS) - Math::Abs(s));
				Assert::IsTrue(diffS < 1.0e-05);

				which = 3;
                p = 0.067347;
                q = 1.0 - p;
                s = 5.0;
                xn = -1.0;
                pr = 0.875;
                ompr = 1.0 - pr;
				target->cdfbin(which, p, q, s, xn, pr, ompr, status, bound);
				double expectedXN = 8.0;
				double diffXN = Math::Abs(Math::Abs(expectedXN) - Math::Abs(xn));
				Assert::IsTrue(diffXN < 1.0e-05);

				which = 4;
                p = 0.067347;
                q = 1.0 - p;
                s = 5.0;
                xn = 8.0;
                pr = -1.0;
                ompr = -1.0;
				target->cdfbin(which, p, q, s, xn, pr, ompr, status, bound);
				double expectedPR = 0.875;
				double expectedOMPR = 0.125;
				double diffPR = Math::Abs(Math::Abs(expectedPR) - Math::Abs(pr));
				double diffOMPR = Math::Abs(Math::Abs(expectedOMPR) - Math::Abs(ompr));
				Assert::IsTrue(diffPR < 1.0e-05);
				Assert::IsTrue(diffOMPR < 1.0e-05);

			}
			/// <summary>
			///A test for test03
			///</summary>
	public: [TestMethod]
			void test03()
			{
				StatClass^  target = (gcnew StatClass()); // TODO: Initialize to an appropriate value
				double bound = 0;
				double df = 0;
				double p = 0;
				double q = 0;
				int status = 0;
				int which = 0;
				double x = 0;

                which = 1;
				p = -1.0;
                q = -1.0;
                x = 5.0;
                df = 8.0;
				target->cdfchi(which, p, q, x, df, status, bound);
				double expectedP = 0.242424;
				double expectedQ = 0.757576;
				double diffP = Math::Abs(Math::Abs(expectedP) - Math::Abs(p));
				double diffQ = Math::Abs(Math::Abs(expectedQ) - Math::Abs(q));
				Assert::IsTrue(diffP < 1.0e-05);
				Assert::IsTrue(diffQ < 1.0e-05);

				which = 2;
				p = 0.242424;
				q = 1.0 - p;
				x = -1.0;
				df = 8.0;
				target->cdfchi(which, p, q, x, df, status, bound);
				double expectedX = 5.0;
				double diffX = Math::Abs(Math::Abs(expectedX) - Math::Abs(x));
				Assert::IsTrue(diffX < 1.0e-05);

                which = 3;
                p = 0.242424;
                q = 1.0 - p;
                x = 5.0;
                df = -1.0;
				target->cdfchi(which, p, q, x, df, status, bound);
				double expectedDF = 8.0;
				double diffDF = Math::Abs(Math::Abs(expectedDF) - Math::Abs(df));
				Assert::IsTrue(diffDF < 1.0e-05);

			}
			
			/// <summary>
			///A test for test04
			///</summary>
	public: [TestMethod]
			void test04()
			{
				StatClass^  target = gcnew StatClass();
				double bound = 0;
				double df = 0;
				double p = 0;
				double pnonc = 0;
				double q = 0;
				int status = 0;
				int which = 0;
				double x = 0;

                which = 1;
				p = -1.0;
                q = -1.0;
                x = 5.0;
                df = 8.0;
                pnonc = 0.5;
                target->cdfchn(which, p, q, x, df, pnonc, status, bound);
				double expectedP = 0.210747;
				double expectedQ = 0.789253;
				double diffP = Math::Abs(Math::Abs(expectedP) - Math::Abs(p));
				double diffQ = Math::Abs(Math::Abs(expectedQ) - Math::Abs(q));
				Assert::IsTrue(diffP < 1.0e-05);
				Assert::IsTrue(diffQ < 1.0e-05);

				which = 2;
				p = 0.211040;
                q = 1.0 - p;
                x = -1.0;
                df = 8.0;
                pnonc = 0.5;
                target->cdfchn(which, p, q, x, df, pnonc, status, bound);
				double expectedX = 5.0;
				double diffX = Math::Abs(Math::Abs(expectedX) - Math::Abs(x));
				Assert::IsTrue(diffX < 1.0e-02);

				which = 3;
                p = 0.211040;
                q = 1.0 - p;
                x = 5.0;
                df = -1.0;
                pnonc = 0.5;
                target->cdfchn(which, p, q, x, df, pnonc, status, bound);
				double expectedDF = 8.0;
				double diffDF = Math::Abs(Math::Abs(expectedDF) - Math::Abs(df));
				Assert::IsTrue(diffDF < 1.0e-02);

				which = 4;
                p = 0.211040;
                q = 1.0 - p;
                x = 5.0;
                df = 8.0;
                pnonc = -1.0;
                target->cdfchn(which, p, q, x, df, pnonc, status, bound);
				double expectedpnonc = 0.5;
				double diffpnonc = Math::Abs(Math::Abs(expectedpnonc) - Math::Abs(pnonc));
				Assert::IsTrue(diffpnonc < 1.0e-02);

			}
			/// <summary>
			///A test for test05
			///</summary>
	public: [TestMethod]
			void test05()
			{
				StatClass^  target = (gcnew StatClass());
				double bound = 0;
				double dfd = 0;
				double dfn = 0;
				double f = 0;
				double p = 0;
				double q = 0;
				int status = 0;
				int which = 0;

				which = 1;
				p = -1.0;
				q = -1.0;
				f = 5.0;
				dfn = 8.0;
				dfd = 3.0;
				target->cdff(which, p, q, f, dfn, dfd, status, bound);
				double expectedP = 0.89351;
				double expectedQ = 0.10649;
				double diffP = Math::Abs(Math::Abs(expectedP) - Math::Abs(p));
				double diffQ = Math::Abs(Math::Abs(expectedQ) - Math::Abs(q));
				Assert::IsTrue(diffP < 1.0e-05);
				Assert::IsTrue(diffQ < 1.0e-05);

				which = 2;
				p = 0.893510;
                q = 1.0 - p;
                f = -1.0;
                dfn = 8.0;
                dfd = 3.0;
				target->cdff(which, p, q, f, dfn, dfd, status, bound);
				double expectedF = 5.0;
				double diffF = Math::Abs(Math::Abs(expectedF) - Math::Abs(f));
				Assert::IsTrue(diffF < 1.0e-05);

				which = 3;
				p = 0.893510;
                q = 1.0 - p;
                f = 5.0;
                dfn = -1.0;
                dfd = 3.0;
				target->cdff(which, p, q, f, dfn, dfd, status, bound);
				double expecteddfn = 8.0;
				double diffdfn = Math::Abs(Math::Abs(expecteddfn) - Math::Abs(dfn));
				Assert::IsTrue(diffdfn < 1.0e-02);

				which = 4;
				p = 0.893510;
                q = 1.0 - p;
                f = 5.0;
                dfn = 8.0;
                dfd = -1.0;
				target->cdff(which, p, q, f, dfn, dfd, status, bound);
				double expecteddfd = 3.0;
				double diffdfd = Math::Abs(Math::Abs(expecteddfd) - Math::Abs(dfd));
				Assert::IsTrue(diffdfd < 1.0e-05);

			}
			/// <summary>
			///A test for test06
			///</summary>
	public: [TestMethod]
			void test06()
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

				StatClass^  target = (gcnew StatClass());

                which = 1;
				p = -1.0;
                q = -1.0;
                f = 5.0;
                dfn = 8.0;
                dfd = 3.0;
                pnonc = 17.648016;
                target->cdffnc(which, p, q, f, dfn, dfd, pnonc, status, bound);
				double expectedP = 0.6;
				double expectedQ = 0.4;
				double diffP = Math::Abs(Math::Abs(expectedP) - Math::Abs(p));
				double diffQ = Math::Abs(Math::Abs(expectedQ) - Math::Abs(q));
				Assert::IsTrue(diffP < 1.0e-05);
				Assert::IsTrue(diffQ < 1.0e-05);

                which = 2;
				p = 0.60;
                q = 1.0 - p;
                f = -1.0;
                dfn = 8.0;
                dfd = 3.0;
                pnonc = 17.648016;
                target->cdffnc(which, p, q, f, dfn, dfd, pnonc, status, bound);
				double expectedF = 5.0;
				double diffF = Math::Abs(Math::Abs(expectedF) - Math::Abs(f));
				Assert::IsTrue(diffF < 1.0e-05);

                which = 3;
				p = 0.60;
                q = 1.0 - p;
                f = 5.0;
                dfn = -1.0;
                dfd = 3.0;
                pnonc = 17.648016;
                target->cdffnc(which, p, q, f, dfn, dfd, pnonc, status, bound);
				double expecteddfn = 8.0;
				double diffdfn = Math::Abs(Math::Abs(expecteddfn) - Math::Abs(dfn));
				Assert::IsTrue(diffdfn < 1.0e-02);

                which = 4;
				p = 0.60;
                q = 1.0 - p;
                f = 5.0;
                dfn = 8.0;
                dfd = -1.0;
                pnonc = 17.648016;
                target->cdffnc(which, p, q, f, dfn, dfd, pnonc, status, bound);
				double expecteddfd = 3.0;
				double diffdfd = Math::Abs(Math::Abs(expecteddfd) - Math::Abs(dfd));
				Assert::IsTrue(diffdfd < 1.0e-05);

                which = 5;
				p = 0.60;
                q = 1.0 - p;
                f = 5.0;
                dfn = 8.0;
                dfd = 3.0;
                pnonc = -1.0;
                target->cdffnc(which, p, q, f, dfn, dfd, pnonc, status, bound);
				double expectedpnonc = 17.64801798;
				double diffdpnonc = Math::Abs(Math::Abs(expectedpnonc) - Math::Abs(pnonc));
				Assert::IsTrue(diffdpnonc < 1.0e-05);
			}
			/// <summary>
			///A test for test07
			///</summary>
	public: [TestMethod]
			void test07()
			{
				double bound = 0;
				double p = 0;
				double q = 0;
				double scale = 0;
				double shape = 0;
				int status = 0;
				int which = 0;
				double x = 0;

				StatClass^  target = gcnew StatClass();
                which = 1;
				p = -1.0;
                q = -1.0;
                x = 5.0;
                shape = 8.0;
                scale = 3.0;
                target->cdfgam(which, p, q, x, shape, scale, status, bound);
				double expectedP = 0.981998;
				double expectedQ = 0.0180022;
				double diffP = Math::Abs(Math::Abs(expectedP) - Math::Abs(p));
				double diffQ = Math::Abs(Math::Abs(expectedQ) - Math::Abs(q));
				Assert::IsTrue(diffP < 2.0e-07);
				Assert::IsTrue(diffQ < 2.0e-07);

				which = 2;
                p = 0.981998;
                q = 1.0 - p;
                x = -1.0;
                shape = 8.0;
                scale = 3.0;
                target->cdfgam(which, p, q, x, shape, scale, status, bound);
				double expectedX = 5.0;
				double diffX = Math::Abs(Math::Abs(expectedX) - Math::Abs(x));
				Assert::IsTrue(diffX < 7.0e-06);

				which = 3;
				p = 0.981998;
                q = 1.0 - p;
                x = 5.0;
                shape = -1.0;
                scale = 3.0;
                target->cdfgam(which, p, q, x, shape, scale, status, bound);
				double expectedshape = 8.0;
				double diffshape = Math::Abs(Math::Abs(expectedshape) - Math::Abs(shape));
				Assert::IsTrue(diffshape < 2.0e-05);

				which = 4;
				p = 0.981998;
                q = 1.0 - p;
                x = 5.0;
                shape = 8.0;
                scale = -1.0;
                target->cdfgam(which, p, q, x, shape, scale, status, bound);
				double expectedscale = 3.0;
				double diffscale = Math::Abs(Math::Abs(expectedscale) - Math::Abs(scale));
				Assert::IsTrue(diffscale < 4.0e-06);

			}
			/// <summary>
			///A test for test08
			///</summary>
	public: [TestMethod]
			void test08()
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
				StatClass^  target = (gcnew StatClass());

				which = 1;
				p = -1.0;
                q = -1.0;
                f = 3.0;
                s = 5.0;
                pr = 0.875;
                ompr = 1.0 - pr;
                target->cdfnbn(which, p, q, f, s, pr, ompr, status, bound);
				double expectedP = 0.988752;
				double expectedQ = 0.0112478;
				double diffP = Math::Abs(Math::Abs(expectedP) - Math::Abs(p));
				double diffQ = Math::Abs(Math::Abs(expectedQ) - Math::Abs(q));
				Assert::IsTrue(diffP < 2.0e-07);
				Assert::IsTrue(diffQ < 2.0e-07);

				which = 2;
                p = 0.988752;
                q = 1.0 - p;
                f = -1.0;
                s = 5.0;
                pr = 0.875;
                ompr = 1.0 - pr;
                target->cdfnbn(which, p, q, f, s, pr, ompr, status, bound);
				double expectedF = 3.0;
				double diffF = Math::Abs(Math::Abs(expectedF) - Math::Abs(f));
				Assert::IsTrue(diffF < 2.0e-05);

				which = 3;
				p = 0.988752;
                q = 1.0 - p;
                f = 3.0;
                s = -1.0;
                pr = 0.875;
                ompr = 1.0 - pr;
                target->cdfnbn(which, p, q, f, s, pr, ompr, status, bound);
				double expectedS = 5.0;
				double diffS = Math::Abs(Math::Abs(expectedS) - Math::Abs(s));
				Assert::IsTrue(diffS < 5.0e-05);

				which = 4;
                p = 0.988752;
                q = 1.0 - p;
                f = 3.0;
                s = 5.0;
                pr = -1.0;
                ompr = -1.0;
                target->cdfnbn(which, p, q, f, s, pr, ompr, status, bound);
				double expectedPR = 0.874999;
				double expectedOMPR = 0.125001;
				double diffPR = Math::Abs(Math::Abs(expectedPR) - Math::Abs(pr));
				double diffOMPR = Math::Abs(Math::Abs(expectedOMPR) - Math::Abs(ompr));
				Assert::IsTrue(diffPR < 5.0e-07);
				Assert::IsTrue(diffOMPR < 5.0e-07);


			}
			/// <summary>
			///A test for test09
			///</summary>
	public: [TestMethod]
			void test09()
			{
				double bound = 0.0;
				double mean = 0.0;
				double p = 0.0;
				double q = 0.0;
				double sd = 0.0;
				int status = 0;
				int which = 0;
				double x = 0.0;

				StatClass^  target = (gcnew StatClass()); // TODO: Initialize to an appropriate value
            
                which = 1;
				p = -1.0;
                q = -1.0;
                x = 3.0;
                mean = 5.0;
                sd = 0.875;
                target->cdfnor (which, p, q, x, mean, sd, status, bound );
				double expectedP = 0.011135;
				double expectedQ = 0.988865;
				double diffP = Math::Abs(Math::Abs(expectedP) - Math::Abs(p));
				double diffQ = Math::Abs(Math::Abs(expectedQ) - Math::Abs(q));
				Assert::IsTrue(diffP < 9.0e-07);
				Assert::IsTrue(diffQ < 9.0e-07);

                which = 2;
				p = 0.011135;
                q = 1.0 - p;
                x = -1.0;
                mean = 5.0;
                sd = 0.875;
                target->cdfnor (which, p, q, x, mean, sd, status, bound );
				double expectedX = 3.0;
				double diffX = Math::Abs(Math::Abs(expectedX) - Math::Abs(x));
				Assert::IsTrue(diffX < 2.0e-05);

                which = 3;
				p = 0.011135;
                q = 1.0 - p;
                x = 3.0;
                mean = -1.0;
                sd = 0.875;
                target->cdfnor (which, p, q, x, mean, sd, status, bound );
				double expectedMean = 5.0;
				double diffMean = Math::Abs(Math::Abs(expectedMean) - Math::Abs(mean));
				Assert::IsTrue(diffMean < 2.0e-05);

                which = 4;
				p = 0.011135;
                q = 1.0 - p;
                x = 3.0;
                mean = 5.0;
                sd = -1.0;
                target->cdfnor (which, p, q, x, mean, sd, status, bound );
				double expectedSD = 0.875;
				double diffSD = Math::Abs(Math::Abs(expectedSD) - Math::Abs(sd));
				Assert::IsTrue(diffSD < 2.0e-05);

			}
			/// <summary>
			///A test for test10
			///</summary>
	public: [TestMethod]
			void test10()
			{
				double bound = 0.0;
				double p = 0.0;
				double q = 0.0;
				double s = 0.0;
				int status = 0;
				int which = 0;
				double xlam = 0.0;

				StatClass^ target = gcnew StatClass();

                which = 1;
                p = -1.0;
                q = -1.0;
                s = 3.0;
                xlam = 5.0;
                target->cdfpoi ( which, p, q, s, xlam, status, bound );
				double expectedP = 0.265026;
				double expectedQ = 0.734974;
				double diffP = Math::Abs(Math::Abs(expectedP) - Math::Abs(p));
				double diffQ = Math::Abs(Math::Abs(expectedQ) - Math::Abs(q));
				Assert::IsTrue(diffP < 9.0e-07);
				Assert::IsTrue(diffQ < 9.0e-07);

				which = 2;
                p = 0.265026;
                q = 1.0 - p;
                s = -1.0;
                xlam = 5.0;
                target->cdfpoi ( which, p, q, s, xlam, status, bound );
				double expectedS = 3.0;
				double diffS = Math::Abs(Math::Abs(expectedS) - Math::Abs(s));
				Assert::IsTrue(diffS < 9.0e-07);

				which = 3;
                p = 0.265026;
                q = 1.0 - p;
                s = 3.0;
                xlam = -1.0;
                target->cdfpoi ( which, p, q, s, xlam, status, bound );
				double expectedxlam = 5.0;
				double diffxlam = Math::Abs(Math::Abs(expectedxlam) - Math::Abs(xlam));
				Assert::IsTrue(diffxlam < 9.0e-07);

			}
			/// <summary>
			///A test for test11
			///</summary>
	public: [TestMethod]
			void test11()
			{
				double bound = 0.0;
				double df = 0.0;
				double p = 0.0;
				double q = 0.0;
				int status = 0;
				double t = 0.0;
				int which = 0;
				StatClass^ target = gcnew StatClass();
                
				which = 1;
                p = -1.0;
                q = -1.0;
                t = 3.0;
                df = 5.0;
                target->cdft ( which, p, q, t, df, status, bound );
				double expectedP = 0.984950;
				double expectedQ = 1.0 - expectedP;
				double diffP = Math::Abs(Math::Abs(expectedP) - Math::Abs(p));
				double diffQ = Math::Abs(Math::Abs(expectedQ) - Math::Abs(q));
				Assert::IsTrue(diffP < 9.0e-07);
				Assert::IsTrue(diffQ < 9.0e-07);

				which = 2;
                p = 0.984950;
                q = 1.0 - p;
                t = -1.0;
                df = 5.0;
                target->cdft ( which, p, q, t, df, status, bound );
				double expectedT = 3.0;
				double diffT = Math::Abs(Math::Abs(expectedT) - Math::Abs(t));
				Assert::IsTrue(diffT < 3.0e-05);
				
				which = 3;
                p = 0.984950;
                q = 1.0 - p;
                t = 3.0;
                df = -1.0;
                target->cdft ( which, p, q, t, df, status, bound );
				double expectedDF = 5.0;
				double diffDF = Math::Abs(Math::Abs(expectedDF) - Math::Abs(df));
				Assert::IsTrue(diffDF < 1.0e-04);

			}
			/// <summary>
			///A test for test12
			///</summary>
	public: [TestMethod]
			void test12()
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

				StatClass^  target = gcnew StatClass();
				n_data = 0;

				for ( ; ; )
				{
					target->beta_inc_values ( n_data, a, b, x, cdf_lookup );
					if ( n_data == 0 )
						break;
					y = 1.0 - x;
					target->cumbet ( x, y, a, b, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(cdf_compute) - Math::Abs(cdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-06);
				}

				n_data = 0;

				for ( ; ; )
				{
					target->beta_inc_values ( n_data, a, b, x, cdf_lookup );
					if ( n_data == 0 )
						break;
					ccdf_lookup = 1.0 - cdf_lookup;
					y = 1.0 - x;
					target->cumbet ( x, y, a, b, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(ccdf_compute) - Math::Abs(ccdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-06);
				}
			}
			/// <summary>
			///A test for test13
			///</summary>
	public: [TestMethod]
			void test13()
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

				StatClass^  target = gcnew StatClass();

				n_data = 0;

				for ( ; ; )
				{
					target->binomial_cdf_values(n_data, x, pr, s, cdf_lookup);
					if ( n_data == 0 )
						break;
					ompr = 1.0 - pr;
					s_double = ( double ) s;
					x_double = ( double ) x;
					target->cumbin (s_double, x_double, pr, ompr, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(cdf_compute) - Math::Abs(cdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-05);
				}

				n_data = 0;

				for ( ; ; )
				{
					target->binomial_cdf_values ( n_data, x, pr, s, cdf_lookup );
					if ( n_data == 0 )
						break;
					ccdf_lookup = 1.0 - cdf_lookup;
					ompr = 1.0 - pr;
					s_double = ( double ) s;
					x_double = ( double ) x;
					target->cumbin (s_double, x_double, pr, ompr, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(ccdf_compute) - Math::Abs(ccdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-05);
				}
			}
			/// <summary>
			///A test for test14
			///</summary>
	public: [TestMethod]
			void test14()
			{
				double ccdf_compute = 0.0;
				double ccdf_lookup = 0.0;
				double cdf_compute = 0.0;
				double cdf_lookup = 0.0;
				int df = 0;
				double df_double = 0.0;
				int n_data = 0;
				double x = 0.0;

				StatClass^ target = gcnew StatClass();

				n_data = 0;

				for ( ; ; )
				{
					target->chi_square_cdf_values (n_data, df, x, cdf_lookup );
					if ( n_data == 0 )
						break;
					df_double = ( double ) df;
					target->cumchi (x, df_double, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(cdf_compute) - Math::Abs(cdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-06);
				}

				n_data = 0;

				for ( ; ; )
				{
					target->chi_square_cdf_values (n_data, df, x, cdf_lookup );
					if ( n_data == 0 )
						break;
					ccdf_lookup = 1.0 - cdf_lookup;
					df_double = ( double ) df;
					target->cumchi (x, df_double, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(ccdf_compute) - Math::Abs(ccdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-06);
				}
			}
			/// <summary>
			///A test for test15
			///</summary>
	public: [TestMethod]
			void test15()
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

				StatClass^  target = gcnew StatClass();

				n_data = 0;

				for ( ; ; )
				{
					target->chi_noncentral_cdf_values (n_data, x, lambda, df, cdf_lookup );
					if ( n_data == 0 )
						break;
					df_double = ( double ) df;
					target->cumchn(x, df_double, lambda, cdf_compute, ccdf_compute);
					double diffDF = Math::Abs(Math::Abs(cdf_compute) - Math::Abs(cdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-02);
				}

				n_data = 0;

				for ( ; ; )
				{
					target->chi_noncentral_cdf_values (n_data, x, lambda, df, cdf_lookup );
					if ( n_data == 0 )
						break;
					ccdf_lookup = 1.0 - cdf_lookup;
					df_double = ( double ) df;
					target->cumchn (x, df_double, lambda, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(ccdf_compute) - Math::Abs(ccdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-02);
				}
			}
			/// <summary>
			///A test for test16
			///</summary>
	public: [TestMethod]
			void test16()
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

				StatClass^  target = gcnew StatClass();

				n_data = 0;

				for ( ; ; )
				{
					target->f_cdf_values (n_data, dfn, dfd, x, cdf_lookup );
					if ( n_data == 0 )
						break;
					dfn_double = ( double ) dfn;
					dfd_double = ( double ) dfd;
					target->cumf (x, dfn_double, dfd_double, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(cdf_compute) - Math::Abs(cdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-06);

				}

				n_data = 0;

				for ( ; ; )
				{
					target->f_cdf_values(n_data, dfn, dfd, x, cdf_lookup );
					if ( n_data == 0 )
						break;
					ccdf_lookup = 1.0 - cdf_lookup;
					dfn_double = ( double ) dfn;
					dfd_double = ( double ) dfd;
					target->cumf(x, dfn_double, dfd_double, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(ccdf_compute) - Math::Abs(ccdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-07);
				}
			}
			/// <summary>
			///A test for test17
			///</summary>
	public: [TestMethod]
			void test17()
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

				StatClass^  target = gcnew StatClass();

				n_data = 0;
				for ( ; ; )
				{
					target->f_noncentral_cdf_values (n_data, dfn, dfd, lambda, x, cdf_lookup );
					if ( n_data == 0 )
						break;
					dfn_double = ( double ) dfn;
					dfd_double = ( double ) dfd;
					target->cumfnc (x, dfn_double, dfd_double, lambda, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(cdf_compute) - Math::Abs(cdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-06);
				}

				n_data = 0;
				for ( ; ; )
				{
					target->f_noncentral_cdf_values (n_data, dfn, dfd, lambda, x, cdf_lookup );
					if ( n_data == 0 )
						break;
					ccdf_lookup = 1.0 - cdf_lookup;
					dfn_double = ( double ) dfn;
					dfd_double = ( double ) dfd;
					target->cumfnc (x, dfn_double, dfd_double, lambda, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(ccdf_compute) - Math::Abs(ccdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-06);
				}
			}
			/// <summary>
			///A test for test18
			///</summary>
	public: [TestMethod]
			void test18()
			{
				double a = 0.0;
				double ccdf_compute = 0.0;
				double ccdf_lookup = 0.0;
				double cdf_compute = 0.0;
				double cdf_lookup = 0.0;
				int n_data = 0;
				double x = 0.0;

				StatClass^  target = gcnew StatClass();

				n_data = 0;
				for ( ; ; )
				{
					target->gamma_inc_values(n_data, a, x, cdf_lookup );
					if ( n_data == 0 )
						break;
					target->cumgam(x, a, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(cdf_compute) - Math::Abs(cdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-06);
				}

				n_data = 0;
				for ( ; ; )
				{
					target->gamma_inc_values(n_data, a, x, cdf_lookup );
					if ( n_data == 0 )
						break;
					ccdf_lookup = 1.0 - cdf_lookup;
					target->cumgam (x, a, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(ccdf_compute) - Math::Abs(ccdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-06);
				}
		
				// Special values
				n_data = -2;
				target->gamma_inc_values(n_data, a, x, cdf_lookup );
				Assert::AreEqual(1, n_data);


			}
			/// <summary>
			///A test for test19
			///</summary>
	public: [TestMethod]
			void test19()
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

				StatClass^  target = gcnew StatClass();

				n_data = 0;
				for ( ; ; )
				{
					target->negative_binomial_cdf_values (n_data, f, s, pr, cdf_lookup );
					if ( n_data == 0 )
						break;
					ompr = 1.0 - pr;
					f_double = ( double ) f;
					s_double = ( double ) s;
					target->cumnbn (f_double, s_double, pr, ompr, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(cdf_compute) - Math::Abs(cdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-05);
				}

				n_data = 0;
				for ( ; ; )
				{
					target->negative_binomial_cdf_values (n_data, f, s, pr, cdf_lookup );
					if ( n_data == 0 )
						break;
					ccdf_lookup = 1.0 - cdf_lookup;
					ompr = 1.0 - pr;
					f_double = ( double ) f;
					s_double = ( double ) s;
					target->cumnbn (f_double, s_double, pr, ompr, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(ccdf_compute) - Math::Abs(ccdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-05);
				}

				n_data = -2;
				target->negative_binomial_cdf_values (n_data, f, s, pr, cdf_lookup );
				Assert::AreEqual(1, n_data);
			}
			/// <summary>
			///A test for test20
			///</summary>
	public: [TestMethod]
			void test20()
			{
				double ccdf_compute = 0.0;
				double ccdf_lookup = 0.0;
				double cdf_compute = 0.0;
				double cdf_lookup = 0.0;
				int n_data = 0;
				double x = 0.0;

				StatClass^ target = gcnew StatClass();

				n_data = 0;
				for ( ; ; )
				{
					target->normal_cdf_values(n_data, x, cdf_lookup);
					if ( n_data == 0 )
						break;
					target->cumnor (x, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(cdf_compute) - Math::Abs(cdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-06);
				}

				n_data = 0;
				for ( ; ; )
				{
					target->normal_cdf_values(n_data, x, cdf_lookup);
					if ( n_data == 0 )
						break;
					ccdf_lookup = 1.0 - cdf_lookup;
					target->cumnor (x, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(ccdf_compute) - Math::Abs(ccdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-06);
				}

				// Special values
				n_data = -2;
				target->normal_cdf_values(n_data, x, cdf_lookup);
				Assert::AreEqual(1, n_data);

			}
			/// <summary>
			///A test for test21
			///</summary>
	public: [TestMethod]
			void test21()
			{
				double ccdf_compute = 0.0;
				double ccdf_lookup = 0.0;
				double cdf_compute = 0.0;
				double cdf_lookup = 0.0;
				double lambda = 0.0;
				int n_data = 0;
				int x = 0;
				double x_double = 0.0;

				StatClass^ target = gcnew StatClass();

				n_data = 0;
				for ( ; ; )
				{
					target->poisson_cdf_values(n_data, lambda, x, cdf_lookup );
					if ( n_data == 0 )
						break;
					x_double = ( double ) x;
					target->cumpoi(x_double, lambda, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(cdf_compute) - Math::Abs(cdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-03);
				}

				n_data = 0;
				for ( ; ; )
				{
					target->poisson_cdf_values(n_data, lambda, x, cdf_lookup );
					if ( n_data == 0 )
						break;
					x_double = ( double ) x;
					ccdf_lookup = 1.0 - cdf_lookup;
					target->cumpoi(x_double, lambda, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(ccdf_compute) - Math::Abs(ccdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-03);
				}
				
				// Special values
				n_data = -2;
				target->poisson_cdf_values(n_data, lambda, x, cdf_lookup );
				Assert::AreEqual(1, n_data);
			}

			/// <summary>
			///A test for test22
			///</summary>
	public: [TestMethod]
			void test22()
			{
				double ccdf_compute = 0.0;
				double ccdf_lookup = 0.0;
				double cdf_compute = 0.0;
				double cdf_lookup = 0.0;
				int df = 0;
				double df_double = 0.0;
				int n_data = 0;
				double x = 0.0;

				StatClass^ target = gcnew StatClass();

				n_data = 0;
				for ( ; ; )
				{
					target->student_cdf_values (n_data, df, x, cdf_lookup );
					if ( n_data == 0 )
						break;
					df_double = ( double ) df;
					target->cumt (x,df_double,cdf_compute,ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(cdf_compute) - Math::Abs(cdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-03);
				}

				n_data = 0;
				for (; ; )
				{
					target->student_cdf_values (n_data, df, x, cdf_lookup );
					if ( n_data == 0 )
						break;
					ccdf_lookup = 1.0 - cdf_lookup;
					df_double = ( double ) df;
					target->cumt (x, df_double, cdf_compute, ccdf_compute );
					double diffDF = Math::Abs(Math::Abs(ccdf_compute) - Math::Abs(ccdf_lookup));
					Assert::IsTrue(diffDF < 9.0e-03);
				}

				// special values
				n_data = -2;
				target->student_cdf_values (n_data, df, x, cdf_lookup );
				Assert::AreEqual(1, n_data);

			}
			/// <summary>
			///A test for test23
			///</summary>
	public: [TestMethod]
			void test23()
			{
				double a = 2.2;
				double b = 3.7;
				double apb = a + b;
				double beta1 = 0.0;
				double beta2 = 0.0;

				StatClass^ target = gcnew StatClass();

				beta1 = target->beta ( a, b );
				beta2 = target->gamma_x ( a ) * target->gamma_x ( b ) / target->gamma_x ( apb );
				Assert::AreEqual(a, 2.2);
				Assert::AreEqual(b, 3.7);
				double diff = Math::Abs(Math::Abs(beta1) - Math::Abs(0.0454));
				Assert::IsTrue(diff < 9.0e-03);
				diff = Math::Abs(Math::Abs(beta2) - Math::Abs(0.0454));
				Assert::IsTrue(diff < 9.0e-03);

				a = -5.0;
				double y = target->gamma_x(a);
				diff = 0.041629017338776220;
				diff = Math::Abs(diff - Math::Abs(y));
				Assert::IsTrue(diff < 9.0e-09);

				a = 1001.7;
				y = target->gamma_x(a);
				Assert::IsTrue(y < 2.22e-16);

				a = 1.0e-31;
				y = target->gamma_x(a);
				Assert::AreEqual(9.9999999999999996e+030, y);
			
				a = -17.1;
				y = target->gamma_x(a);
				Assert::AreEqual(2.1461766635544670e-014, y);

				return;
			}
			/// <summary>
			///A test for test24
			///</summary>
	public: [TestMethod]
			void test24()
			{
				double erf_compute = 0.0;
				double erf_lookup = 0.0;
				double erfc_compute = 0.0;
				double erfc_lookup = 0.0;
				int ind = 0;
				int n_data = 0;
				double x = 0.0;

				StatClass^ target = gcnew StatClass();

				n_data = 0;
				for ( ; ; )
				{
					target->erf_values (n_data, x, erf_lookup );
					if ( n_data == 0 )
						break;
					erf_compute =target->error_f ( x );
					double diffDF = Math::Abs(Math::Abs(erf_compute) - Math::Abs(erf_lookup));
					Assert::IsTrue(diffDF < 9.0e-06);
				}

				ind = 0;
				n_data = 0;
				for ( ; ; )
				{
					target->erf_values (n_data, x, erf_lookup );
					if ( n_data == 0 )
						break;
					erfc_lookup = 1.0 - erf_lookup;
					erfc_compute = target->error_fc (ind, x );
					double diffDF = Math::Abs(Math::Abs(erfc_lookup) - Math::Abs(erfc_compute));
					Assert::IsTrue(diffDF < 9.0e-06);
				}
			}
			/// <summary>
			///A test for test25
			///</summary>
	public: [TestMethod]
			void test25()
			{
				double gamma_compute = 0.0;
				double gamma_lookup = 0.0;
				int n_data = 0;
				double x = 0.0;

				StatClass^ target = gcnew StatClass();

				for ( ; ; )
				{
					target->gamma_values ( n_data, x, gamma_lookup );
					if ( n_data == 0 )
						break;
					gamma_compute = target->gamma_x(x);
					double diffDF = Math::Abs(Math::Abs(gamma_lookup) - Math::Abs(gamma_compute))/Math::Abs(gamma_lookup);
					Assert::IsTrue(diffDF < 9.0e-06);
				}
				
				n_data = -2;
				target->gamma_values ( n_data, x, gamma_lookup );
				Assert::AreEqual(1, n_data);

			}
			/// <summary>
			///A test for test26
			///</summary>
	public: [TestMethod]
			void test26()
			{
//				double a = 3.0;
				int i = 0;
				int ierror = 0;
				int ind = 1;
				double p = 0.0;
				double q = 0.0;
				int test_num = 30;
				double x = 0.0;
				double x0 = 0.0;
				double x2 = 0.0;

				StatClass^ target = gcnew StatClass();

				for (double a = 0.1; a <= 50; a += 0.1)
				{
					for (i = 0; i <= test_num; i++)
					{
						// x = ( double ) i / ( double ) test_num;
						x = ( double ) i / 10.0;
						target->gamma_inc (a, x, p, q, ind );
						target->gamma_inc_inv (a, x2, x0, p, q, ierror);
						double diffpq = Math::Abs(Math::Abs(1.0) - Math::Abs(p+q));
						double diffxy = Math::Abs(Math::Abs(x) - Math::Abs(x2));
						Assert::IsTrue(diffpq < 9.0e-06);
						Assert::IsTrue(diffxy < 9.0e-06);
						x0 = x2;
						x2 = 0.0;
						target->gamma_inc_inv (a, x2, x0, p, q, ierror);
						x0 = 0;
						diffxy = Math::Abs(Math::Abs(x) - Math::Abs(x2));
						Assert::IsTrue(diffxy < 9.0e-04);
					}
				}
				
				// Special values
				double a = -2;
				target->gamma_inc_inv (a, x2, x0, p, q, ierror);
				Assert::AreEqual(-2, ierror);

				p = 1.05;
				q = 1.05;
				a = 1.0;
				target->gamma_inc_inv (a, x2, x0, p, q, ierror);
				Assert::AreEqual(-4, ierror);

				p = 1.0;
				q = 0.0;
				target->gamma_inc_inv (a, x2, x0, p, q, ierror);
				Assert::AreEqual(double::MaxValue, x2);

				p = 0.05;
				q = 0.95;
				a = 1.0;
				target->gamma_inc_inv (a, x2, x0, p, q, ierror);
				Assert::AreEqual(0.051293294387550543, x2);

				p = 0.15;
				q = 0.85;
				a = 1.0;
				target->gamma_inc_inv (a, x2, x0, p, q, ierror);
				Assert::AreEqual(0.16251892949777494, x2);

			}
			/// <summary>
			///A test for test27
			///</summary>
	public: [TestMethod]
			void test27()
			{
				double psi_compute = 0.0;
				double psi_lookup = 0.0;
				int n_data = 0;
				double x = 0.0;

				StatClass^ target = gcnew StatClass();

				n_data = 0;
				for (; ; )
				{
					target->psi_values (n_data, x, psi_lookup );
					if ( n_data == 0 )
						break;
					psi_compute = target->psi (x );
					double diffpsi = Math::Abs(Math::Abs(psi_lookup) - Math::Abs(psi_compute));
					Assert::IsTrue(diffpsi < 9.0e-06);
				}

				// special values
				x = 5.0e-10;
				psi_compute = target->psi (x );
				double diffpsi = Math::Abs(Math::Abs(psi_compute) - Math::Abs(-2000000000.5772154))/Math::Abs(psi_compute);
				Assert::IsTrue(diffpsi < 9.0e-06);

				x = 1.e-6;
				psi_compute = target->psi (x );
				diffpsi = Math::Abs(Math::Abs(psi_compute) - Math::Abs(-1000000.5772140201))/Math::Abs(psi_compute);
				Assert::IsTrue(diffpsi < 9.0e-06);

				x = 2147483625;
				psi_compute = target->psi (x );
				diffpsi = Math::Abs(Math::Abs(psi_compute) - Math::Abs(21.487562586415265))/Math::Abs(psi_compute);
				Assert::IsTrue(diffpsi < 9.0e-06);

				n_data = -2;
				target->psi_values (n_data, x, psi_lookup );
				Assert::AreEqual(1, n_data);
			}

			/// <summary>
			///A test for rcomp
			///</summary>
	public: [TestMethod]
			void rcomp()
			{
				double a = 30.0;
				double x = 0.0;

				StatClass^ target = gcnew StatClass();
				double y = target->rcomp(a, x);
				Assert::AreEqual(0.0, y);
			}

			///<summary>
			///A test for beta_asym
			///</summary>
	public: [TestMethod]
			void beta_asym()
			{
				double a = 15.5;
				double b = 15.8;
				double lamda = 6.1;
				double e = 1.e-6;

				StatClass^ target = gcnew StatClass();
				double y = target->beta_asym(a, b, lamda, e);
				Assert::AreEqual(0.011954744586439752, y);

				a = 17.5;
				b = 16.1;
				lamda = 9.75;
				e = 1.e-7;
				y = target->beta_asym(a, b, lamda, e);
				Assert::AreEqual(0.00014718137581028794, y);

				a = 32.7;
				b = 32.7;
				lamda = 6.5;
				e = 1.e-7;
				y = target->beta_asym(a, b, lamda, e);
				Assert::AreEqual(0.052877945173738956, y);
			}

			///<summary>
			///A test for dexpm1
			///</summary>
	public: [TestMethod]
			void dexpm1()
			{
				double x = 0.12;
				StatClass^ target = gcnew StatClass();
				
				double y =target->dexpm1(x);
				Assert::AreEqual(0.12749685157937568, y);

				x = -0.05;
				y = target->dexpm1(x);
				Assert::AreEqual(-0.048770575499285998, y);

				x = -2.3;
				y = target->dexpm1(x);
				Assert::AreEqual(-0.89974115627719620, y);

				x = 3.1;
				y = target->dexpm1(x);
				Assert::AreEqual(21.197951281441636, y);
			}

			///<summary>
			///A test for dstrem
			///</summary>
	public: [TestMethod]
			void dstrem()
			{

				double z = 4.35;
				StatClass^ target = gcnew StatClass();

				double y = target->dstrem(z);
				Assert::AreEqual(0.019123832248980843, y);

				z = 7.2;
				y = target->dstrem(z);
				Assert::AreEqual(0.011566672336202791, y);

				bool exceptionThrown = false;
				try
				{
					z = -2.1;
					y = target->dstrem(z);
				}

				catch (ArgumentException^)
				{
					exceptionThrown = true;
				}

				Assert::IsTrue(exceptionThrown);
			}

			///<summary>
			///A test for dlanor
			///</summary>
	public: [TestMethod]
			void dlanor()
			{

				double z = 5.35;
				StatClass^ target = gcnew StatClass();

				double y = target->dlanor(z);
				Assert::AreEqual(-16.939595936729798, y);

				z = 11.2;
				y = target->dlanor(z);
				Assert::AreEqual(-66.062671286925976, y);

				bool exceptionThrown = false;
				try
				{
					z = 4.2;
					y = target->dlanor(z);
				}

				catch (ArgumentException^)
				{
					exceptionThrown = true;
				}

				Assert::IsTrue(exceptionThrown);
			}

			///<summary>
			///A test for dbetrm
			///</summary>
	public: [TestMethod]
			void dbetrm()
			{

				double a = 5.35;
				double b = 2.87;
				StatClass^ target = gcnew StatClass();

				double y = target->dbetrm(a, b);
				Assert::AreEqual(0.034347722314398754, y);

				a = 2.87;
				b = 3.87;
				y = target->dbetrm(a, b);
				Assert::AreEqual(0.038053384639413830, y);

				bool exceptionThrown = false;
				try
				{
					a = -1.0;
					b = 0.0;
					y = target->dbetrm(a, b);
				}

				catch (ArgumentException^)
				{
					exceptionThrown = true;
				}

				Assert::IsTrue(exceptionThrown);
			}

			///<summary>
			///A test for esum
			///</summary>
	public: [TestMethod]
			void esum()
			{

				int mu = -2;
				double x = -0.05;
				StatClass^ target = gcnew StatClass();

				double y = target->esum(mu, x);
				Assert::AreEqual(0.12873490358780423, y);

				mu = -2;
				x = 0.05;
				y = target->esum(mu, x);
				Assert::AreEqual(0.14227407158651359, y);

				mu = 1;
				x = 2.1;
				y = target->esum(mu, x);
				Assert::AreEqual(22.197951281441636, y);
			}

	};
}
namespace DCDFlibUTest {
    
}
