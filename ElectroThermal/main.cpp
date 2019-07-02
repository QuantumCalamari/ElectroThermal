# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS

using namespace std;

int main();
void dtable_data_write(ofstream& output, int m, int n, double table[]);
void dtable_write(string output_filename, int m, int n, double table[],
	bool header);
void f(double a, double b, double t0, double t, int n, double x[],
	double value[]);
int r83_np_fa(int n, double a[]);
double* r83_np_sl(int n, double a_lu[], double b[], int job);
void timestamp();
void u0(double a, double b, double t0, int n, double x[], double value[]);
double ua(double a, double b, double t0, double t);
double ub(double a, double b, double t0, double t);
double sqr(double a);

//****************************************************************************80

int main()

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FD1D_HEAT_IMPLICIT.
//
//  Discussion:
//
//    FD1D_HEAT_IMPLICIT solves the 1D heat equation with an implicit method.
//
//    This program solves
//
//      dUdT - k * d2UdX2 = F(X,T)
//
//    over the interval [A,B] with boundary conditions
//
//      U(A,T) = UA(T),
//      U(B,T) = UB(T),
//
//    over the time interval [T0,T1] with initial conditions
//
//      U(X,T0) = U0(X)
//
//    The code uses the finite difference method to approximate the
//    second derivative in space, and an implicit backward Euler approximation
//    to the first derivative in time.
//
//    The finite difference form can be written as
//
//      U(X,T+dt) - U(X,T)                  ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) )
//      ------------------  = F(X,T) + k *  ------------------------------------
//               dt                                   dx * dx
//
//    so that we have the following linear system for the values of U at time T+dt:
//
//            -     k * dt / dx / dx   * U(X-dt,T+dt)
//      + ( 1 + 2 * k * dt / dx / dx ) * U(X,   T+dt)
//            -     k * dt / dx / dx   * U(X+dt,T+dt)
//      =               dt             * F(X,   T+dt)
//      +                                U(X,   T)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 May 2009
//
//  Author:
//
//    John Burkardt
//
{
	double* a;
	double* b;
	double* fvec;
	bool header;
	int i;
	int info;
	int j;
	int job;
	double k;
	double* t;
	double t_delt;
	string t_file;
	double t_max;
	double t_min;
	int t_num;
	double* u;
	string u_file;
	double w;
	double* x;
	double x_delt;
	string x_file;
	double x_max;
	double x_min;
	int x_num;

	timestamp();
	cout << "\n";
	cout << "FD1D_HEAT_IMPLICIT\n";
	cout << "  C++ version\n";
	cout << "\n";
	cout << "  Finite difference solution of\n";
	cout << "  the time dependent 1D heat equation\n";
	cout << "\n";
	cout << "    Ut - k * Uxx = F(x,t)\n";
	cout << "\n";
	cout << "  for space interval A <= X <= B with boundary conditions\n";
	cout << "\n";
	cout << "    U(A,t) = UA(t)\n";
	cout << "    U(B,t) = UB(t)\n";
	cout << "\n";
	cout << "  and time interval T0 <= T <= T1 with initial condition\n";
	cout << "\n";
	cout << "    U(X,T0) = U0(X).\n";
	cout << "\n";
	cout << "  A second order difference approximation is used for Uxx.\n";
	cout << "\n";
	cout << "  A first order backward Euler difference approximation\n";
	cout << "  is used for Ut.\n";

	clock_t begin = clock();

	//Lorenz Constant
	double lorenz = 2.45E-8;

	//kappa
	k = 5.0E-05;

	//Boltzmann Constant
	double kb = 1.38064852E-23;

	//applied current
	double i_app = 140E-6;
	//Ic(0)
	double i_c_abs = 320E-6;
	//double i_c = i_c_abs;
	//critical current
	double* i_c;
	double i_c_test;
	
	//resistance
	double r = 3.0;
	//normal state resistance
	double rn;

	//physical wire dimensions
	double b_width = 140E-9;
	double b_length = 15E-9;
	double len_seg;
	double d = 40E-9;

	//current density
	double j_den = i_app / (b_width * d);
	double j_c;// = i_c / (b_width * d);

	//resistance per square
	double rho = r * (b_width * d);

	//substrate temperature
	double t_sub = 2.0;
	//critical temperature
	double t_c = 4;
	//constants alpha and beta
	double alpha = 2E-103;
	double beta = alpha;


	//specific heat constants
	double gamma = 240.0;
	double cen;
	double ces;
	double delta = 2.1E-3 * 1.602176634E-19; //http://www.jetp.ac.ru/cgi-bin/dn/e_053_06_1270.pdf
	double c;
	double A_prop = 2.43 * gamma * t_c / exp(-delta / (kb * t_c));
	double nb_den = 8570;

	//radiative transfer term
	double rad_sub = 0;
	//heat generated from Joule heating
	double joule = 0;


	
	
	//
	//  Set X values.
	//
	x_min = 0.0;
	x_max = 0.3;
	x_num = 101;
	x_delt = (x_max - x_min) / (double)(x_num - 1);
	
	len_seg = b_length / x_num;

	i_c = new double[x_num];
	x = new double[x_num];

	for (i = 0; i < x_num; i++)
	{
		x[i] = ((double)(x_num - i - 1) * x_min
			+ (double)(i)* x_max)
			/ (double)(x_num - 1);
	}
	// 
	//  Set T values.
	//
	t_min = 0.0;
	t_max = 0.00001;
	t_num = 10001;
	t_delt = (t_max - t_min) / (double)(t_num - 1);

	//double istep = 10E-3 / (2 * (t_num - 1));

	t = new double[t_num];

	for (j = 0; j < t_num; j++)
	{
		t[j] = ((double)(t_num - j - 1) * t_min
			+ (double)(j)* t_max)
			/ (double)(t_num - 1);
	}
	//
	//  Set the initial data, for time T_MIN.
	//
	u = new double[x_num * t_num];

	u0(x_min, x_max, t_min, x_num, x, u);
	//
	//  The matrix A does not change with time.  We can set it once,
	//  factor it once, and solve repeatedly.
	//
	if (j == 0) {
		k = lorenz * u[i] / rho;



		if (u[i] < t_c)
			i_c[i] = i_c_abs * sqr(1 - (sqr(u[i] / t_c)));
		else
			i_c = 0;

		j_den = i_app / (b_width * d);
		
		if (u[i + j * x_num] < t_c || i_app < i_c[i]) {
			rho = 0;
			k *= (u[i] / t_c);
			c = A_prop * exp(-delta / (kb * u[i + j * x_num]));
		}
		else {
			rho = r * (b_width * d);
			c = gamma * u[i + j * x_num];
		}
	}
	else
		if (u[i + (j - 1) * x_num] > t_c) {
			//this is the normal state
			k = lorenz * u[i + (j - 1) * x_num] / rho;
			c = gamma * u[i + j * x_num];
		}
		else {
			//this is the superconducting state
			k = lorenz * u[i + (j - 1) * x_num] / rho * (u[i + (j - 1) * x_num] / t_c);
			//c = A_prop * exp(-delta / (kb * u[i + j * x_num]));
			//below is a temporary value for c
			c = 9800 + 2400;
		}

	k = 5E-1;

	w = k * t_delt / x_delt / x_delt;

	a = new double[3 * x_num];

	a[0 + 0 * 3] = 0.0;

	a[1 + 0 * 3] = 1.0;
	a[0 + 1 * 3] = 0.0;

	for (i = 1; i < x_num - 1; i++)
	{
		a[2 + (i - 1) * 3] = -w;
		a[1 + i * 3] = 1.0 + 2.0 * w;
		a[0 + (i + 1) * 3] = -w;
	}

	a[2 + (x_num - 2) * 3] = 0.0;
	a[1 + (x_num - 1) * 3] = 1.0;

	a[2 + (x_num - 1) * 3] = 0.0;
	//
	//  Factor the matrix.
	//
	info = r83_np_fa(x_num, a);

	b = new double[x_num];
	fvec = new double[x_num];

	for (j = 1; j < t_num; j++)
	{
		//
		//  Set the right hand side B.
		//
		b[0] = ua(x_min, x_max, t_min, t[j]);

		f(x_min, x_max, t_min, t[j - 1], x_num, x, fvec);

		for (i = 1; i < x_num - 1; i++)
		{
			b[i] = u[i + (j - 1) * x_num] + t_delt * fvec[i];
		}

		b[x_num - 1] = ub(x_min, x_max, t_min, t[j]);

		delete[] fvec;

		job = 0;
		fvec = r83_np_sl(x_num, a, b, job);

		//don't touch anything above this line ------------------------------------------------------------------------------------------------------------------------------------------------------------------

		for (i = 0; i < x_num; i++)
		{
			
			if (j == 0) {
				i_c[i] = 0;
			}
			else {
				//checks to see if superconductor or normal
				if (u[i + (j-1) * x_num] < t_c)
					//normal state
					i_c[i] = i_c_abs * sqr(1 - (sqr(u[i + (j-1) * x_num] / t_c)));
				else
					//superconducting state
					i_c[i] = 0;

				j_den = i_app / (b_width * d);
				//k = lorenz * 
				if (u[i + (j-1) * x_num] < t_c && i_app < i_c[i])
					rho = 0;
				else
					rho = r*(b_width * d)/ len_seg;
			}

			alpha = beta * sqr(u[i + (j - 1) * x_num]);
			
			//calculation of joule and radiative transfer
			
			joule = sqr(j_den) * rho * t_delt * (b_width * d * len_seg) / (c * nb_den * (b_width * d * len_seg));
			//joule = i_app * r * t_delt;
			rad_sub = alpha / d * (u[i + (j - 1) * x_num] - t_sub) * t_delt;

			u[i + j * x_num] = fvec[i] - rad_sub + joule;
		}
	}

	//don't touch anything below here -----------------------------------------------------------------------------------------------------------------------------------------------------

	x_file = "x.txt";
	header = false;
	dtable_write(x_file, 1, x_num, x, header);

	cout << "\n";
	cout << "  X data written to \"" << x_file << "\".\n";

	t_file = "t.txt";
	header = false;
	dtable_write(t_file, 1, t_num, t, header);

	cout << "  T data written to \"" << t_file << "\".\n";

	u_file = "u.txt";
	header = false;
	dtable_write(u_file, x_num, t_num, u, header);

	cout << "  U data written to \"" << u_file << "\".\n";

	cout << "\n";
	cout << "FD1D_HEAT_IMPLICIT\n";
	cout << "  Normal end of execution.\n";
	cout << "\n";
	timestamp();

	delete[] a;
	delete[] b;
	delete[] fvec;
	delete[] t;
	delete[] u;
	delete[] x;

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	std::cout << "this took " << elapsed_secs << " seconds" << std::endl;

	return 0;
}
//****************************************************************************80

void dtable_data_write(ofstream& output, int m, int n, double table[])
{
	int i;
	int j;

	for (j = 0; j < n; j++)
	{
		for (i = 0; i < m; i++)
		{
			output << setw(10) << table[i + j * m] << "  ";
		}
		output << "\n";
	}

	return;
}
//****************************************************************************80

void dtable_write(string output_filename, int m, int n, double table[],
	bool header)
{
	ofstream output;

	output.open(output_filename.c_str());

	if (!output)
	{
		cerr << "\n";
		cerr << "DTABLE_WRITE - Fatal error!\n";
		cerr << "  Could not open the output file.\n";
		return;
	}

	if (header)
	{
		//  dtable_header_write ( output_filename, output, m, n );
	}

	dtable_data_write(output, m, n, table);

	output.close();

	return;
}
//****************************************************************************80

void f(double a, double b, double t0, double t, int n, double x[],
	double value[])
{
	int i;

	for (i = 0; i < n; i++)
	{
		value[i] = 0.0;
	}
	return;
}
//****************************************************************************80

int r83_np_fa(int n, double a[])
{
	int i;

	for (i = 1; i <= n - 1; i++)
	{
		if (a[1 + (i - 1) * 3] == 0.0)
		{
			cout << "\n";
			cout << "R83_NP_FA - Fatal error!\n";
			cout << "  Zero pivot on step " << i << "\n";
			return i;
		}
		//
		//  Store the multiplier in L.
		//
		a[2 + (i - 1) * 3] = a[2 + (i - 1) * 3] / a[1 + (i - 1) * 3];
		//
		//  Modify the diagonal entry in the next column.
		//
		a[1 + i * 3] = a[1 + i * 3] - a[2 + (i - 1) * 3] * a[0 + i * 3];
	}

	if (a[1 + (n - 1) * 3] == 0.0)
	{
		cout << "\n";
		cout << "R83_NP_FA - Fatal error!\n";
		cout << "  Zero pivot on step " << n << "\n";
		return n;
	}

	return 0;
}
//****************************************************************************80

double* r83_np_sl(int n, double a_lu[], double b[], int job)
{
	int i;
	double* x;

	x = new double[n];

	for (i = 0; i < n; i++)
	{
		x[i] = b[i];
	}

	if (job == 0)
	{
		//
		//  Solve L * Y = B.
		//
		for (i = 1; i < n; i++)
		{
			x[i] = x[i] - a_lu[2 + (i - 1) * 3] * x[i - 1];
		}
		//
		//  Solve U * X = Y.
		//
		for (i = n; 1 <= i; i--)
		{
			x[i - 1] = x[i - 1] / a_lu[1 + (i - 1) * 3];
			if (1 < i)
			{
				x[i - 2] = x[i - 2] - a_lu[0 + (i - 1) * 3] * x[i - 1];
			}
		}
	}
	else
	{
		//
		//  Solve U' * Y = B
		//
		for (i = 1; i <= n; i++)
		{
			x[i - 1] = x[i - 1] / a_lu[1 + (i - 1) * 3];
			if (i < n)
			{
				x[i] = x[i] - a_lu[0 + i * 3] * x[i - 1];
			}
		}
		//
		//  Solve L' * X = Y.
		//
		for (i = n - 1; 1 <= i; i--)
		{
			x[i - 1] = x[i - 1] - a_lu[2 + (i - 1) * 3] * x[i];
		}
	}

	return x;
}
//****************************************************************************80

void timestamp()
{
# define TIME_SIZE 40

	static char time_buffer[TIME_SIZE];
	const struct tm* tm;
	size_t len;
	time_t now;

	now = time(NULL);
	tm = localtime(&now);

	len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

	cout << time_buffer << "\n";

	return;
# undef TIME_SIZE
}
//****************************************************************************80

void u0(double a, double b, double t0, int n, double x[], double value[])
{
	int i;

	for (i = 0; i < n; i++)
	{
		value[i] = 2.0;
	}

	value[50] = 12;
	
	return;
}
//****************************************************************************80

double ua(double a, double b, double t0, double t)
{
	double value;

	value = 2.0;

	return value;
}
//****************************************************************************80

double ub(double a, double b, double t0, double t)
{
	double value;

	value = 2.0;
	
	return value;
}
//****************************************************************************80

double sqr(double a) {

	return a * a;
}