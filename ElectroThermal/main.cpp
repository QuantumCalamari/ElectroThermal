# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <string>
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
double cube(double a);

//****************************************************************************80

int main()

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
	string r_file;
	double* rt;

	timestamp();
	cout << "\n";
	cout << "Electro-Thermal simulation of 1D superconducting wire\n";
	cout << "\n";

	clock_t begin = clock();

	//Lorenz Constant
	double lorenz = 2.45E-8;

	//kappa
	k = 5.0E-05;

	//Boltzmann Constant
	double kb = 1.38064852E-23;

	//applied current
	double i_app = 20E-6;
	//Ic(0)
	double i_c_abs = 320E-6;
	//double i_c = i_c_abs;
	//critical current
	double i_c;
	double j_c;
	double* i_str;
	string i_file;
	double delta_i = 0;
	double delta_i_max = 0.2;

	double i_wire = i_app;
	double i_prev;
	
	double z0 = 50;

	//resistance per square
	double r = 600;
	//normal state resistance
	double rn;

	//physical wire dimensions
	double b_width = 100E-9;
	double b_length = 10E-6;
	double len_seg;
	double d = 4E-9;
	double* Lk;
	string l_file;

	//photon parameters
	int p_t = 50;
	double p_lambda = 1550E-9;
	double p_e = 0;

	//current density
	double j_den = i_app / (b_width * d);
	//double j_c;// = i_c / (b_width * d);
	double js;
	double jn;

	//resistance per square
	double rho = r * (b_width * d);
	double r0 = 600;
	double t0 = 300;
	double a0 = -0.0013986;

	//substrate temperature
	double t_sub = 2.0;
	//critical temperature
	double t_c = 4;
	//constants alpha and beta
	double alpha = 2E3;
	double beta = alpha;

	//specific heat constants
	double gamma = 240.0;
	double sc_band_gap = 2E-3;
	double ev_con = 1.602176634E-19;
	double cen;
	double ces;
	double delta = 2.1E-3 * 1.602176634E-19; //http://www.jetp.ac.ru/cgi-bin/dn/e_053_06_1270.pdf
	double c;
	double A_prop = 2.43 * gamma * t_c / exp(-delta / (kb * t_c));
	double nb_den = 8570;
	double B = 8E2; //calculated from the Yang paper, alpha = 8E5 at 10 K

	//radiative transfer term
	double rad_sub = 0;
	//heat generated from Joule heating
	double joule = 0;


	
	
	//
	//  Set X values.
	//
	x_min = 0.0;
	x_max = 0.3;
	x_num = 201;
	x_delt = (x_max - x_min) / (double)(x_num - 1);
	
	len_seg = b_length / x_num;

	//i_c = new double[x_num];
	x = new double[x_num];

	for (i = 0; i < x_num; i++)
	{
		x[i] = ((double)(x_num - i - 1) * x_min + (double)(i)* x_max) / (double)(x_num - 1);
	}
	// 
	//  Set T values.
	//
	t_min = 0.0;
	t_max = 1E-5;
	t_num = 10001;
	t_delt = (t_max - t_min) / (double)(t_num - 1);

	//double istep = 10E-3 / (2 * (t_num - 1));

	t = new double[t_num];
	rt = new double[t_num];
	i_str = new double[t_num];
	Lk = new double[t_num];

	for (j = 0; j < t_num; j++)
	{
		t[j] = ((double)(t_num - j - 1) * t_min
			+ (double)(j)* t_max)
			/ (double)(t_num - 1);
		rt[j] = 0;
		i_str[j] = 0;
		Lk[j] = 0;
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
			i_c = i_c_abs * sqr(1 - (sqr(u[i] / t_c)));
		else
			i_c = 0;

		//I should move this into a function -------------
		rn = 0;
		i_wire = 0;
		

		for (int r_sum = 0; r_sum < x_num; r_sum++) {
			//need to include rho term
			rn += r * (b_width * d);
		}

		if (rn == 0) {
			i_wire = i_app;
		}
		else {
			i_wire = z0 / rn * i_app * (1 / (1 + z0 / rn));
		}
		// --------- down to here, all function

		j_den = i_wire / (b_width * d);

		if (u[i + j * x_num] < t_c && i_wire < i_c) {
			rho = 0;
			k = lorenz * sqr(u[i + (j - 1) * x_num]) / (r0 * (1 + a0 * (u[i + (j - 1) * x_num] - t0)) * d * (u[i + (j - 1) * x_num] / t_c));
			c = A_prop * exp(-delta / (kb * u[i]));
		}
		else {
			rho = r0 * (1 + a0 * (u[i + (j - 1) * x_num] - t0)) * d;
			c = gamma * u[i + j * x_num];
			k = lorenz * u[i] / rho;
		}
	}
	else
		if (u[i] > t_c) {
			//this is the normal state
			k = lorenz * u[i] / rho;
			c = gamma * u[i + j * x_num];
		}
		else {
			//this is the superconducting state
			k = lorenz * u[i] / rho;

			//k = lorenz * sqr(u[i + (j - 1) * x_num]) / (r * (b_width * d) * (u[i + (j - 1) * x_num] / t_c));
			c = A_prop * exp(-delta / (kb * u[i]));
		}

	k = 2E2;

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
				for (i = 0; i < x_num; i++) {
					i_c = 0;
				}
			}

			if (j > 0) {

				if (j == p_t && i == (x_num - 1) / 2) {
					//p_e = 6.626E-34 / p_lambda / (A_prop * exp(-delta / (kb * u[i + (j - 1) * x_num])) * b_width * d * len_seg);
					p_e = 12;
				}
				else {
					p_e = 0;
				}



				if (rt[j - 1] == 0)
					i_wire = i_app;
				else {

					i_wire = z0 / rt[j - 1] * i_app * (1 / (1 + z0 / rt[j - 1]));
					delta_i = i_wire - i_str[j - 1];

					if (i_str[j - 1] > i_wire) {
						if (i_str[j - 1] * (1 - delta_i) > i_wire)
							i_wire = i_str[j - 1] * (1 - delta_i);
					}

					if (i_str[j - 1] < i_wire) {
						if (i_str[j - 1] * (1 + delta_i) < i_wire)
							i_wire = i_str[j - 1] * (1 + delta_i);
					}
					
					if (i_str[j - 1] == 0) {
						i_wire = i_app * 0.5;
					}
					/*
					else if (delta_i > 0 && delta_i > abs(i_str[j - 1] * (1 + delta_i_max))) {
						i_wire = i_str[j - 1] * (1 + delta_i_max);
					}
					else if (delta_i < 0 && abs(delta_i) > abs(i_str[j - 1]) * (1 - delta_i_max)) {
						i_wire = i_str[j - 1] * (1 - delta_i_max);
					}*/
				}

				if (i_wire > i_app)
					i_wire = i_app;

				i_str[j] = i_wire;
				//temp delta_i code
				delta_i = i_wire - i_str[j - 1];

				if (u[i + (j - 1) * x_num] < t_c)
					i_c = i_c_abs * sqr(1 - (sqr(u[i + (j - 1) * x_num] / t_c)));
				else
					i_c = 0;

				j_c = i_c / (b_width * d);

				

				if (u[i + (j - 1) * x_num] < t_c && i_wire < i_c) {
					//superconducting state
					c = A_prop * exp(-delta / (kb * u[i + (j - 1) * x_num]));
					rho = 0;

					jn = j_den * sqr(sqr((u[i + (j - 1) * x_num] / t_c)));
					js = j_den - jn;

					Lk[j] += 9.10938356E-31 / (2 * js * 1.60217662E-19) * (len_seg / (d * b_width));
				}
				else {
					//normal state
					c = gamma * u[i + (j - 1) * x_num];
					
					if (u[i + (j - 1) * x_num] < t0)
						rho = r0 * (1 + a0 * (u[i + (j - 1) * x_num] - t0)) / (b_width/len_seg);
					else
						rho = r0;

					if (rho < 0) {

						std::cout << "seriously but, how the fuck is it negative" << std::endl;
						rho = r0 * (1 + a0 * (u[i + (j - 1) * x_num] - t0));
					}

					rt[j] += rho;

					jn = j_den;
					js = 0;

					Lk[j] += 9.10938356E-31 / (jn * 1.60217662E-19) * (len_seg / (d * b_width));
				}
			}


			j_den = i_wire / (b_width * d);

			joule = (sqr(j_den) * rho * t_delt) / c;

			if (joule > u[i + (j - 1) * x_num] * 1.1)
				joule = u[i + (j - 1) * x_num] * 0.1;

			if (joule > 0.03)
				joule = 0.03;

			//c = 12200;

			alpha = B * cube(u[i + (j - 1) * x_num]);
			//rad_sub = alpha / d * (u[i + (j - 1) * x_num] - t_sub) * t_delt / (c * nb_den); //needs a heat capacitance and a mass term somewhere
			rad_sub = (alpha / d * (fvec[i] - t_sub)) / (c + 9.8 * cube(fvec[i])) * t_delt;
			rad_sub = 0;
			

			//if (abs(rad_sub) > abs(u[i + (j - 1) * x_num] * 1.8))
			//	rad_sub = u[i + (j - 1) * x_num] * 1.8;

			u[i + j * x_num] = fvec[i] - rad_sub + joule + p_e;

			if (u[i + j * x_num] < 0) {
				u[i + j * x_num] = u[i + (j - 1) * x_num];
			}
		}
	}

	//don't touch anything below here -----------------------------------------------------------------------------------------------------------------------------------------------------

	l_file = "l.txt";
	header = false;
	dtable_write(l_file, 1, t_num, Lk, header);

	cout << "\n";
	cout << "  Inductance data written to \"" << l_file << "\".\n";

	r_file = "r.txt";
	header = false;
	dtable_write(r_file, 1, t_num, rt, header);

	cout << "  Resistance data written to \"" << r_file << "\".\n";

	i_file = "i.txt";
	header = false;
	dtable_write(i_file, 1, t_num, i_str, header);

	cout << "  Current data written to \"" << i_file << "\".\n";

	x_file = "x.txt";
	header = false;
	dtable_write(x_file, 1, x_num, x, header);

	cout << "  Positional data written to \"" << x_file << "\".\n";

	t_file = "t.txt";
	header = false;
	dtable_write(t_file, 1, t_num, t, header);

	cout << "  Time data written to \"" << t_file << "\".\n";

	u_file = "u.txt";
	header = false;
	dtable_write(u_file, x_num, t_num, u, header);

	cout << "  Thermal data written to \"" << u_file << "\".\n";

	cout << "\n";
	cout << "  Normal end of execution.\n";
	cout << "\n";
	timestamp();

	delete[] a;
	delete[] b;
	delete[] fvec;
	delete[] t;
	delete[] u;
	delete[] x;
	delete[] rt;
	delete[] i_str;
	delete[] Lk;

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	std::cout << "Calculation time was " << elapsed_secs << " seconds" << std::endl;

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

double cube(double a) {
	return a * a * a;
}