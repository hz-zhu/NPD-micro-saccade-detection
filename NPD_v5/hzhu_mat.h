#pragma once

#ifndef hzhu_mat_h
#define hzhu_mat_h

#include <gsl/gsl_matrix.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sort.h>

#ifndef HZHU_WARNING_DISP
#define HZHU_WARNING_DISP 1
#endif // !HZHU_WARNING_DISP

#ifndef HZHU_WARNING_EXIT
#define HZHU_WARNING_EXIT 0
#endif // !HZHU_WARNING_EXIT


class hzhu_mat;
class hzhu_mat_LU_inv;
class hzhu_qr_flat;

class hzhu_mat
{
public:
	hzhu_mat();
	hzhu_mat(double);
	hzhu_mat(int);
	hzhu_mat(int, int);
	hzhu_mat(size_t, size_t);
	hzhu_mat(double, int, int);
	hzhu_mat(double, int);
	hzhu_mat(std::string, std::string);
	hzhu_mat(std::string);
	hzhu_mat(hzhu_mat &);
	hzhu_mat(const hzhu_mat &);
	hzhu_mat(const gsl_matrix *);
	
public://constructors for random number generation
	hzhu_mat(gsl_rng *r);// one standard normal variable
	hzhu_mat(gsl_rng *r, int size1, int size2);// A matrix of standard normal variables
	hzhu_mat(gsl_rng *r, int size1, int size2, double mean, double sd);// A matrix of normal variables N(mean,sd^2)

	// A matrix of type variables, supported type include uniform, gaussian, log normal, exponential, Laplace, Exponential Power, Cauchy, Rayleigh, Gamma 
	hzhu_mat(gsl_rng *r, std::string type, int size1, int size2, double p1, double p2);
	

	~hzhu_mat();

public:
	gsl_matrix *data = NULL;

public:
	void disp();
	void disp(int);

	// Overload operators
public:
	void operator*=(hzhu_mat a);
	void operator/=(hzhu_mat a);
	void operator+=(hzhu_mat a);
	void operator-=(hzhu_mat a);
	void operator^=(hzhu_mat a);
	void operator*=(double a);
	void operator/=(double a);
	void operator+=(double a);
	void operator-=(double a);
	void operator^=(double a);
	
	hzhu_mat operator*(hzhu_mat a);
	hzhu_mat operator/(hzhu_mat a);
	hzhu_mat operator+(hzhu_mat a);
	hzhu_mat operator-(hzhu_mat a);
	hzhu_mat operator^(hzhu_mat a);

	hzhu_mat operator*(double a);
	hzhu_mat operator/(double a);
	hzhu_mat operator+(double a);
	hzhu_mat operator-(double a);
	hzhu_mat operator^(double a);

	double operator[](int);
	double operator[](hzhu_mat);

	// Operations
public:
	hzhu_mat transpose();
	void transpose_self();
	hzhu_mat inv();// the inverse is very inaccurate
	hzhu_mat pinv();// the left inverse: a.pinv()*a=I
	hzhu_mat product(hzhu_mat a);

	// IO
public:
	std::string print_to_string(int precision);
	int save_to_file(std::string name, int precision);
	int save_to_file(std::string name, std::string path, int precision);
	int read_from_file(std::string name);
	int read_from_file(std::string name, std::string path);

	// Statisitcs
public:
	hzhu_mat max_index_value_row();
	hzhu_mat max_index_row();
	hzhu_mat max_value_row();
	hzhu_mat max_index_value_col();
	hzhu_mat max_index_col();
	hzhu_mat max_value_col();

	hzhu_mat max_index_value_all();
	hzhu_mat max_value_all();
	hzhu_mat max_index_all();

	hzhu_mat min_index_value_row();
	hzhu_mat min_index_row();
	hzhu_mat min_value_row();
	hzhu_mat min_index_value_col();
	hzhu_mat min_index_col();
	hzhu_mat min_value_col();
	hzhu_mat min_index_value_all();
	hzhu_mat min_value_all();
	hzhu_mat min_index_all();

	hzhu_mat sum_row();
	hzhu_mat sum_col();
	hzhu_mat sum_all();

	hzhu_mat mean_row();
	hzhu_mat mean_col();
	hzhu_mat mean_all();

	hzhu_mat var_row();
	hzhu_mat var_col();
	hzhu_mat var_all();

	hzhu_mat uvar_row();
	hzhu_mat uvar_col();
	hzhu_mat uvar_all();

	hzhu_mat median_row();
	hzhu_mat median_col();
	hzhu_mat median_all();

	hzhu_mat quantile_row(double q);
	hzhu_mat quantile_col(double q);
	hzhu_mat quantile_all(double q);

	// Access elements
	hzhu_mat get_row(int);
	hzhu_mat get_col(int);
	hzhu_mat get_sub(int, int, int, int);

	// trace
	double tr();

	// Dimension change append
	void append_left(double a);
	void append_right(double a);
	void append_top(double a);
	void append_bottom(double a);
	void append_left(hzhu_mat&);
	void append_right(hzhu_mat&);
	void append_top(hzhu_mat&);
	void append_bottom(hzhu_mat&);

	// Generate random number
	//void rand_gaussian(gsl_rng *);
	//void rand_gaussian(gsl_rng *, double mean, double var)

	// Other functions
	hzhu_mat diff(int n);
	hzhu_mat diff();
};

class hzhu_mat_LU_inv
{
public:
	int signum;
	double det;
	double lndet;
	gsl_permutation *p = NULL;
	gsl_matrix *LU = NULL;
	gsl_matrix *A = NULL;
	gsl_matrix *A_inv = NULL;

	hzhu_mat_LU_inv(hzhu_mat &);
	~hzhu_mat_LU_inv();
	hzhu_mat_LU_inv(const hzhu_mat_LU_inv &);

public:
	void compute_LU();
	void compute_det();
	void compute_lndet();
	double get_det();
	double get_lndet();
	void compute_inv();
	void compute_inv_fine(int);
	hzhu_mat get_inv();
};

class hzhu_qr_flat
{
public:
	hzhu_qr_flat(int);
	hzhu_qr_flat();
	~hzhu_qr_flat();

public:
	int d;
	gsl_qrng * qr;

public:
	hzhu_mat generate(int length);
};

// Other Function Declearation =============================================================================================================
hzhu_mat hzhu_line_space(double, double, int);
hzhu_mat hzhu_equal_space(double, double, double);
hzhu_mat hzhu_mat_index(int x, int y);
hzhu_mat hzhu_mat_repmat(hzhu_mat, int, int);

// hzhu_qr_flat class function realization ============================================================================================================

hzhu_qr_flat::hzhu_qr_flat(int n)
{
	d = n;
	qr = gsl_qrng_alloc(gsl_qrng_reversehalton, d);
}

hzhu_qr_flat::hzhu_qr_flat()
{
	d = 1;
	qr = gsl_qrng_alloc(gsl_qrng_reversehalton, d);
}

hzhu_qr_flat::~hzhu_qr_flat()
{
	gsl_qrng_free(qr);
}

hzhu_mat hzhu_qr_flat::generate(int length)
{
	hzhu_mat re(length, d);
	double *holder = new double[d];

	int j;
	for (int i = 0; i < length; i++)
	{
		gsl_qrng_get(qr, holder);
		for (j = 0; j < d; j++) re.data->data[i*d + j] = holder[j];
	}
	delete holder;

	return re;
}

// hzhu_mat_inv Class Function Realization =============================================================================================================
hzhu_mat_LU_inv::hzhu_mat_LU_inv(hzhu_mat &m)
{
	if (m.data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL matrix @ hzhu_mat_LU_inv::hzhu_mat_LU_inv(hzhu_mat &m)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	if (m.data->size1 != m.data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Non-square Matrix @ hzhu_mat_LU_inv::hzhu_mat_LU_inv(hzhu_mat &m)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	int n = m.data->size1;
	A = gsl_matrix_alloc(n, n);
	gsl_matrix_memcpy(A, m.data);

	signum = GSL_NAN;
	det = GSL_NAN;
	lndet = GSL_NAN;
}

hzhu_mat_LU_inv::~hzhu_mat_LU_inv()
{
	if (p != NULL) gsl_permutation_free(p);
	if (LU != NULL) gsl_matrix_free(LU);
	if (A != NULL) gsl_matrix_free(A);
	if (A_inv != NULL) gsl_matrix_free(A_inv);
}

hzhu_mat_LU_inv::hzhu_mat_LU_inv(const hzhu_mat_LU_inv &s)
{
	if (s.p != NULL)
	{
		p = gsl_permutation_alloc(s.p->size);
		gsl_permutation_memcpy(p, s.p);
	}
	if (s.LU != NULL)
	{
		LU = gsl_matrix_alloc(s.LU->size1, s.LU->size2);
		gsl_matrix_memcpy(LU, s.LU);
	}
	if (s.A != NULL)
	{
		A = gsl_matrix_alloc(s.A->size1, s.A->size2);
		gsl_matrix_memcpy(A, s.A);
	}
	if (s.A_inv != NULL)
	{
		A_inv = gsl_matrix_alloc(s.A_inv->size1, s.A_inv->size2);
		gsl_matrix_memcpy(A_inv, s.A_inv);
	}
	det = s.det;
	signum = s.signum;
	lndet = s.lndet;
}

void hzhu_mat_LU_inv::compute_LU()
{
	if (A == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL matrix @ void hzhu_mat_LU_inv::compute_LU()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	if (A->size1!=A->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Non-square matrix @ void hzhu_mat_LU_inv::compute_LU()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	if (LU != NULL && p != NULL) return;
	if (LU != NULL) gsl_matrix_free(LU);
	if (p != NULL) gsl_permutation_free(p);
	LU = gsl_matrix_alloc(A->size1, A->size2);
	gsl_matrix_memcpy(LU, A);
	p = gsl_permutation_alloc(A->size1);
	gsl_linalg_LU_decomp(LU, p, &signum);
}

void hzhu_mat_LU_inv::compute_det()
{
	if (A == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL matrix @ void hzhu_mat_LU_inv::compute_LU()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	if (A->size1 != A->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Non-square matrix @ void hzhu_mat_LU_inv::compute_LU()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	if (LU == NULL || p == NULL)
	{
		this->compute_LU();
	}
	det = gsl_linalg_LU_det(LU, signum);
}

double hzhu_mat_LU_inv::get_det()
{
	if (A == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL matrix @ double hzhu_mat_LU_inv::get_det()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return GSL_NAN;
	}
	if (A->size1 != A->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Non-square matrix @ double hzhu_mat_LU_inv::get_det()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return GSL_NAN;
	}

	if (det != GSL_NAN && p != NULL && LU != NULL) return det;
	else
	{
		this->compute_det();
		return det;
	}
}

void hzhu_mat_LU_inv::compute_lndet()
{
	if (A == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL matrix @ void hzhu_mat_LU_inv::compute_lndet()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	if (A->size1 != A->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Non-square matrix @ void hzhu_mat_LU_inv::compute_lndet()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	if (LU == NULL || p == NULL)
	{
		this->compute_LU();
	}
	lndet = gsl_linalg_LU_lndet(LU);
}

double hzhu_mat_LU_inv::get_lndet()
{
	if (A == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL matrix @ double hzhu_mat_LU_inv::get_lndet()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return GSL_NAN;
	}
	if (A->size1 != A->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Non-square matrix @ double hzhu_mat_LU_inv::get_lndet()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return GSL_NAN;
	}

	if (det != GSL_NAN && p != NULL && LU != NULL) return lndet;
	else
	{
		this->compute_lndet();
		return lndet;
	}
}

void hzhu_mat_LU_inv::compute_inv()
{
	if (A == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL matrix @ void hzhu_mat_LU_inv::compute_inv()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	if (A->size1 != A->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Non-square matrix @ void hzhu_mat_LU_inv::compute_inv()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}

	if (A_inv != NULL) return;

	if (p == NULL || LU == NULL)
	{
		this->compute_LU();
	}
	int n = A->size1;
	A_inv = gsl_matrix_alloc(n, n);

	gsl_vector *x = gsl_vector_alloc(n);
	for (int i = 0; i < n; i++)
	{
		gsl_vector_set_all(x, 0.0);
		x->data[i] = 1.0;
		gsl_linalg_LU_svx(LU, p, x);
		gsl_matrix_set_col(A_inv, i, x);
	}
	gsl_vector_free(x);
}

hzhu_mat hzhu_mat_LU_inv::get_inv()
{
	if (A == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL matrix @ void hzhu_mat_LU_inv::compute_inv()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	if (A->size1 != A->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Non-square matrix @ void hzhu_mat_LU_inv::compute_inv()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	if (A_inv != NULL) return hzhu_mat(A_inv);
	else
	{
		this->compute_inv();
		return hzhu_mat(A_inv);
	}
}

void hzhu_mat_LU_inv::compute_inv_fine(int steps = 5)
{
	if (A == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL matrix @ void hzhu_mat_LU_inv::compute_inv()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	if (A->size1 != A->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Non-square matrix @ void hzhu_mat_LU_inv::compute_inv()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}

	if (A_inv != NULL) return;

	if (p == NULL || LU == NULL)
	{
		this->compute_LU();
	}
	int n = A->size1;
	A_inv = gsl_matrix_alloc(n, n);

	gsl_vector *x = gsl_vector_alloc(n);
	gsl_vector *work = gsl_vector_alloc(n);
	gsl_vector *b = gsl_vector_alloc(n);
	int j;
	for (int i = 0; i < n; i++)
	{
		gsl_vector_set_all(b, 0.0);
		b->data[i] = 1.0;
		gsl_vector_set_all(x, 0.0);
		x->data[i] = 1.0;
		gsl_linalg_LU_svx(LU, p, x);
		for (j = 0; j < steps; j++) gsl_linalg_LU_refine(A, LU, p, b, x, work);
		gsl_matrix_set_col(A_inv, i, x);
	}

	gsl_vector_free(x);
	gsl_vector_free(work);
	gsl_vector_free(b);
}



// hzhu_mat Class Function Realization =============================================================================================================

// Constructor ------------------------------------------------------------------------------------------------------------
hzhu_mat::hzhu_mat()
{

}

hzhu_mat::hzhu_mat(std::string name, std::string path)
{
	if (read_from_file(name, path)==0)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Cannot read file @ hzhu_mat::hzhu_mat(std::string name, std::string path)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
}

hzhu_mat::hzhu_mat(std::string name)
{
	if (read_from_file(name) == 0)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Cannot read file @ hzhu_mat::hzhu_mat(std::string name)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
}

hzhu_mat::hzhu_mat(double a)
{
	data = gsl_matrix_alloc(1, 1);
	gsl_matrix_set_all(data, a);
}

hzhu_mat::hzhu_mat(double a, int n1, int n2)
{
	data = gsl_matrix_alloc(abs(n1), abs(n2));
	gsl_matrix_set_all(data, a);
}

hzhu_mat::hzhu_mat(double a, int n)
{
	int N = abs(n);
	data = gsl_matrix_alloc(N, N);
	gsl_matrix_set_all(data, a);
}

hzhu_mat::hzhu_mat(int n)
{
	int N = abs(n);
	data = gsl_matrix_alloc(N, N);
}

hzhu_mat::hzhu_mat(int n1, int n2)
{
	data = gsl_matrix_alloc(abs(n1), abs(n2));
}

hzhu_mat::hzhu_mat(size_t n1, size_t n2)
{
	data = gsl_matrix_alloc(abs((int)n1), abs((int)n2));
}

hzhu_mat::~hzhu_mat()
{
	if (data != NULL)
	{
		gsl_matrix_free(data);
	}
}

void hzhu_mat::disp()
{
	if (data != NULL)
	{
		std::cout << "  A " << data->size1 << " x " << data->size2 << " Matrix" << std::endl;
		for (int i = 0; i < data->size1; i++)
		{
			for (int j = 0; j < data->size2; j++)
				std::cout << gsl_matrix_get(data, i, j) << "\t";
			std::cout << std::endl;
		}
	}
	else
	{
		std::cout << "NULL Matrix" << std::endl;
	}
}

void hzhu_mat::disp(int a)
{
	if (data != NULL)
	{
		std::cout << "  A " << data->size1 << " x " << data->size2 << " Matrix" << std::endl;
		std::cout.precision(a);
		for (int i = 0; i < data->size1; i++)
		{
			for (int j = 0; j < data->size2; j++)
				std::cout << gsl_matrix_get(data, i, j) << "\t";
			std::cout << std::endl;
		}
	}
	else
	{
		std::cout << "NULL Matrix" << std::endl;
	}
	std::cout.precision(6);
}

hzhu_mat::hzhu_mat(hzhu_mat &s)
{
	if (s.data != NULL)
	{
		data = gsl_matrix_alloc(s.data->size1, s.data->size2);
		gsl_matrix_memcpy(data, s.data);
	}
}

hzhu_mat::hzhu_mat(const hzhu_mat & s)
{
	if (s.data != NULL)
	{
		data = gsl_matrix_alloc(s.data->size1, s.data->size2);
		gsl_matrix_memcpy(data, s.data);
	}
}

hzhu_mat::hzhu_mat(const gsl_matrix *m)
{
	if (m == NULL)
		data = NULL;
	else
	{
		data = gsl_matrix_alloc(m->size1, m->size2);
		gsl_matrix_memcpy(data, m);
	}
}

// Operator ------------------------------------------------------------------------------------------------------------

void hzhu_mat::operator*=(hzhu_mat a)
{
	if (data == NULL || a.data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ void hzhu_mat::operator*=(hzhu_mat& a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}

	if (data->size1 != a.data->size1 || data->size2 != a.data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension Mismatch @ void hzhu_mat::operator*=(hzhu_mat& a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}

	for (int i = 0; i < data->size1; i++)
		for (int j = 0; j < data->size2; j++)
			data->data[j + i * data->size2] = data->data[j + i * data->size2] * a.data->data[j + i * data->size2];
}


void hzhu_mat::operator/=(hzhu_mat a)
{
	if (data == NULL || a.data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ void hzhu_mat::operator/=(hzhu_mat& a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}

	if (data->size1 != a.data->size1 || data->size2 != a.data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension Mismatch @ void hzhu_mat::operator/=(hzhu_mat& a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}

	for (int i = 0; i < data->size1; i++)
		for (int j = 0; j < data->size2; j++)
			data->data[j + i * data->size2] = data->data[j + i * data->size2] / a.data->data[j + i * data->size2];
}

void hzhu_mat::operator+=(hzhu_mat a)
{
	if (data == NULL || a.data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ void hzhu_mat::operator+=(hzhu_mat& a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}

	if (data->size1 != a.data->size1 || data->size2 != a.data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension Mismatch @ void hzhu_mat::operator+=(hzhu_mat& a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	int j;
	for (int i = 0; i < data->size1; i++)
		for (j = 0; j < data->size2; j++)
			data->data[j + i * data->size2] = data->data[j + i * data->size2] + a.data->data[j + i * data->size2];
}

void hzhu_mat::operator-=(hzhu_mat a)
{
	if (data == NULL || a.data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ void hzhu_mat::operator-=(hzhu_mat& a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}

	if (data->size1 != a.data->size1 || data->size2 != a.data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension Mismatch @ void hzhu_mat::operator-=(hzhu_mat& a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	int j;
	for (int i = 0; i < data->size1; i++)
		for (j = 0; j < data->size2; j++)
			data->data[j + i * data->size2] = data->data[j + i * data->size2] - a.data->data[j + i * data->size2];
}

void hzhu_mat::operator^=(hzhu_mat a)
{
	if (data == NULL || a.data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ void hzhu_mat::operator^=(hzhu_mat& a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}

	if (data->size1 != a.data->size1 || data->size2 != a.data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension Mismatch @ void hzhu_mat::operator^=(hzhu_mat& a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	int j;
	for (int i = 0; i < data->size1; i++)
		for (j = 0; j < data->size2; j++)
			data->data[j + i * data->size2] = pow(data->data[j + i * data->size2], a.data->data[j + i * data->size2]);
}

hzhu_mat hzhu_mat::operator*(hzhu_mat a)
{
	if (data == NULL || a.data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::operator*(hzhu_mat &a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	if (data->size1 != a.data->size1 || data->size2 != a.data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension Mismatch @ hzhu_mat hzhu_mat::operator*(hzhu_mat &a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(this[0]);
	re *= a;
	return re;
}

hzhu_mat hzhu_mat::operator/(hzhu_mat a)
{
	if (data == NULL || a.data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::operator/(hzhu_mat &a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	if (data->size1 != a.data->size1 || data->size2 != a.data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension Mismatch @ hzhu_mat hzhu_mat::operator/(hzhu_mat &a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	hzhu_mat re(this[0]);
	re /= a;
	return re;
}

hzhu_mat hzhu_mat::operator+(hzhu_mat a)
{
	if (data == NULL || a.data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::operator+(hzhu_mat &a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	if (data->size1 != a.data->size1 || data->size2 != a.data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension Mismatch @ hzhu_mat hzhu_mat::operator+(hzhu_mat &a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	hzhu_mat re(this[0]);
	re += a;
	return re;
}

hzhu_mat hzhu_mat::operator-(hzhu_mat a)
{
	if (data == NULL || a.data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::operator-(hzhu_mat &a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	if (data->size1 != a.data->size1 || data->size2 != a.data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension Mismatch @ hzhu_mat hzhu_mat::operator-(hzhu_mat &a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	hzhu_mat re(this[0]);
	re -= a;
	return re;
}

hzhu_mat hzhu_mat::operator^(hzhu_mat a)
{
	if (data == NULL || a.data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::operator^(hzhu_mat &a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	if (data->size1 != a.data->size1 || data->size2 != a.data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension Mismatch @ hzhu_mat hzhu_mat::operator^(hzhu_mat &a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	int j;
	hzhu_mat re(this[0]);
	for (int i = 0; i < data->size1; i++)
		for (j = 0; j < data->size2; j++)
			re.data->data[j + i * data->size2] = pow(data->data[j + i * data->size2], a.data->data[j + i * data->size2]);
	return re;
}

void hzhu_mat::operator^=(double a)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ void hzhu_mat::operator^=(double a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	int j;
	for (int i = 0; i < data->size1; i++)
		for (j = 0; j < data->size2; j++)
			data->data[j + i * data->size2] = pow(data->data[j + i * data->size2], a);
}

void hzhu_mat::operator*=(double a)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ void hzhu_mat::operator*=(double& a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	gsl_matrix_scale(data, a);
}

void hzhu_mat::operator/=(double a)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ void hzhu_mat::operator/=(double& a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	gsl_matrix_scale(data, 1.0/a);
}

void hzhu_mat::operator+=(double a)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ void hzhu_mat::operator+=(double& a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	gsl_matrix_add_constant(data, a);
}

void hzhu_mat::operator-=(double a)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ void hzhu_mat::operator-=(double& a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	gsl_matrix_add_constant(data, -a);
}

hzhu_mat hzhu_mat::operator*(double a)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat operator*(double a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(this[0]);
	re *= a;
	return re;
}

hzhu_mat hzhu_mat::operator/(double a)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat operator/(double a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(this[0]);
	re /= a;
	return re;
}

hzhu_mat hzhu_mat::operator+(double a)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat operator+(double a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(this[0]);
	re += a;
	return re;
}

hzhu_mat hzhu_mat::operator-(double a)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat operator-(double a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(this[0]);
	re -= a;
	return re;
}

hzhu_mat hzhu_mat::operator^(double a)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat operator^(double a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(this[0]);
	int j;
	for (int i = 0; i < data->size1; i++)
		for (j = 0; j < data->size2; j++)
			re.data->data[j + i * data->size2] = pow(data->data[j + i * data->size2], a);
	return re;
}

double hzhu_mat::operator[](int a)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ double hzhu_mat::operator[](int a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return GSL_NAN;
	}
	if (a == -1) return data->size1;
	if (a == -2) return data->size2;
	if (a<0 || a>(data->size1*data->size2 - 1))
	{
		if (HZHU_WARNING_DISP)	std::cout << "Index out of bound @ double hzhu_mat::operator[](int a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return GSL_NAN;
	}
	return data->data[a];
}

double hzhu_mat::operator[](hzhu_mat a)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ double hzhu_mat::operator[](int a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return GSL_NAN;
	}
	if (a.data == NULL || a.data->size1!= 1 || a.data->size2 != 2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Invalid input @ double hzhu_mat::operator[](int a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return GSL_NAN;
	}
	if (a.data->data[0] < 0 || a.data->data[0] >= data->size1 || a.data->data[1] < 0 || a.data->data[1] >= data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Index out of bound @ double hzhu_mat::operator[](int a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return GSL_NAN;
	}
	return gsl_matrix_get(data, a.data->data[0], a.data->data[1]);
}

// Matrix operation ------------------------------------------------------------------------------------------------------------.

hzhu_mat hzhu_mat::transpose()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::transpose()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(data->size2, data->size1);
	int j;
	for (int i = 0; i < data->size1; i++)
		for (j = 0; j < data->size2; j++)
			re.data->data[i + j * data->size1] = data->data[j + i * data->size2];
	return re;
}

void hzhu_mat::transpose_self()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::transpose_self()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}
	hzhu_mat tmp = this->transpose();
	gsl_matrix_free(data);
	data = gsl_matrix_alloc(tmp.data->size1, tmp.data->size2);
	gsl_matrix_memcpy(data, tmp.data);
}

hzhu_mat hzhu_mat::inv()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::inv()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	if (data->size1 != data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Non-square Matrix @ hzhu_mat hzhu_mat::inv()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	int n = (int)data->size1;
	gsl_permutation * p = gsl_permutation_alloc(n);
	int signum;

	hzhu_mat re(n);
	hzhu_mat tmp(this[0]);
	gsl_linalg_LU_decomp(tmp.data, p, &signum);

	if (gsl_linalg_LU_det(tmp.data, signum) == 0)
	{
		gsl_permutation_free(p);
		if (HZHU_WARNING_DISP)	std::cout << "Matrix Singular @ hzhu_mat hzhu_mat::inv()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	gsl_linalg_LU_invert(tmp.data, p, re.data);
	gsl_permutation_free(p);
	return re;
}

hzhu_mat hzhu_mat::pinv()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::pinv()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	if (data->size1 < data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "size1<size2 @ hzhu_mat hzhu_mat::pinv()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	int n = (int)data->size2;
	hzhu_mat U(this[0]);
	hzhu_mat V(n);
	gsl_vector *s = gsl_vector_alloc(n);
	gsl_vector *work = gsl_vector_alloc(n);
	gsl_linalg_SV_decomp(U.data, V.data, s, work);

	hzhu_mat S(0.0, n);
	for (int i = 0; i < n; i++) gsl_matrix_set(S.data, i, i, 1.0 / s->data[i]);

	gsl_vector_free(s);
	gsl_vector_free(work);

	hzhu_mat tmp = V.product(S);
	hzhu_mat tmp2 = U.transpose();
	return tmp.product(tmp2);
}

hzhu_mat hzhu_mat::product(hzhu_mat a)
{
	if (data == NULL || a.data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::product(hzhu_mat &a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	if (data->size2 != a.data->size1)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension Mismatch @ hzhu_mat hzhu_mat::product(hzhu_mat &a)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(data->size1, a.data->size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, data, a.data, 0, re.data);
	return re;
}

// IO ------------------------------------------------------------------------------------------------------------

std::string hzhu_mat::print_to_string(int precision = 8)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension Mismatch @ std::string hzhu_mat::print_to_string()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return std::string();
	}
	std::stringstream stream;
	int i, j;
	for (i = 0; i < data->size1; i++)
	{
		for (j = 0; j < data->size2; j++)
		{
			//re.append(std::to_string(data->data[j + i * data->size2]));
			
			stream << std::fixed << std::setprecision(precision) << data->data[j + i * data->size2];
			if (j != data->size2 - 1) stream << ",";
		}
		stream << "\n";
	}
	return stream.str();
}

int hzhu_mat::save_to_file(std::string name, int precision = 8)
{
	std::ofstream myfile;
	myfile.open(name, std::ios::out);
	if (myfile.is_open())
	{
		myfile<< this->print_to_string(precision);
		myfile.close();
		return 1;
	}
	else
	{
		if (HZHU_WARNING_DISP)	std::cout << "Cannot save file @ int hzhu_mat::save_to_file(std::string name)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return 0;
	}
}

int hzhu_mat::save_to_file(std::string name, std::string path, int precision = 8)
{
	std::string dir(path);
	dir.append("/");
	dir.append(name);
	return save_to_file(dir, precision);
}


int hzhu_mat::read_from_file(std::string name)
{
	std::ifstream myfile;
	myfile.open(name, std::ios::in);
	if (myfile.is_open())
	{
		int line_counter = 0;
		std::string line_string;
		while (std::getline(myfile, line_string))
		{
			line_counter++;
		}
		int i;

		myfile.clear();
		myfile.seekg(0, std::ios::beg);

		std::getline(myfile, line_string);
		int col_counter = 0;
		for (i = 0; i < line_string.length(); i++)
		{
			if (line_string.at(i) == ',') col_counter++;
		}

		//std::cout << line_string;
		//std::cout << "line_counter = " << line_counter << std::endl;
		//std::cout << "line_string.length() = " << line_string.length() << std::endl;

		if (line_string.length() <= 0 || line_counter <= 0)
		{
			if (HZHU_WARNING_DISP)	std::cout << "Empty file @ int hzhu_mat::read_from_file(std::string name)" << std::endl;
			if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
			return 0;
		}

		if (data == NULL)
		{
			data = gsl_matrix_alloc(line_counter, col_counter + 1);
		}
		else
		{
			gsl_matrix_free(data);
			data = gsl_matrix_alloc(line_counter, col_counter + 1);
		}

		myfile.clear();
		myfile.seekg(0, std::ios::beg);

		int j;
		for (int i = 0; i < line_counter; i++)
		{
			std::string tmp;
			for (j = 0; j < col_counter; j++)
			{
				std::getline(myfile, tmp, ',');
				//std::cout << i << "\t" << j << "\t" << tmp << std::endl;
				gsl_matrix_set(data, i, j, std::stod(tmp));
			}
			std::getline(myfile, tmp, '\n');
			//std::cout << i << "\t" << j << "\t" << tmp << std::endl;
			gsl_matrix_set(data, i, j, std::stod(tmp));
		}

		myfile.close();
		return 1;
	}
	else
	{
		if (HZHU_WARNING_DISP)	std::cout << "Cannot open file @ int hzhu_mat::read_from_file(std::string name)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return 0;
	}
}

int hzhu_mat::read_from_file(std::string name, std::string path)
{
	std::string dir(path);
	dir.append("/");
	dir.append(name);
	return read_from_file(dir);
}

// Statistics ------------------------------------------------------------------------------------------------------------

hzhu_mat hzhu_mat::quantile_row(double q)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::quantile_row(double q)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	int N1 = data->size1;
	int N2 = data->size2;
	gsl_vector *v = gsl_vector_alloc(N2);

	hzhu_mat re(N1, 1);

	int i, j;
	for (i = 0; i < N1; i++)
	{
		for (j = 0; j < N2; j++)
		{
			v->data[j] = data->data[i*N2 + j];
		}
		gsl_sort_vector(v);
		re.data->data[i] =  gsl_stats_quantile_from_sorted_data(v->data, 1, N2, q);
	}

	gsl_vector_free(v);
	return re;
}

hzhu_mat hzhu_mat::quantile_col(double q)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::quantile_col(double q)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	int N1 = data->size1;
	int N2 = data->size2;
	gsl_vector *v = gsl_vector_alloc(N1);

	hzhu_mat re(1,N2);

	int i, j;
	for (i = 0; i < N2; i++)
	{
		for (j = 0; j < N1; j++)
		{
			v->data[j] = data->data[i + j * N2];
		}
		gsl_sort_vector(v);
		re.data->data[i] = gsl_stats_quantile_from_sorted_data(v->data, 1, N1, q);
	}

	gsl_vector_free(v);
	return re;
}

hzhu_mat hzhu_mat::quantile_all(double q)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::quantile_all(double q)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	hzhu_mat re(1);

	int N = data->size1*data->size2;
	hzhu_mat tmp(*this);
	gsl_sort(tmp.data->data, 1, N);

	re.data->data[0] = gsl_stats_quantile_from_sorted_data(tmp.data->data, 1, N, q);

	return re;
}

hzhu_mat hzhu_mat::median_row()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::median_row()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	return quantile_row(0.5);
}

hzhu_mat hzhu_mat::median_col()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::median_col()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	return quantile_col(0.5);
}

hzhu_mat hzhu_mat::median_all()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::median_all()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	return quantile_all(0.5);
}

hzhu_mat hzhu_mat::max_index_value_row()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::max_index_value_row()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re((int)data->size1, 3);
	int i;
	gsl_vector *v = gsl_vector_alloc(data->size2);
	for (i = 0; i < data->size1; i++)
	{
		gsl_matrix_get_row(v, data, i);
		size_t index = gsl_vector_max_index(v);
		gsl_matrix_set(re.data, i, 0, i);
		gsl_matrix_set(re.data, i, 1, index);
		gsl_matrix_set(re.data, i, 2, v->data[index]);
	}
	gsl_vector_free(v);
	return re;
}


hzhu_mat hzhu_mat::max_index_row()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::max_index_row()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re((int)data->size1, 2);
	int i;
	gsl_vector *v = gsl_vector_alloc(data->size2);
	for (i = 0; i < data->size1; i++)
	{
		gsl_matrix_get_row(v, data, i);
		size_t index = gsl_vector_max_index(v);
		gsl_matrix_set(re.data, i, 0, i);
		gsl_matrix_set(re.data, i, 1, index);
	}
	gsl_vector_free(v);
	return re;
}

hzhu_mat hzhu_mat::max_value_row()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::max_value_row()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re((int)data->size1, 1);
	int i;
	gsl_vector *v = gsl_vector_alloc(data->size2);
	for (i = 0; i < data->size1; i++)
	{
		gsl_matrix_get_row(v, data, i);
		size_t index = gsl_vector_max_index(v);
		gsl_matrix_set(re.data, i, 0, v->data[index]);
	}
	gsl_vector_free(v);
	return re;
}


hzhu_mat hzhu_mat::max_index_value_col()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::max_index_value_col()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(3, (int)data->size2);
	int i;
	gsl_vector *v = gsl_vector_alloc(data->size1);
	for (i = 0; i < data->size2; i++)
	{
		gsl_matrix_get_col(v, data, i);
		size_t index = gsl_vector_max_index(v);
		gsl_matrix_set(re.data, 1, i, i);
		gsl_matrix_set(re.data, 0, i, index);
		gsl_matrix_set(re.data, 2, i, v->data[index]);
	}
	gsl_vector_free(v);
	return re;
}


hzhu_mat hzhu_mat::max_index_col()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::max_index_col()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(2, (int)data->size2);
	int i;
	gsl_vector *v = gsl_vector_alloc(data->size1);
	for (i = 0; i < data->size2; i++)
	{
		gsl_matrix_get_col(v, data, i);
		size_t index = gsl_vector_max_index(v);
		gsl_matrix_set(re.data, 1, i, i);
		gsl_matrix_set(re.data, 0, i, index);
	}
	gsl_vector_free(v);
	return re;
}

hzhu_mat hzhu_mat::max_value_col()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::max_value_col()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(1, (int)data->size2);
	int i;
	gsl_vector *v = gsl_vector_alloc(data->size1);
	for (i = 0; i < data->size2; i++)
	{
		gsl_matrix_get_col(v, data, i);
		size_t index = gsl_vector_max_index(v);
		gsl_matrix_set(re.data, 0, i, v->data[index]);
	}
	gsl_vector_free(v);
	return re;
}


hzhu_mat hzhu_mat::max_index_value_all()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::max_value_index_all()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(1, 4);
	size_t index = gsl_stats_max_index(data->data, 1, data->size1*data->size2);
	gsl_matrix_set(re.data, 0, 3, index);
	gsl_matrix_set(re.data, 0, 2, data->data[index]);
	gsl_matrix_set(re.data, 0, 1, index%data->size2);
	gsl_matrix_set(re.data, 0, 0, (int)index / (int)data->size2);

	return re;
}

hzhu_mat hzhu_mat::max_value_all()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::max_value_all()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(1, 1);
	size_t index = gsl_stats_max_index(data->data, 1, data->size1*data->size2);
	gsl_matrix_set(re.data, 0, 0, data->data[index]);

	return re;
}

hzhu_mat hzhu_mat::max_index_all()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::max_index_all()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(1, 2);
	size_t index = gsl_stats_max_index(data->data, 1, data->size1*data->size2);
	gsl_matrix_set(re.data, 0, 1, index%data->size2);
	gsl_matrix_set(re.data, 0, 0, (int)index / (int)data->size2);

	return re;
}

hzhu_mat hzhu_mat::min_index_value_row()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::min_index_value_row()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re((int)data->size1, 3);
	int i;
	gsl_vector *v = gsl_vector_alloc(data->size2);
	for (i = 0; i < data->size1; i++)
	{
		gsl_matrix_get_row(v, data, i);
		size_t index = gsl_vector_min_index(v);
		gsl_matrix_set(re.data, i, 0, i);
		gsl_matrix_set(re.data, i, 1, index);
		gsl_matrix_set(re.data, i, 2, v->data[index]);
	}
	gsl_vector_free(v);
	return re;
}


hzhu_mat hzhu_mat::min_index_row()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::min_index_row()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re((int)data->size1, 2);
	int i;
	gsl_vector *v = gsl_vector_alloc(data->size2);
	for (i = 0; i < data->size1; i++)
	{
		gsl_matrix_get_row(v, data, i);
		size_t index = gsl_vector_min_index(v);
		gsl_matrix_set(re.data, i, 0, i);
		gsl_matrix_set(re.data, i, 1, index);
	}
	gsl_vector_free(v);
	return re;
}

hzhu_mat hzhu_mat::min_value_row()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::min_value_row()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re((int)data->size1, 1);
	int i;
	gsl_vector *v = gsl_vector_alloc(data->size2);
	for (i = 0; i < data->size1; i++)
	{
		gsl_matrix_get_row(v, data, i);
		size_t index = gsl_vector_min_index(v);
		gsl_matrix_set(re.data, i, 0, v->data[index]);
	}
	gsl_vector_free(v);
	return re;
}


hzhu_mat hzhu_mat::min_index_value_col()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::min_index_value_col()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(3, (int)data->size2);
	int i;
	gsl_vector *v = gsl_vector_alloc(data->size1);
	for (i = 0; i < data->size2; i++)
	{
		gsl_matrix_get_col(v, data, i);
		size_t index = gsl_vector_min_index(v);
		gsl_matrix_set(re.data, 1, i, i);
		gsl_matrix_set(re.data, 0, i, index);
		gsl_matrix_set(re.data, 2, i, v->data[index]);
	}
	gsl_vector_free(v);
	return re;
}


hzhu_mat hzhu_mat::min_index_col()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::min_index_col()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(2, (int)data->size2);
	int i;
	gsl_vector *v = gsl_vector_alloc(data->size1);
	for (i = 0; i < data->size2; i++)
	{
		gsl_matrix_get_col(v, data, i);
		size_t index = gsl_vector_min_index(v);
		gsl_matrix_set(re.data, 1, i, i);
		gsl_matrix_set(re.data, 0, i, index);
	}
	gsl_vector_free(v);
	return re;
}

hzhu_mat hzhu_mat::min_value_col()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::min_value_col()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(1, (int)data->size2);
	int i;
	gsl_vector *v = gsl_vector_alloc(data->size1);
	for (i = 0; i < data->size2; i++)
	{
		gsl_matrix_get_col(v, data, i);
		size_t index = gsl_vector_min_index(v);
		gsl_matrix_set(re.data, 0, i, v->data[index]);
	}
	gsl_vector_free(v);
	return re;
}


hzhu_mat hzhu_mat::min_index_value_all()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::min_value_index_all()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(1, 4);
	size_t index = gsl_stats_min_index(data->data, 1, data->size1*data->size2);
	gsl_matrix_set(re.data, 0, 3, index);
	gsl_matrix_set(re.data, 0, 2, data->data[index]);
	gsl_matrix_set(re.data, 0, 1, index%data->size2);
	gsl_matrix_set(re.data, 0, 0, (int)index / (int)data->size2);

	return re;
}

hzhu_mat hzhu_mat::min_value_all()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::min_value_all()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(1, 1);
	size_t index = gsl_stats_min_index(data->data, 1, data->size1*data->size2);
	gsl_matrix_set(re.data, 0, 0, data->data[index]);

	return re;
}

hzhu_mat hzhu_mat::min_index_all()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::min_index_all()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(1, 2);
	size_t index = gsl_stats_min_index(data->data, 1, data->size1*data->size2);
	gsl_matrix_set(re.data, 0, 1, index%data->size2);
	gsl_matrix_set(re.data, 0, 0, (int)index / (int)data->size2);

	return re;
}

hzhu_mat hzhu_mat::sum_row()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::sum_row()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re((int)data->size1, 1);
	int i, j;
	for (i = 0; i < data->size1; i++)
	{
		double tmp = 0.0;
		for (j = 0; j < data->size2; j++)
		{
			tmp += data->data[j + i * data->size2];
		}
		re.data->data[i] = tmp;
	}

	return re;
}

hzhu_mat hzhu_mat::sum_col()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::sum_col()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(1, (int)data->size2);
	int i, j;
	for (j = 0; j < data->size2; j++)
	{
		double tmp = 0.0;
		for (i = 0; i < data->size1; i++)
		{
			tmp += data->data[j + i * data->size2];
		}
		re.data->data[j] = tmp;
	}

	return re;
}

hzhu_mat hzhu_mat::sum_all()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::sum_all()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	hzhu_mat re(1, 1);
	double tmp = 0.0;
	for (int i = 0; i < data->size1*data->size2; i++)
		tmp += data->data[i];
	re.data->data[0] = tmp;
	return re;
}



hzhu_mat hzhu_mat::mean_row()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::mean_row()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	
	return this->sum_row() / (double)data->size2;
}

hzhu_mat hzhu_mat::mean_col()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::mean_col()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	return this->sum_col() / (double)data->size1;
}

hzhu_mat hzhu_mat::mean_all()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::mean_col()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	return this->sum_all() / (double)data->size1 / (double)data->size2;
}


hzhu_mat hzhu_mat::var_row()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::var_row()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	
	hzhu_mat M1 = this->mean_row();
	hzhu_mat M2 = (this[0] * this[0]).mean_row();
	M1 *= M1;

	return M2 - M1;
}

hzhu_mat hzhu_mat::var_col()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::var_col()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	hzhu_mat M1 = this->mean_col();
	hzhu_mat M2 = (this[0] * this[0]).mean_col();
	M1 *= M1;

	return M2 - M1;
}

hzhu_mat hzhu_mat::var_all()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::var_all()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	hzhu_mat M1 = this->mean_all();
	hzhu_mat M2 = (this[0]*this[0]).mean_all();
	M1 *= M1;

	return M2 - M1;
}

hzhu_mat hzhu_mat::uvar_row()
{
	double s = data->size2;
	return this->var_row()*s / (s - 1.0);
}

hzhu_mat hzhu_mat::uvar_col()
{
	double s = data->size1;
	return this->var_col()*s / (s - 1.0);
}

hzhu_mat hzhu_mat::uvar_all()
{
	double s = data->size2*data->size1;
	return this->var_all()*s / (s - 1.0);
}

hzhu_mat hzhu_mat::get_row(int x)
{
	if (data == NULL)
	{ 
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::get_row(int x)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	if (x < 0 || x >= data->size1)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Invalid input @ hzhu_mat hzhu_mat::get_row(int x)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	return get_sub(x, 0, x, data->size2 - 1);
}

hzhu_mat hzhu_mat::get_col(int x)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::get_col(int x)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	if (x < 0 || x >= data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Invalid input @ hzhu_mat hzhu_mat::get_col(int x)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	return get_sub(0, x, data->size1 - 1, x);
}


hzhu_mat hzhu_mat::get_sub(int x1, int y1, int x2, int y2)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ hzhu_mat hzhu_mat::get_sub(int x1, int y1, int x2, int y2)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	if (x1 > x2 || y1 > y2 || x1<0 || y1<0 || x2<0 || y2<0 || x1>=data->size1 || x2>=data->size1 || y1>=data->size2 || y2>=data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Invalid input @ hzhu_mat hzhu_mat::get_sub(int x1, int y1, int x2, int y2)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	int x = x2 - x1 + 1;
	int y = y2 - y1 + 1;
	hzhu_mat re(x, y);

	for (int i = x1; i <= x2; i++)
		for (int j = y1; j <= y2; j++)
			re.data->data[(i - x1)*y + (j - y1)] = data->data[i*data->size2 + j];
	return re;
}

double hzhu_mat::tr()
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ double hzhu_mat::tr()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return GSL_NAN;
	}
	if (data->size1 != data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Non-square matrix @ double hzhu_mat::tr()" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return GSL_NAN;
	}
	double re = 0.0;
	for (int i = 0; i < data->size1; i++)
		re += data->data[i + i * data->size1];
	return re;
}

void hzhu_mat::append_left(double a)
{
	if (data == NULL)
	{
		data = gsl_matrix_alloc(1, 1);
		data->data[0] = a;
	}
	else
	{
		gsl_matrix *re = gsl_matrix_alloc(data->size1, data->size2 + 1);
		int i, j;
		for (i = 0; i < data->size1; i++) re->data[re->size2*i] = a;
		for (i = 0; i < data->size1; i++)
		{
			for(j = 0; j < data->size2; j++)
			{
				re->data[j + 1 + i * re->size2] = data->data[j + i * data->size2];
			}
		}
		gsl_matrix_free(data);
		data = re;
	}
}

void hzhu_mat::append_right(double a)
{
	if (data == NULL)
	{
		data = gsl_matrix_alloc(1, 1);
		data->data[0] = a;
	}
	else
	{
		gsl_matrix *re = gsl_matrix_alloc(data->size1, data->size2 + 1);
		int i, j;
		for (i = 0; i < re->size1; i++) re->data[re->size2*i + re->size2 - 1] = a;
		for (i = 0; i < data->size1; i++)
		{
			for (j = 0; j < data->size2; j++)
			{
				re->data[j + i * re->size2] = data->data[j + i * data->size2];
			}
		}
		gsl_matrix_free(data);
		data = re;
	}
}

void hzhu_mat::append_top(double a)
{
	if (data == NULL)
	{
		data = gsl_matrix_alloc(1, 1);
		data->data[0] = a;
	}
	else
	{
		gsl_matrix *re = gsl_matrix_alloc(data->size1 + 1, data->size2);
		int i, j;
		for (i = 0; i < data->size2; i++) re->data[i] = a;
		for (i = 0; i < data->size1; i++)
		{
			for (j = 0; j < data->size2; j++)
			{
				re->data[(i + 1)*data->size2 + j] = data->data[i*data->size2 + j];
			}
		}
		gsl_matrix_free(data);
		data = re;
	}
}

void hzhu_mat::append_bottom(double a)
{
	if (data == NULL)
	{
		data = gsl_matrix_alloc(1, 1);
		data->data[0] = a;
	}
	else
	{
		gsl_matrix *re = gsl_matrix_alloc(data->size1 + 1, data->size2);
		int i, j;
		for (i = 0; i < data->size2; i++) re->data[i + data->size1*data->size2] = a;
		for (i = 0; i < data->size1; i++)
		{
			for (j = 0; j < data->size2; j++)
			{
				re->data[i*data->size2 + j] = data->data[i*data->size2 + j];
			}
		}
		gsl_matrix_free(data);
		data = re;
	}
}

void hzhu_mat::append_left(hzhu_mat& m)
{
	if (data == NULL)
	{
		if (m.data == NULL)
		{
			if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ void hzhu_mat::append_left(hzhu_mat& m)" << std::endl;
			if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		}
		else
		{
			data = gsl_matrix_alloc(m.data->size1, m.data->size2);
			gsl_matrix_memcpy(data, m.data);
		}
		return;
	}
	if (data->size1 != m.data->size1)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension mismatch @ void hzhu_mat::append_left(hzhu_mat& m)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}

	int i, j;
	gsl_matrix *re = gsl_matrix_alloc(data->size1, data->size2 + m.data->size2);
	for(i = 0;i<m.data->size1;i++)
		for (j = 0; j < m.data->size2; j++)
		{
			re->data[j + i * re->size2] = m.data->data[j + i * m.data->size2];
		}
	for (i = 0; i < data->size1; i++)
		for (j = 0; j < data->size2; j++)
			re->data[j + m.data->size2 + i * re->size2] = data->data[j + i * data->size2];
	gsl_matrix_free(data);
	data = re;
}

void hzhu_mat::append_right(hzhu_mat& m)
{
	if (data == NULL)
	{
		if (m.data == NULL)
		{
			if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ void hzhu_mat::append_left(hzhu_mat& m)" << std::endl;
			if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		}
		else
		{
			data = gsl_matrix_alloc(m.data->size1, m.data->size2);
			gsl_matrix_memcpy(data, m.data);
		}
		return;
	}
	if (data->size1 != m.data->size1)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension mismatch @ void hzhu_mat::append_left(hzhu_mat& m)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}

	int i, j;
	gsl_matrix *re = gsl_matrix_alloc(data->size1, data->size2 + m.data->size2);
	for (i = 0; i < m.data->size1; i++)
		for (j = 0; j < m.data->size2; j++)
		{
			re->data[j + data->size2 + i * re->size2] = m.data->data[j + i * m.data->size2];
		}
	for (i = 0; i < data->size1; i++)
		for (j = 0; j < data->size2; j++)
			re->data[j + i * re->size2] = data->data[j + i * data->size2];

	gsl_matrix_free(data);
	data = re;
}


void hzhu_mat::append_top(hzhu_mat& m)
{
	if (data == NULL)
	{
		if (m.data == NULL)
		{
			if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ void hzhu_mat::append_top(hzhu_mat& m)" << std::endl;
			if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		}
		else
		{
			data = gsl_matrix_alloc(m.data->size1, m.data->size2);
			gsl_matrix_memcpy(data, m.data);
		}
		return;
	}
	if (data->size2 != m.data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension mismatch @ void hzhu_mat::append_top(hzhu_mat& m)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}

	gsl_matrix *re = gsl_matrix_alloc(m.data->size1 + data->size1, data->size2);
	int i, j;
	for (i = 0; i < data->size1; i++)
		for (j = 0; j < data->size2; j++)
			re->data[(i + m.data->size1)*re->size2 + j] = data->data[j + i * data->size2];

	for (i = 0; i < m.data->size1; i++)
		for (j = 0; j < m.data->size2; j++)
			re->data[i*re->size2 + j] = m.data->data[j + i * m.data->size2];

	gsl_matrix_free(data);
	data = re;
}

void hzhu_mat::append_bottom(hzhu_mat& m)
{
	if (data == NULL)
	{
		if (m.data == NULL)
		{
			if (HZHU_WARNING_DISP)	std::cout << "NULL Matrix @ void hzhu_mat::append_top(hzhu_mat& m)" << std::endl;
			if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		}
		else
		{
			data = gsl_matrix_alloc(m.data->size1, m.data->size2);
			gsl_matrix_memcpy(data, m.data);
		}
		return;
	}
	if (data->size2 != m.data->size2)
	{
		if (HZHU_WARNING_DISP)	std::cout << "Dimension mismatch @ void hzhu_mat::append_top(hzhu_mat& m)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return;
	}

	gsl_matrix *re = gsl_matrix_alloc(m.data->size1 + data->size1, data->size2);
	int i, j;
	for (i = 0; i < data->size1; i++)
		for (j = 0; j < data->size2; j++)
			re->data[i *re->size2 + j] = data->data[j + i * data->size2];

	for (i = 0; i < m.data->size1; i++)
		for (j = 0; j < m.data->size2; j++)
			re->data[(i + data->size1)*re->size2 + j] = m.data->data[j + i * m.data->size2];

	gsl_matrix_free(data);
	data = re;
}

hzhu_mat::hzhu_mat(gsl_rng *r)// one standard normal variable
{
	data = gsl_matrix_alloc(1, 1);
	data->data[0] = gsl_ran_gaussian_ziggurat(r, 1.0);
}

hzhu_mat::hzhu_mat(gsl_rng *r, int size1, int size2)// A matrix of standard normal variables
{
	int m = abs(size1);
	int n = abs(size2);
	data = gsl_matrix_alloc(m, n);
	int N = m * n;
	for (int i = 0; i < N; i++)
		data->data[i] = gsl_ran_gaussian_ziggurat(r, 1.0);
}

hzhu_mat::hzhu_mat(gsl_rng *r, int size1, int size2, double mean, double sd)// A matrix of normal variables N(mean,sd^2)
{
	int m = abs(size1);
	int n = abs(size2);
	data = gsl_matrix_alloc(m, n);
	int N = m * n;
	for (int i = 0; i < N; i++)
		data->data[i] = gsl_ran_gaussian_ziggurat(r, sd) + mean;
}

// A matrix of type variables, supported type include: uniform, gaussian, log normal, exponential, Laplace, Cauchy, Rayleigh, Gamma 
hzhu_mat::hzhu_mat(gsl_rng *r, std::string type, int size1, int size2, double p1, double p2 = 0)
{
	int m = abs(size1);
	int n = abs(size2);
	data = gsl_matrix_alloc(m, n);
	int N = m * n;

	std::string Dis[10];
	int i = 0;
	Dis[i++].append("Gaussian,gaussian,normal,Normal");//0
	Dis[i++].append("flat,uniform,Flat,Uniform");//1
	Dis[i++].append("Gamma,gamma");//2
	Dis[i++].append("exponential,Exponential");//3
	Dis[i++].append("Laplace,laplace");//4
	Dis[i++].append("Rayleigh,rayleigh");//5
	Dis[i++].append("log-normal,Log-normal,log normal,Log normal,Log Normal,log Normal,lognormal,Lognormal");//6
	Dis[i++].append("Cauchy,cauchy");//7

	int value = -1;
	for (i = 0; i < 10; i++)
	{
		if (Dis[i].find(type) != std::string::npos)
		{
			value = i;
			break;
		}
	}

	switch (value)
	{
	case -1:
		if (HZHU_WARNING_DISP) std::cout << "Invalid type @ hzhu_mat::hzhu_mat(gsl_rng *r, std::string type, int size1, int size2, double p1, double p2)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		gsl_matrix_free(data);
		data = NULL;
		return;
	case 0:
		for (i = 0; i < N; i++) data->data[i] = gsl_ran_gaussian_ziggurat(r, p1) + p2;
		return;
	case 1:
		for (i = 0; i < N; i++) data->data[i] = gsl_ran_flat(r, p1, p2);
		return;
	case 2:
		for (i = 0; i < N; i++) data->data[i] = gsl_ran_gamma(r, p1, p2);
		return;
	case 3:
		for (i = 0; i < N; i++) data->data[i] = gsl_ran_exponential(r, p1);
		return;
	case 4:
		for (i = 0; i < N; i++) data->data[i] = gsl_ran_laplace(r, p1) + p2;
		return;
	case 5:
		for (i = 0; i < N; i++) data->data[i] = gsl_ran_rayleigh(r, p1);
		return;
	case 6:
		for (i = 0; i < N; i++) data->data[i] = gsl_ran_lognormal(r, p1, p2);
		return;
	case 7:
		for (i = 0; i < N; i++) data->data[i] = gsl_ran_cauchy(r, p1);
		return;
	}
}

hzhu_mat hzhu_mat::diff(int n)
{
	if (data == NULL)
	{
		if (HZHU_WARNING_DISP) std::cout << "Null matrix @ hzhu_mat hzhu_mat::diff(int n)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}
	int n1 = data->size1;
	int n2 = data->size2;

	if (n1 <= n || n <= 0)
	{
		if (HZHU_WARNING_DISP) std::cout << " Invalid input @ hzhu_mat hzhu_mat::diff(int n)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	hzhu_mat re(n1 - n, n2);
	int j;
	for (int i = 0; i < n2; i++)
	{
		for (int j = 0; j < n1 - n; j++)
			re.data->data[j*n2 + i] = data->data[(j + n)*n2 + i] - data->data[j*n2 + i];
	}

	return re;
}

hzhu_mat hzhu_mat::diff()
{
	return diff(1);
}

// Other Function Realization ============================================================================================================
hzhu_mat hzhu_line_space(double a, double b, int n)
{
	int N = abs(n);
	hzhu_mat re(N, 1);
	double gap = (b - a) / (N - 1);
	for (int i = 0; i < N; i++)
		re.data->data[i] = a + gap * i;
	return re;
}

hzhu_mat hzhu_equal_space(double a, double b, double c)
{
	if ((c - a)*b > 0)
	{
		int N = (int)((c - a) / b + 1);
		hzhu_mat re(N, 1);
		for (int i = 0; i < N; i++)
			re.data->data[i] = a + i * b;
		return re;
	}
	else
	{
		if (HZHU_WARNING_DISP) std::cout << "NULL Matrix returned @ hzhu_mat hzhu_equal_space(double a, double b, double c)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

}

hzhu_mat hzhu_mat_index(int x, int y)
{
	hzhu_mat re(1, 2);
	re.data->data[0] = x;
	re.data->data[1] = y;
	return re;
}

hzhu_mat hzhu_mat_repmat(hzhu_mat m, int a, int b)
{
	if (m.data == NULL)
	{
		if (HZHU_WARNING_DISP) std::cout << "NULL Matrix @ hzhu_mat hzhu_mat_repmat(hzhu_mat &m, int a, int b)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	if (a <= 0 || b <= 0)
	{
		if (HZHU_WARNING_DISP) std::cout << "Invalid input @ hzhu_mat hzhu_mat_repmat(hzhu_mat &m, int a, int b)" << std::endl;
		if (HZHU_WARNING_EXIT) exit(EXIT_FAILURE);
		return hzhu_mat();
	}

	int x = m.data->size1*a;
	int y = m.data->size2*b;

	hzhu_mat re(x, y);
	int i, j;
	int p, q;
	for (i = 0; i < a; i++)
	{
		for (j = 0; j < b; j++)
		{
			for (p = 0; p < m.data->size1; p++)
				for (q = 0; q < m.data->size2; q++)
					re.data->data[(i*m.data->size1 + p)*y + q + j * m.data->size2] = m.data->data[q + p * m.data->size2];
		}
	}
	return re;
}

#endif
