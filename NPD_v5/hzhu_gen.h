#pragma once

#ifndef hzhu_gen_h
#define hzhu_gen_h

#include "hzhu_mat.h"
#include <gsl/gsl_rng.h>
#include <omp.h>
#include <string>
#include <iostream>
#include <ctime>

class hzhu_gen
{
public:
	hzhu_gen();
	hzhu_gen(int);
	~hzhu_gen();

	std::time_t Time;
	std::tm* Now;
	std::string ID;
	double start_time;
	double counter;
	bool enable_disp = true;

	int max_possible_threads;
	gsl_rng **r;
	unsigned long int *seeds;

public:
	void out();
	void out(double);
	void out(std::string s, double);
	void out(const char *, double);
	void out(std::string s);
	void out(const char *);
	void out(int);
	void out(double, double);
	void out(double, double, double);
	std::string get_time_string();
	void disp_runtime();
	void disp_counter_time();
	void start_counter();

	void operator<<(double);
	void operator<<(int);
	void operator<<(const char*);
	void operator<<(std::string);

	void set_enable_disp(bool);
};

// General function ===================================================================================================
std::string hzhu_gen_second_to_string(double);
void hzhu_gen_free(hzhu_mat *);
void hzhu_gen_free(double *);
void hzhu_gen_free(int *);
void hzhu_gen_free(gsl_vector *);
void hzhu_gen_free(gsl_matrix *);

std::string hzhu_gen_int_to_string(int);
std::string hzhu_gen_int_to_string(int, int);

// Class function realization =========================================================================================
hzhu_gen::hzhu_gen()
{
	start_time = omp_get_wtime();
	ID.append(this->get_time_string());
	max_possible_threads = omp_get_max_threads();
	r = new gsl_rng*[max_possible_threads];
	seeds = new unsigned long int[max_possible_threads];
	Time = std::time(0);
	Now = std::localtime(&Time);
	for (int i = 0; i < max_possible_threads; i++)
	{
		r[i] = gsl_rng_alloc(gsl_rng_taus2);
		seeds[i] = (Now->tm_sec + 1)*(Now->tm_min + 1)*(Now->tm_hour + 1)*(Now->tm_mday)*(Now->tm_mon + 1) - 1 + 17 * i;
		gsl_rng_set(r[i], seeds[i]);
	}
}

hzhu_gen::hzhu_gen(int n)
{
	start_time = omp_get_wtime();
	ID.append(this->get_time_string());
	max_possible_threads = omp_get_max_threads();
	r = new gsl_rng*[max_possible_threads];
	seeds = new unsigned long int[max_possible_threads];
	Time = std::time(0);
	Now = std::localtime(&Time);
	for (int i = 0; i < max_possible_threads; i++)
	{
		r[i] = gsl_rng_alloc(gsl_rng_taus2);
		seeds[i] = (Now->tm_sec + 1)*(Now->tm_min + 1)*(Now->tm_hour + 1)*(Now->tm_mday)*(Now->tm_mon + 1) - 1 + 17 * i;
		gsl_rng_set(r[i], seeds[i]);
	}
	
	if (n > 0)
	{
		int N = n > max_possible_threads ? max_possible_threads : n;
		omp_set_num_threads(N);
	}
	else
	{
		if (-n < max_possible_threads) omp_set_num_threads(max_possible_threads + n);
	}
}

hzhu_gen::~hzhu_gen()
{
	for (int i = 0; i < max_possible_threads; i++)
	{
		gsl_rng_free(r[i]);
	}
	delete r;
	delete seeds;
}

std::string hzhu_gen::get_time_string()
{
	Time = std::time(0);
	Now = std::localtime(&Time);
	std::string re;
	re.append(std::to_string(Now->tm_year + 1900));
	re.append("_");
	if (Now->tm_mon + 1 < 10) re.append("0");
	re.append(std::to_string(Now->tm_mon + 1));
	re.append("_");
	if (Now->tm_mday < 10) re.append("0");
	re.append(std::to_string(Now->tm_mday));
	re.append("_");
	if (Now->tm_hour < 10) re.append("0");
	re.append(std::to_string(Now->tm_hour));
	re.append("_");
	if (Now->tm_min < 10) re.append("0");
	re.append(std::to_string(Now->tm_min));
	re.append("_");
	if (Now->tm_sec < 10) re.append("0");
	re.append(std::to_string(Now->tm_sec));

	return re;
}

void hzhu_gen::out()
{
	if (enable_disp == false) return;
	std::cout << std::endl;
}

void hzhu_gen::out(double a)
{
	if (enable_disp == false) return;
	std::cout << a << std::endl;
}

void hzhu_gen::out(std::string s, double a)
{
	if (enable_disp == false) return;
	std::cout << s.c_str() << "\t" << a << std::endl;
}

void hzhu_gen::out(const char *s, double a)
{
	if (enable_disp == false) return;
	std::cout << s << "\t" << a << std::endl;
}

void hzhu_gen::out(int a)
{
	if (enable_disp == false) return;
	std::cout << a << std::endl;
}

void hzhu_gen::out(std::string s)
{
	if (enable_disp == false) return;
	std::cout << s.c_str() << std::endl;
}

void hzhu_gen::out(const char *s)
{
	if (enable_disp == false) return;
	std::cout << s << std::endl;
}

void hzhu_gen::out(double a, double b)
{
	if (enable_disp == false) return;
	std::cout << a << "\t" << b << std::endl;
}
void hzhu_gen::out(double a, double b, double c)
{
	if (enable_disp == false) return;
	std::cout << a << "\t" << b << "\t" << c << std::endl;
}

void hzhu_gen::disp_runtime()
{
	if (enable_disp == false) return;
	std::cout << "Running time = " << hzhu_gen_second_to_string(omp_get_wtime() - start_time) << std::endl;
}

void hzhu_gen::disp_counter_time()
{
	if (enable_disp == false) return;
	std::cout << "Counter time = " << hzhu_gen_second_to_string(omp_get_wtime() - counter) << std::endl;
}

void hzhu_gen::start_counter()
{
	counter = omp_get_wtime();
}

void hzhu_gen::operator<<(double a)
{
	if (enable_disp == false) return;
	std::cout << a;
}

void hzhu_gen::operator<<(int a) 
{
	if (enable_disp == false) return;
	std::cout << a;
}

void hzhu_gen::operator<<(const char* a)
{
	if (enable_disp == false) return;
	std::cout << a;
}

void hzhu_gen::operator<<(std::string a)
{
	if (enable_disp == false) return;
	std::cout << a;
}

void hzhu_gen::set_enable_disp(bool f)
{
	if (enable_disp == false) return;
	enable_disp = f;
}

// General function realization ===========================================================================
std::string hzhu_gen_second_to_string(double t)
{
	double sec;
	int min, hour, day;
	day = t / 86400.0;
	hour = (t - day * 86400.0) / 3600.0;
	min = (t - day * 86400.0 - hour * 3600.0) / 60.0;
	sec = t - day * 86400.0 - hour * 3600.0 - min * 60.0;
	std::string re;
	re.append(std::to_string(day));
	re.append(" day, ");
	re.append(std::to_string(hour));
	re.append(" hour, ");
	re.append(std::to_string(min));
	re.append(" min, ");
	re.append(std::to_string(sec));
	re.append(" sec");
	return re;
}

void hzhu_gen_free(hzhu_mat *s)
{
	if (s != NULL) delete s;
}

void hzhu_gen_free(double *s)
{
	if (s != NULL) delete s;
}

void hzhu_gen_free(int *s)
{
	if (s != NULL) delete s;
}

void hzhu_gen_free(gsl_vector *s)
{
	if (s != NULL) gsl_vector_free(s);
}

void hzhu_gen_free(gsl_matrix *s)
{
	if (s != NULL) gsl_matrix_free(s);
}

std::string hzhu_gen_int_to_string(int n)
{
	return std::to_string(n);
}

std::string hzhu_gen_int_to_string(int n, int pad)
{
	std::string re = std::to_string(n);
	int l = re.length();
	std::string tmp;
	for (int i = 0; i < pad - l; i++)
	{
		tmp.append("0");
	}

	return tmp.append(re);
}

#endif // !hzhu_gen_h
