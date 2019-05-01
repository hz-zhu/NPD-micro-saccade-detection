#pragma once

#include "hzhu_mat.h"
#include "hzhu_gen.h"
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_filter.h>
#include <gsl/gsl_fit.h>

// Functions =================================================================

hzhu_mat hzhu_npd_func_h_s(hzhu_mat &t, double theta1, double theta2, double theta3, double dm);
double hzhu_npd_ML_func(const gsl_vector * x, void * params);

hzhu_mat hzhu_npd_detect(hzhu_mat &result, double alpha, double sd, double amp, int strengh);

hzhu_mat hzhu_npd_event(hzhu_mat &detect);
hzhu_mat hzhu_npd_events_compare(hzhu_mat &est, hzhu_mat &truth);

hzhu_mat hzhu_npd_results(hzhu_mat &detect, hzhu_mat &result_X, hzhu_mat &result_Y);

hzhu_mat hzhu_npd_noise_est(hzhu_mat &data, int n1, int n2, int chunk_size);
hzhu_mat hzhu_npd_abnormal(hzhu_mat &data);
hzhu_mat hzhu_npd_incomplete(hzhu_mat &data);

hzhu_mat hzhu_npd_SIGMA_inv(int n, double var_c, double var_v);

// Class =====================================================================

class hzhu_npd
{
public:
	hzhu_npd(hzhu_mat &source, hzhu_mat &SIGMA_inverse);
	~hzhu_npd();

	void init(double DM, int max_interations, double error_tolerance, double step);
	int solve(gsl_vector *x);
	hzhu_mat test(double alpha);
	hzhu_mat test();
	hzhu_mat test(hzhu_mat alphas);
	hzhu_mat get_theta();

public:
	hzhu_mat *data;
	hzhu_mat *SIGMA_inv;

	double *theta_est = NULL;
	double *NPD_T = NULL;
	double *NPD_T_var = NULL;
	double *min_value = NULL;

	int flag;// -1: not init; 0: ready for minimization; 1: minimized but not optimal; 2: mimimized and optimal
	double dm;
	int max_inter;
	double error_tol;
	double step_size;

	gsl_multimin_fminimizer *minimizer = NULL;
	gsl_multimin_function *min_func = NULL;
};

// Class function realization ====================================================
hzhu_npd::hzhu_npd(hzhu_mat &source, hzhu_mat &SIGMA_inverse)
{
	data = new hzhu_mat(source);
	SIGMA_inv = new hzhu_mat(SIGMA_inverse);
	flag = -1;
}

hzhu_npd::~hzhu_npd()
{
	hzhu_gen_free(data);
	hzhu_gen_free(SIGMA_inv);
	hzhu_gen_free(theta_est);
	hzhu_gen_free(NPD_T);
	hzhu_gen_free(NPD_T_var);
	hzhu_gen_free(min_value);

	if (min_func != NULL) delete min_func;
	if (minimizer != NULL) gsl_multimin_fminimizer_free(minimizer);
}

void hzhu_npd::init(double DM, int max_interations, double error_tolerance, double step)
{
	if (min_func == NULL) min_func = new gsl_multimin_function;
	if (minimizer == NULL) minimizer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, 2);
	min_func->n = 2;
	min_func->params = this;
	min_func->f = &hzhu_npd_ML_func;

	step_size = step;
	dm = DM;
	max_inter = max_interations;
	error_tol = error_tolerance;

	flag = 0;
}

int hzhu_npd::solve(gsl_vector *x)
{
	if (flag == -1)
	{
		if (HZHU_WARNING_DISP) std::cout << "Need to call init() first @ void hzhu_npd::solve()" << std::endl;
		return 0;
	}

	gsl_vector *step_size_vec = gsl_vector_alloc(2);
	gsl_vector_set_all(step_size_vec, step_size);
	gsl_multimin_fminimizer_set(minimizer, min_func, x, step_size_vec);

	gsl_vector_free(step_size_vec);

	int iter = 0;
	int status;
	double size;

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(minimizer);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size(minimizer);
		status = gsl_multimin_test_size(size, error_tol);

		if (status == GSL_SUCCESS)
		{
			//printf("converged to minimum at\n");
			theta_est = new double[3];
			theta_est[1] = gsl_vector_get(minimizer->x, 0);
			theta_est[2] = gsl_vector_get(minimizer->x, 1);

			int n_data = this->data->data->size1;
			hzhu_mat t_tmp = hzhu_equal_space(0.0, 1.0, n_data);
			hzhu_mat h_tmp = hzhu_npd_func_h_s(t_tmp, 1.0, theta_est[1], theta_est[2], dm);

			hzhu_mat mu_bar = h_tmp.diff();

			hzhu_mat tmp = mu_bar.transpose().product(*this->SIGMA_inv);
			double numerator = tmp.product(*this->data).data->data[0];
			double denominator = tmp.product(mu_bar).data->data[0];

			theta_est[0] = numerator / denominator;

			NPD_T = new double;
			NPD_T[0] = -minimizer->fval;
			NPD_T_var = new double;
			NPD_T_var[0] = NPD_T[0];
			min_value = new double;
			min_value[0] = minimizer->fval;

			flag = 2;
			return 2;
		}

		/*
		printf("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
			iter,
			gsl_vector_get(minimizer->x, 0),
			gsl_vector_get(minimizer->x, 1),
			minimizer->fval, size);
		*/

	} while (status == GSL_CONTINUE && iter < max_inter);

	theta_est = new double[3];
	theta_est[1] = gsl_vector_get(minimizer->x, 0);
	theta_est[2] = gsl_vector_get(minimizer->x, 1);

	int n_data = this->data->data->size1;
	hzhu_mat t_tmp = hzhu_equal_space(0.0, 1.0, n_data);
	hzhu_mat h_tmp = hzhu_npd_func_h_s(t_tmp, 1.0, theta_est[1], theta_est[2], dm);

	hzhu_mat mu_bar = h_tmp.diff();

	hzhu_mat tmp = mu_bar.transpose().product(*this->SIGMA_inv);
	double numerator = (tmp.product(*this->data)).data->data[0];
	double denominator = (tmp.product(mu_bar)).data->data[0];

	theta_est[0] = numerator / denominator;

	NPD_T = new double;
	NPD_T[0] = -minimizer->fval;
	NPD_T_var = new double;
	NPD_T_var[0] = NPD_T[0];
	min_value = new double;
	min_value[0] = minimizer->fval;

	flag = 1;
	return 1;
}

hzhu_mat hzhu_npd::test(double alpha)
{
	if (flag == -1 || flag == 0)
	{
		if (HZHU_WARNING_DISP) std::cout << "NULL matrix returned @ hzhu_mat hzhu_npd::test(double alpha)" << std::endl;
		return hzhu_mat();
	}

	hzhu_mat re(2, 1);
	re.data->data[1] = gsl_cdf_ugaussian_Pinv(1.0 - alpha)*pow(NPD_T_var[0], 0.5);
	if (NPD_T[0] > re.data->data[1]) re.data->data[0] = 1.0;// It is a saccadic event
	else re.data->data[0] = 0.0;// It is not a saccadic event

	return re;
}

hzhu_mat hzhu_npd::test()
{
	return test(0.05);
}

hzhu_mat hzhu_npd::test(hzhu_mat alphas)
{
	if (alphas.data == NULL)
	{
		if (HZHU_WARNING_DISP) std::cout << "NULL matrix returned @ hzhu_mat hzhu_npd::test(hzhu_mat alphas)" << std::endl;
		return hzhu_mat();
	}

	if (flag == -1 || flag == 0)
	{
		if (HZHU_WARNING_DISP) std::cout << "NULL matrix returned @ hzhu_mat hzhu_npd::test(hzhu_mat alphas)" << std::endl;
		return hzhu_mat();
	}

	int N = alphas.data->size1*alphas.data->size2;
	hzhu_mat re((int)2, N);
	double tmp = pow(NPD_T_var[0], 0.5);

	for (int i = 0; i < N; i++)
	{
		re.data->data[i + N] = gsl_cdf_ugaussian_Pinv(1.0 - alphas.data->data[i])*tmp;
		if (NPD_T[0] > re.data->data[i + N]) re.data->data[i] = 1.0;// It is a saccadic event
		else re.data->data[i] = 0.0;// It is not a saccadic event
	}

	return re;
}

hzhu_mat hzhu_npd::get_theta()
{
	if (theta_est == NULL)
	{
		if (HZHU_WARNING_DISP) std::cout << "Theta not calculated yet @ hzhu_mat hzhu_npd::get_theta()" << std::endl;
		return hzhu_mat();
	}

	hzhu_mat re(1, 3);
	re.data->data[0] = theta_est[0];
	re.data->data[1] = theta_est[1];
	re.data->data[2] = theta_est[2];

	return re;
}

// Function realization ===========================================================

hzhu_mat hzhu_npd_func_h_s(hzhu_mat &t, double theta1, double theta2, double theta3, double dm)
{
	if (t.data == NULL)
	{
		if (HZHU_WARNING_DISP) std::cout << "Null matrix t @ hzhu_mat hzhu_npd_func_h(hzhu_mat &t, double theta1, double theta2, double theta3, double dm)" << std::endl;
		return hzhu_mat();
	}

	int n1 = t.data->size1;
	int n2 = t.data->size2;
	int N = n1 * n2;
	hzhu_mat re(n1, n2);

	double theta_3_inv = 1.0 / theta3;
	double t1 = theta2 * pow(-log(0.5 + 0.5*dm), theta_3_inv);

	for (int i = 0; i < N; i++)
	{
		double x = t.data->data[i] + t1;
		if (x <= 0.0) re.data->data[i] = 0.0;
		else re.data->data[i] = theta1 * (1.0 - exp(-pow(x / theta2, theta3)));
	}

	return re;
}

double hzhu_npd_ML_func(const gsl_vector * x, void * params)
{
	hzhu_npd *p = (hzhu_npd *)params;
	double theta2 = x->data[0];
	double theta3 = x->data[1];
	double dm = p->dm;

	double theta_3_inv = 1.0 / theta3;
	double theta_3_f = 1.0 - theta_3_inv;

	double t1 = theta2 * pow(-log(0.5 + 0.5*dm), theta_3_inv);
	double t2 = theta2 * pow(-log(0.5 - 0.5*dm), theta_3_inv);
	double d = t2 - t1;
	double t_vmax = theta2 * pow(theta_3_f, theta_3_inv) - t1;

	int n_data = p->data->data->size1;

	if (theta2 < 1.0 || theta3 <= 0.0) return abs(1.0 - theta2) + abs(theta3) + 1.0e99;
	if (d <= 0.0) return abs(1.0 - theta2) + abs(theta3) + 3.0 - d + 1.0e99;
	if (d >= n_data) return abs(1.0 - theta2) + abs(theta3) + 3.0 + d + 1.0e99;
	if (t_vmax > d) return abs(1.0 - theta2) + abs(theta3) + 3.0 + d + t_vmax + 1.0e99;

	hzhu_mat t_tmp = hzhu_equal_space(0.0, 1.0, n_data);
	hzhu_mat h_tmp = hzhu_npd_func_h_s(t_tmp, 1.0, theta2, theta3, dm);

	hzhu_mat mu_bar = h_tmp.diff();

	hzhu_mat tmp = mu_bar.transpose().product(*p->SIGMA_inv);
	double numerator = tmp.product(*p->data).data->data[0];
	double denominator = tmp.product(mu_bar).data->data[0];

	double re = -numerator * numerator / denominator;
	if (isfinite(re)) return re;
	else return 2.0e99;
}

hzhu_mat hzhu_npd_detect(hzhu_mat &result, double alpha, double sd, double amp, int strengh)
{
	if (result.data == NULL)
	{
		if (HZHU_WARNING_DISP) std::cout << "Null matrix result @ hzhu_mat hzhu_npd_detect(hzhu_mat &result, double alpha, double sd, double amp, int strengh)" << std::endl;
		return hzhu_mat();
	}

	alpha = alpha < 0.0 ? 0.0 : alpha;
	alpha = alpha > 1.0 ? 1.0 : alpha;

	strengh = strengh <= 0 ? 1 : strengh;

	int N1 = result.data->size1;
	int N2 = result.data->size2;

	gsl_vector *tmp = gsl_vector_alloc(N1);
	double alpha_factor = gsl_cdf_ugaussian_Pinv(1.0 - alpha);

	switch (gsl_isinf(alpha_factor))
	{
	case +1: alpha_factor = 1.0e299;
		break;
	case -1: alpha_factor = -1.0e299;
		break;
	}

	double amp_factor = sd * amp;

	hzhu_mat re(N1, 1);

	for (int i = 0; i < N1; i++)
	{
		tmp->data[i] = result.data->data[i*N2 + 3] > alpha_factor * pow(result.data->data[i*N2 + 4], 0.5) ? 1.0 : 0.0;
		tmp->data[i] = abs(result.data->data[i*N2]) > amp_factor ? tmp->data[i] : 0.0;
	}

	gsl_filter_median_workspace *filter = gsl_filter_median_alloc(strengh);
	gsl_filter_median(GSL_FILTER_END_PADZERO, tmp, tmp, filter);

	for (int i = 0; i < N1; i++) re.data->data[i] = tmp->data[i];

	gsl_filter_median_free(filter);
	gsl_vector_free(tmp);
	return re;
}

hzhu_mat hzhu_npd_event(hzhu_mat &detect)
{
	if (detect.data == NULL)
	{
		if (HZHU_WARNING_DISP) std::cout << "Null matrix detect @ hzhu_mat hzhu_npd_event(hzhu_mat &detect)" << std::endl;
		return hzhu_mat();
	}

	int n1 = detect.data->size1;
	int n2 = detect.data->size2;

	if (n2 != 1)
	{
		if (HZHU_WARNING_DISP) std::cout << "Invalid input @ hzhu_mat hzhu_npd_event(hzhu_mat &detect)" << std::endl;
		return hzhu_mat();
	}

	hzhu_mat local = detect;
	local.append_top(0.0);
	local.append_bottom(0.0);

	hzhu_mat d_local = local.diff();
	hzhu_mat tmp = d_local ^ 2.0;
	int counter = tmp.sum_all().data->data[0] / 2.0;

	if (counter == 0.0) return hzhu_mat(-1.0);

	hzhu_mat re(counter, 2);

	int j = 0;
	for (int i = 0; i < n1 + 1; i++)
	{
		if (d_local.data->data[i] == 1.0)
		{
			re.data->data[j * 2] = i;
			continue;
		}
		if (d_local.data->data[i] == -1.0)
		{
			re.data->data[j * 2 + 1] = i - 1;
			j++;
		}
	}

	return re;
}

hzhu_mat hzhu_npd_noise_est(hzhu_mat &data, int n1, int n2, int chunk_size)
{
	int n = n2 - n1 + 1;
	int N = n >= 6 ? n : 6;
	hzhu_mat re(0.0, 3, N);

	hzhu_mat x = hzhu_equal_space(n1, 1.0, n2);
	hzhu_mat y(n, 1);

	for (int i = 0; i < n; i++)
	{
		hzhu_mat d = data.diff(n1 + i);
		hzhu_mat holder((int)d.data->size1 - chunk_size + 1, 1);
		for (int j = 0; j < d.data->size1 - chunk_size + 1; j++)
		{
			holder.data->data[j] = d.get_sub(j, 0, j + chunk_size - 1, 0).uvar_all().data->data[0];
		}
		y.data->data[i] = holder.min_value_all().data->data[0];
	}

	double var_est_c0, var_est_c1, var_est_cov00, var_est_cov01, var_est_cov11, var_est_chisq;
	gsl_fit_linear(x.data->data, 1, y.data->data, 1, n, &var_est_c0, &var_est_c1, &var_est_cov00, &var_est_cov01, &var_est_cov11, &var_est_chisq);

	double var_v_est = var_est_c0 / 2.0;
	double var_w_est = var_est_c1;
	if (var_v_est < 0.0) var_v_est = 0.0;
	if (var_w_est < 0.0) var_w_est = 0.0;
	double var_c_est = var_v_est * 2.0 + var_w_est;
	if (var_v_est < 0.0 && var_w_est < 0.0) var_c_est = 1.0e-7;

	re.data->data[0] = var_v_est;
	re.data->data[1] = var_w_est;
	re.data->data[2] = var_c_est;

	re.data->data[3] = var_est_c0;
	re.data->data[4] = var_est_c1;
	re.data->data[5] = var_est_chisq;

	for (int i = 0; i < n; i++)
	{
		re.data->data[N + i] = x.data->data[i];
		re.data->data[2 * N + i] = y.data->data[i];
	}

	return re;
}

hzhu_mat hzhu_npd_abnormal(hzhu_mat &data)
{
	int n1 = data.data->size1;
	int n2 = data.data->size2;

	hzhu_mat re(0.0, n2, 1);

	for (int i = 0; i < n2; i++)
	{
		hzhu_mat local = data.get_col(i);
		for (int j = 0; j < n1; j++)
		{
			if (!gsl_finite(local.data->data[j]))
			{
				re.data->data[i] = 1.0;
				break;
			}
		}
		for (int j = 0; j < n1 - 20; j++)
		{
			if (local.get_sub(j, 0, j + 20, 0).var_all().data->data[0] == 0.0)
			{
				re.data->data[i] = 1.0;
				break;
			}
		}
	}
	return re;
}

hzhu_mat hzhu_npd_incomplete(hzhu_mat &data)
{
	int n1 = data.data->size1;
	int n2 = data.data->size2;

	hzhu_mat re(0.0, n2, 1);

	for (int i = 0; i < n2; i++)
	{
		if (data.data->data[i] == 1.0 || data.data->data[((n1 - 1)*n2) + i] == 1.0) re.data->data[i] = 1.0;
	}
	return re;
}

hzhu_mat hzhu_npd_SIGMA_inv(int n, double var_c, double var_v)
{
	hzhu_mat tmp(0.0, n, n);
	for (int i = 0; i < n; i++)
	{
		tmp.data->data[i + i * n] = var_c;
		if (i != n - 1)
		{
			tmp.data->data[i + 1 + i * n] = -var_v;
			tmp.data->data[i + (1 + i)*n] = -var_v;
		}
	}

	hzhu_mat_LU_inv Tmp(tmp);
	Tmp.compute_LU();
	Tmp.compute_inv();

	return Tmp.get_inv();
}

hzhu_mat hzhu_npd_events_compare(hzhu_mat &est, hzhu_mat &truth)
{
	if (est.data == NULL || truth.data == NULL)
	{
		if (HZHU_WARNING_DISP) std::cout << "NULL input data @ hzhu_mat hzhu_npd_event_compare(hzhu_mat &est, hzhu_mat &truth)" << std::endl;
		return hzhu_mat();
	}

	int est_n1 = est.data->size1;
	int est_n2 = est.data->size2;

	int truth_n1 = truth.data->size1;
	int truth_n2 = truth.data->size2;

	hzhu_mat re(4, 1);// #true events, #correctedly detected in all true events, # all detected events, # wrongly detected events in all detected events

	if (est_n2 == 1 && truth_n2 != 1)
	{
		re.data->data[0] = truth_n1;
		re.data->data[1] = 0.0;
		re.data->data[2] = 0.0;
		re.data->data[3] = 0.0;
		return re;
	}

	if (est_n2 == 1 && truth_n2 == 1)
	{
		re.data->data[0] = 0.0;
		re.data->data[1] = 0.0;
		re.data->data[2] = 0.0;
		re.data->data[3] = 0.0;
		return re;
	}

	if (est_n2 != 1 && truth_n2 == 1)
	{
		re.data->data[0] = 0.0;
		re.data->data[1] = 0.0;
		re.data->data[2] = est_n1;
		re.data->data[3] = est_n1;
		return re;
	}

	hzhu_mat flag(1.0, truth_n1, 1);
	int true_events = truth_n1;
	int correctly_detected_events = 0;
	int all_detected_events = est_n1;
	int wrongly_detected_events = 0;

	for (int i = 0; i < est_n1; i++)
	{
		int a = est.data->data[2 * i];
		int b = est.data->data[2 * i + 1];
		for (int j = 0; j < truth_n1; j++)
		{
			if (flag.data->data[j] == 1.0)
			{
				int A = truth.data->data[2 * j];
				int B = truth.data->data[2 * j + 1];

				if (a <= B && a >= A)
				{
					correctly_detected_events++;
					flag.data->data[j] = 0.0;
					break;
				}
				if (b <= B && b >= A)
				{
					correctly_detected_events++;
					flag.data->data[j] = 0.0;
					break;
				}
				if (A <= b && A >= a)
				{
					correctly_detected_events++;
					flag.data->data[j] = 0.0;
					break;
				}
				if (B <= b && B >= a)
				{
					correctly_detected_events++;
					flag.data->data[j] = 0.0;
					break;
				}
			}
		}
	}
	wrongly_detected_events = all_detected_events - correctly_detected_events;

	re.data->data[0] = true_events;
	re.data->data[1] = correctly_detected_events;
	re.data->data[2] = all_detected_events;
	re.data->data[3] = wrongly_detected_events;
	return re;
}

hzhu_mat hzhu_npd_results(hzhu_mat &detect, hzhu_mat &result_X, hzhu_mat &result_Y)
{
	if (detect.data == NULL)
	{
		if (HZHU_WARNING_DISP) std::cout << "Null matrix detect @ hzhu_mat hzhu_npd_event(hzhu_mat &detect)" << std::endl;
		return hzhu_mat();
	}

	int n1 = detect.data->size1;
	int n2 = detect.data->size2;

	if (n2 != 1)
	{
		if (HZHU_WARNING_DISP) std::cout << "Invalid input @ hzhu_mat hzhu_npd_event(hzhu_mat &detect)" << std::endl;
		return hzhu_mat();
	}

	hzhu_mat local = detect;
	local.append_top(0.0);
	local.append_bottom(0.0);

	hzhu_mat d_local = local.diff();
	hzhu_mat tmp = d_local ^ 2.0;
	int counter = tmp.sum_all().data->data[0] / 2.0;

	if (counter == 0.0) return hzhu_mat(-1.0);

	hzhu_mat re(counter, 10);

	int j = 0;
	for (int i = 0; i < n1 + 1; i++)
	{
		if (d_local.data->data[i] == 1.0)
		{
			re.data->data[j * 10] = i;
			continue;
		}
		if (d_local.data->data[i] == -1.0)
		{
			re.data->data[j * 10 + 1] = i - 1;
			j++;
		}
	}

	int N = result_X.data->size2;
	for (int j = 0; j < counter; j++)
	{
		int a = re.data->data[j * 10];// starting point of the detection chunk
		int b = re.data->data[j * 10 + 1]; // Ending point of the detection chunk
		
		hzhu_mat tmp = result_X.get_sub(a, 7, b, 7);
		int max_index = tmp.min_index_all().data->data[0];
		re.data->data[j * 10 + 2] = max_index + a;//index of start point t1 for X (determined by the best fitted data section)
		re.data->data[j * 10 + 3] = result_X.data->data[(max_index + a) * N];//theta1 for X
		re.data->data[j * 10 + 4] = result_X.data->data[(max_index + a) * N + 1];//theta2 for X
		re.data->data[j * 10 + 5] = result_X.data->data[(max_index + a) * N + 2];//theta3 for X

		hzhu_mat Tmp = result_Y.get_sub(a, 7, b, 7);
		max_index = Tmp.min_index_all().data->data[0];
		re.data->data[j * 10 + 6] = max_index + a;//index of start point t1 for Y (determined by the best fitted data section)
		re.data->data[j * 10 + 7] = result_Y.data->data[(max_index + a) * N];//theta1 for Y
		re.data->data[j * 10 + 8] = result_Y.data->data[(max_index + a) * N + 1];//theta2 for Y
		re.data->data[j * 10 + 9] = result_Y.data->data[(max_index + a) * N + 2];//theta3 for Y
	}

	return re;
}