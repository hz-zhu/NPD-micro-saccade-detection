#pragma once

#include "hzhu_npd.h"

class hzhu_npd_trial
{
public:
	hzhu_npd_trial(hzhu_mat x, hzhu_mat y);
	hzhu_npd_trial(hzhu_npd_trial& source);
	~hzhu_npd_trial();

	void init(double stepSize, double errorTol, int maxIter, double dM, int chunkSize, int resultN);
	void noise_est(int noise_chunk_Size, int a, int b);
	void process(double x1, double x2, double y1, double y2);
	void save_all(const char *name, const char *path);
	void save_all(const char *name);
	hzhu_mat detect(double alpha, double magnitude, int strength);
	hzhu_mat median_abs_theta1();
	hzhu_mat detect(double alpha, hzhu_mat median_value, double factor, int strength);

	double step_size = 10;
	double error_tol = 1e-3;
	int max_iter = 2000;
	double dm = 0.95;

	int chunk_size = 35;
	int result_n = 10;
	int noise_chunk_size = 120;

	hzhu_mat *noise_xy = NULL;
	hzhu_mat *noise_x_raw = NULL;
	hzhu_mat *noise_y_raw = NULL;
	hzhu_mat *result_x = NULL;
	hzhu_mat *result_y = NULL;
	hzhu_mat *data_x = NULL;
	hzhu_mat *data_y = NULL;
	hzhu_mat *SIGMA_inv_x = NULL;
	hzhu_mat *SIGMA_inv_y = NULL;
};

// Realization
hzhu_npd_trial::hzhu_npd_trial(hzhu_mat x, hzhu_mat y)
{
	if (x.data == NULL || y.data == NULL)
	{
		if (HZHU_WARNING_DISP) std::cout << "Null input @ hzhu_npd_trial::hzhu_npd_trial(hzhu_mat x, hzhu_mat y)" << std::endl;
		return;
	}

	if (x.data->size1 != y.data->size1 || x.data->size2 != 1 || y.data->size2 != 1)
	{
		if (HZHU_WARNING_DISP) std::cout << "Input invalid (dimension) @ hzhu_npd_trial::hzhu_npd_trial(hzhu_mat x, hzhu_mat y)" << std::endl;
		return;
	}

	data_x = new hzhu_mat(x);
	data_y = new hzhu_mat(y);
}

hzhu_npd_trial::hzhu_npd_trial(hzhu_npd_trial& source)
{
	step_size = source.step_size;
	error_tol = source.error_tol;
	max_iter = source.max_iter;
	dm = source.dm;

	chunk_size = source.chunk_size;
	result_n = source.result_n;
	noise_chunk_size = source.noise_chunk_size;

	if (source.noise_xy != NULL) noise_xy = new hzhu_mat(*source.noise_xy);
	else noise_xy = NULL;

	if (source.noise_x_raw != NULL) noise_x_raw = new hzhu_mat(*source.noise_x_raw);
	else noise_x_raw = NULL;

	if (source.noise_y_raw != NULL) noise_y_raw = new hzhu_mat(*source.noise_y_raw);
	else noise_y_raw = NULL;

	if (source.result_x != NULL) result_x = new hzhu_mat(*source.result_x);
	else result_x = NULL;

	if (source.result_y != NULL) result_y = new hzhu_mat(*source.result_y);
	else result_y = NULL;

	if (source.data_x != NULL) data_x = new hzhu_mat(*source.data_x);
	else data_x = NULL;

	if (source.data_y != NULL) data_y = new hzhu_mat(*source.data_y);
	else data_y = NULL;

	if (source.SIGMA_inv_x != NULL) SIGMA_inv_x = new hzhu_mat(*source.SIGMA_inv_x);
	else SIGMA_inv_x = NULL;

	if (source.SIGMA_inv_y != NULL) SIGMA_inv_y = new hzhu_mat(*source.SIGMA_inv_y);
	else SIGMA_inv_y = NULL;
}

hzhu_npd_trial::~hzhu_npd_trial()
{
	hzhu_gen_free(noise_xy);
	hzhu_gen_free(noise_y_raw);
	hzhu_gen_free(noise_x_raw);
	hzhu_gen_free(result_x);
	hzhu_gen_free(result_y);
	hzhu_gen_free(data_x);
	hzhu_gen_free(data_y);
	hzhu_gen_free(SIGMA_inv_x);
	hzhu_gen_free(SIGMA_inv_y);
}

void hzhu_npd_trial::init(double stepSize, double errorTol, int maxIter, double dM, int chunkSize, int resultN)
{
	if (stepSize <= 0.0 || errorTol <= 0.0 || maxIter <= 0 || dM <= 0.0 || chunkSize <= 0 || resultN <= 0)
	{
		if (HZHU_WARNING_DISP) std::cout << "Invalid input @ void hzhu_npd_trial::init(double stepSize, double errorTol, int maxIter, double dM, int chunkSize, int resultN)" << std::endl;
		return;
	}

	if (stepSize <= 0.0 || errorTol <= 0.0 || maxIter <= 0 || dM >= 1.0 || chunkSize <= 0 || resultN <= 0)
	{
		if (HZHU_WARNING_DISP) std::cout << "Invalid input @ void hzhu_npd_trial::init(double stepSize, double errorTol, int maxIter, double dM, int chunkSize, int resultN)" << std::endl;
		return;
	}

	step_size = stepSize;
	error_tol = errorTol;
	max_iter = maxIter;
	dm = dM;

	chunk_size = chunkSize;
	result_n = resultN;
}

void hzhu_npd_trial::noise_est(int noise_chunk_Size, int a, int b)
{
	if (noise_chunk_Size <= 2 || a <= 0 || b <= 0 || b - a < 2)
	{
		if (HZHU_WARNING_DISP) std::cout << "Invalid input @ void noise_est(int noise_chunk_size, int a, int b)" << std::endl;
		return;
	}
	noise_chunk_size = noise_chunk_Size;

	int n = b - a + 1;
	int N = n >= 6 ? n : 6;
	noise_x_raw = new hzhu_mat(2, n);
	noise_y_raw = new hzhu_mat(2, n);

	hzhu_mat tmpX = hzhu_npd_noise_est(*data_x, a, b, noise_chunk_size);
	hzhu_mat tmpY = hzhu_npd_noise_est(*data_y, a, b, noise_chunk_size);

	for (int i = 0; i < n; i++)
	{
		noise_x_raw->data->data[i] = a + i;
		noise_y_raw->data->data[i] = a + i;
		noise_x_raw->data->data[i + N] = tmpX.data->data[2 * N + i];
		noise_y_raw->data->data[i + N] = tmpY.data->data[2 * N + i];
	}

	noise_xy = new hzhu_mat(2, 6);
	for (int i = 0; i < 6; i++)
	{
		noise_xy->data->data[i] = tmpX.data->data[i];
		noise_xy->data->data[i + 6] = tmpY.data->data[i];
	}

	SIGMA_inv_x = new hzhu_mat(hzhu_npd_SIGMA_inv(chunk_size, tmpX.data->data[2], tmpX.data->data[0]));
	SIGMA_inv_y = new hzhu_mat(hzhu_npd_SIGMA_inv(chunk_size, tmpY.data->data[2], tmpY.data->data[0]));
}

void hzhu_npd_trial::process(double x1, double x2, double y1, double y2)
{
	if (noise_xy == NULL)
	{
		if (HZHU_WARNING_DISP) std::cout << "Call function noise_est() first @ void hzhu_npd_trial::process(double x1, double x2, double y1, double y2)" << std::endl;
		return;
	}

	int data_N = data_x->data->size1;
	int chunk_num = data_N - chunk_size;
	result_x = new hzhu_mat(0.0, chunk_num, result_n);
	result_y = new hzhu_mat(0.0, chunk_num, result_n);

	hzhu_mat dX = data_x->diff();
	hzhu_mat dY = data_y->diff();

	gsl_vector *x_X = gsl_vector_alloc(2);
	x_X->data[0] = x1;
	x_X->data[1] = x2;

	gsl_vector *x_Y = gsl_vector_alloc(2);
	x_Y->data[0] = y1;
	x_Y->data[1] = y2;

	for (int j = 0; j < chunk_num; j++)
	{
		hzhu_mat chunk_X = dX.get_sub(j, 0, j + chunk_size - 1, 0);
		hzhu_mat chunk_Y = dY.get_sub(j, 0, j + chunk_size - 1, 0);

		hzhu_npd *NPD_X = new hzhu_npd(chunk_X, *SIGMA_inv_x);
		hzhu_npd *NPD_Y = new hzhu_npd(chunk_Y, *SIGMA_inv_y);

		NPD_X->init(dm, max_iter, error_tol, step_size);
		NPD_Y->init(dm, max_iter, error_tol, step_size);

		NPD_X->solve(x_X);
		NPD_Y->solve(x_Y);

		result_x->data->data[j * result_n] = NPD_X->theta_est[0];//theta1 1
		result_x->data->data[j * result_n + 1] = NPD_X->theta_est[1];//theta1 2
		result_x->data->data[j * result_n + 2] = NPD_X->theta_est[2];//theta1 3
		result_x->data->data[j * result_n + 3] = NPD_X->NPD_T[0];//test statistics 4
		result_x->data->data[j * result_n + 4] = NPD_X->NPD_T_var[0];//test parameters 5
		result_x->data->data[j * result_n + 5] = NPD_X->flag;//optimal? 6
		result_x->data->data[j * result_n + 6] = NPD_X->min_value[0];//final result value 7

		result_x->data->data[j * result_n + 7] = gsl_cdf_gaussian_Q(pow(NPD_X->NPD_T[0], 0.5), 100.0);//False alarm probability value 8


		result_y->data->data[j * result_n] = NPD_Y->theta_est[0];//theta1 1
		result_y->data->data[j * result_n + 1] = NPD_Y->theta_est[1];//theta1 2
		result_y->data->data[j * result_n + 2] = NPD_Y->theta_est[2];//theta1 3
		result_y->data->data[j * result_n + 3] = NPD_Y->NPD_T[0];//test statistics 4
		result_y->data->data[j * result_n + 4] = NPD_Y->NPD_T_var[0];//test parameters 5
		result_y->data->data[j * result_n + 5] = NPD_Y->flag;//optimal? 6
		result_y->data->data[j * result_n + 6] = NPD_Y->min_value[0];//final result value 7

		result_y->data->data[j * result_n + 7] = gsl_cdf_gaussian_Q(pow(NPD_Y->NPD_T[0], 0.5), 100.0);//False alarm probability value 8

		delete NPD_X;
		delete NPD_Y;
	}

	gsl_vector_free(x_Y);
	gsl_vector_free(x_X);
}

hzhu_mat hzhu_npd_trial::detect(double alpha, double magnitude, int strength)
{
	if (result_y == NULL || result_x == NULL)
	{
		if (HZHU_WARNING_DISP) std::cout << "Call function process() first @ hzhu_mat hzhu_npd_trial::detect(double alpha, double magnitude, int strength)" << std::endl;
		return hzhu_mat();
	}

	int N = result_x->data->size1;

	double chunk_bias = dm * dm * chunk_size;
	hzhu_mat theta1X = result_x->get_col(0);
	hzhu_mat theta1Y = result_y->get_col(0);

	for (int i = 0; i < N; i++)
	{
		theta1X.data->data[i] = theta1X.data->data[i] < 0.0 ? -theta1X.data->data[i] : theta1X.data->data[i];
		theta1Y.data->data[i] = theta1Y.data->data[i] < 0.0 ? -theta1Y.data->data[i] : theta1Y.data->data[i];
	}

	hzhu_mat x_event = hzhu_npd_detect(*result_x, alpha, theta1X.median_all().data->data[0] * 2.0, magnitude, strength);
	hzhu_mat y_event = hzhu_npd_detect(*result_y, alpha, theta1Y.median_all().data->data[0] * 2.0, magnitude, strength);

	hzhu_mat joint_event = x_event + y_event;
	gsl_vector *v = gsl_vector_alloc(joint_event.data->size1);
	for (int j = 0; j < v->size; j++)
	{
		joint_event.data->data[j] = joint_event.data->data[j] == 0.0 ? 0.0 : 1.0;
		v->data[j] = joint_event.data->data[j];
	}
	gsl_filter_median_workspace *filter = gsl_filter_median_alloc(strength);

	for (int j = 0; j < 5; j++) gsl_filter_median(GSL_FILTER_END_PADZERO, v, v, filter);
	for (int j = 0; j < v->size; j++) joint_event.data->data[j] = v->data[j];

	gsl_filter_median_free(filter);
	gsl_vector_free(v);

	return joint_event;
}

hzhu_mat hzhu_npd_trial::median_abs_theta1()
{
	if (result_y == NULL || result_x == NULL)
	{
		if (HZHU_WARNING_DISP) std::cout << "Call function process() first @ hzhu_mat hzhu_npd_trial::median_abs_theta1()" << std::endl;
		return hzhu_mat();
	}
	int N = result_x->data->size1;
	hzhu_mat re(2, 1);

	hzhu_mat theta1X = result_x->get_col(0);
	hzhu_mat theta1Y = result_y->get_col(0);

	for (int i = 0; i < N; i++)
	{
		theta1X.data->data[i] = theta1X.data->data[i] < 0.0 ? -theta1X.data->data[i] : theta1X.data->data[i];
		theta1Y.data->data[i] = theta1Y.data->data[i] < 0.0 ? -theta1Y.data->data[i] : theta1Y.data->data[i];
	}
	re.data->data[0] = theta1X.median_all().data->data[0];
	re.data->data[1] = theta1Y.median_all().data->data[0];
	return re;
}

hzhu_mat hzhu_npd_trial::detect(double alpha, hzhu_mat median_value, double factor, int strength)
{
	if (result_y == NULL || result_x == NULL)
	{
		if (HZHU_WARNING_DISP) std::cout << "Call function process() first @ hzhu_mat hzhu_npd_trial::detect(double alpha, double magnitude, int strength)" << std::endl;
		return hzhu_mat();
	}

	int N = result_x->data->size1;

	hzhu_mat x_event = hzhu_npd_detect(*result_x, alpha, median_value.data->data[0], factor, strength);
	hzhu_mat y_event = hzhu_npd_detect(*result_y, alpha, median_value.data->data[1], factor, strength);

	hzhu_mat joint_event = x_event + y_event;
	gsl_vector *v = gsl_vector_alloc(joint_event.data->size1);
	for (int j = 0; j < v->size; j++)
	{
		joint_event.data->data[j] = joint_event.data->data[j] == 0.0 ? 0.0 : 1.0;
		v->data[j] = joint_event.data->data[j];
	}
	gsl_filter_median_workspace *filter = gsl_filter_median_alloc(strength);

	for (int j = 0; j < 5; j++) gsl_filter_median(GSL_FILTER_END_PADZERO, v, v, filter);
	for (int j = 0; j < v->size; j++) joint_event.data->data[j] = v->data[j];

	gsl_filter_median_free(filter);
	gsl_vector_free(v);

	return joint_event;
}

void hzhu_npd_trial::save_all(const char *name, const char *path)
{
	std::string NAME(path);
	NAME.append("\\");
	NAME.append(name);

	if (SIGMA_inv_y != NULL)
	{
		std::string tmp(NAME);
		tmp.append("SIGMA_inv_y.csv");
		SIGMA_inv_y->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("SIGMA_inv_y.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (SIGMA_inv_x != NULL)
	{
		std::string tmp(NAME);
		tmp.append("SIGMA_inv_x.csv");
		SIGMA_inv_x->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("SIGMA_inv_x.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (data_y != NULL)
	{
		std::string tmp(NAME);
		tmp.append("data_y.csv");
		data_y->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("data_y.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (data_x != NULL)
	{
		std::string tmp(NAME);
		tmp.append("data_x.csv");
		data_x->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("data_x.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (result_y != NULL)
	{
		std::string tmp(NAME);
		tmp.append("result_y.csv");
		result_y->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("result_y.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (result_x != NULL)
	{
		std::string tmp(NAME);
		tmp.append("result_x.csv");
		result_x->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("result_x.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (noise_y_raw != NULL)
	{
		std::string tmp(NAME);
		tmp.append("noise_y_raw.csv");
		noise_y_raw->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("noise_y_raw.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (noise_xy != NULL)
	{
		std::string tmp(NAME);
		tmp.append("noise_xy.csv");
		noise_xy->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("noise_xy.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (noise_x_raw != NULL)
	{
		std::string tmp(NAME);
		tmp.append("noise_x_raw.csv");
		noise_x_raw->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("noise_x_raw.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}
}

void hzhu_npd_trial::save_all(const char *name)
{
	std::string NAME(name);

	if (SIGMA_inv_y != NULL)
	{
		std::string tmp(NAME);
		tmp.append("SIGMA_inv_y.csv");
		SIGMA_inv_y->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("SIGMA_inv_y.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (SIGMA_inv_x != NULL)
	{
		std::string tmp(NAME);
		tmp.append("SIGMA_inv_x.csv");
		SIGMA_inv_x->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("SIGMA_inv_x.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (data_y != NULL)
	{
		std::string tmp(NAME);
		tmp.append("data_y.csv");
		data_y->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("data_y.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (data_x != NULL)
	{
		std::string tmp(NAME);
		tmp.append("data_x.csv");
		data_x->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("data_x.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (result_y != NULL)
	{
		std::string tmp(NAME);
		tmp.append("result_y.csv");
		result_y->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("result_y.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (result_x != NULL)
	{
		std::string tmp(NAME);
		tmp.append("result_x.csv");
		result_x->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("result_x.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (noise_y_raw != NULL)
	{
		std::string tmp(NAME);
		tmp.append("noise_y_raw.csv");
		noise_y_raw->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("noise_y_raw.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (noise_xy != NULL)
	{
		std::string tmp(NAME);
		tmp.append("noise_xy.csv");
		noise_xy->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("noise_xy.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}

	if (noise_x_raw != NULL)
	{
		std::string tmp(NAME);
		tmp.append("noise_x_raw.csv");
		noise_x_raw->save_to_file(tmp);
	}
	else
	{
		std::string tmp(NAME);
		tmp.append("noise_x_raw.csv");
		hzhu_mat(-1.0).save_to_file(tmp);
	}
}