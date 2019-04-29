
#include "hzhu_npd_trial.h"
#include <string>
#include <omp.h>

int main()
{
	hzhu_gen orz(-1);
	orz.out("-- Program strated\r\n");

	// Loading data ==================================================================================

	std::string filenameX = "X.csv";
	std::string filenameY = "Y.csv";
	std::string filenameL = "L.csv";

	orz.out(filenameX);
	orz.out(filenameY);
	orz.out(filenameL);
	orz.out();

	hzhu_mat source_data_X(filenameX);
	hzhu_mat source_data_Y(filenameY);
	hzhu_mat data_use_X = source_data_X.transpose();
	hzhu_mat data_use_Y = source_data_Y.transpose();
	hzhu_mat d_source_X = data_use_X.diff();
	hzhu_mat d_source_Y = data_use_Y.diff();

	hzhu_mat data_label = hzhu_mat(filenameL).transpose();

	orz.out("-- Data loaded: ");
	std::cout << "\t\t" << data_use_X[-1] << " x " << data_use_X[-2] << " matrix\r\n" << std::endl;

	int data_N = data_use_X[-1];
	int trail_N = data_use_X[-2];

	// Pre-processing ==================================================================================

	hzhu_mat x_abnormal = hzhu_npd_abnormal(data_use_X);
	hzhu_mat y_abnormal = hzhu_npd_abnormal(data_use_Y);
	hzhu_mat abnormal = x_abnormal + y_abnormal + hzhu_npd_incomplete(data_label);

	// NPD =============================================================================================
	orz.start_counter();

	int noise_n1 = 1;
	int noise_n2 = 6;
	// int noise_chunk = 100;
	int noise_chunk = 30;
	// int chunk_size = 35;
	int chunk_size = 15;

	double step_size = 10.0;
	double error_tol = 1.0e-3;
	int max_iter = 2000;
	double dm = 0.95;
	//double start_x1 = 5000;
	//double start_x2 = 900;
	double start_x1 = 10.0;
	double start_x2 = 2.0;
	int result_n = 10;

	hzhu_npd_trial **NPD = new hzhu_npd_trial*[trail_N];

#pragma omp parallel for
	for (int i = 0; i < trail_N; i++)
	{
		NPD[i] = new hzhu_npd_trial(data_use_X.get_col(i), data_use_Y.get_col(i));

		if (abnormal.data->data[i] == 0.0)
		{
			NPD[i]->init(step_size, error_tol, max_iter, dm, chunk_size, result_n);
			NPD[i]->noise_est(noise_chunk, noise_n1, noise_n2);
			NPD[i]->process(start_x1, start_x2, start_x1, start_x2);
		}
	}

	for (int i = 0; i < trail_N; i++)
	{
		if (abnormal.data->data[i] == 0.0)
		{
			std::string local_name("Result_" + hzhu_gen_int_to_string(i + 1));
			//NPD[i]->save_all(local_name.c_str());
			hzhu_mat D = NPD[i]->detect(0.0001, NPD[0]->median_abs_theta1(), 4.0, chunk_size);
			//D.save_to_file(local_name + "_detect.csv");
			hzhu_mat R = hzhu_npd_results(D, *NPD[i]->result_x, *NPD[i]->result_y);
			R.save_to_file(local_name + "_detail.csv");
		}
	}

	orz.disp_counter_time();
	orz.out("-- NDP completed\r\n");

	// End ===============================================================================================================

	for (int i = 0; i < trail_N; i++)
	{
		delete NPD[i];
	}
	delete[] NPD;

	orz.disp_runtime();
	system("PAUSE");
	return 1;
}
