# NPD-micro-saccade-detection

What is contained in this repository is the algorithm (C++) that detects (micro)saccades in eye-gaze signal.

  The algorithm is built and tested with visual studio 2017 with intel compiler 2019. The GNU Scientific Library (GSL) version 2.5 and OpenMP are used in our project.

  - For the installation of the GSL, please visit https://www.gnu.org/software/gsl/. (version 2.5 or higher required)
  - For the installation of the intel compiler, please visit https://software.intel.com/en-us/parallel-studio-xe. (If not using the intel compiler, proper configuration would be required before compiling the project due to the use of the header: omp.h)

  In the directory NPD-micro-saccade-detection/NPD_v5, the following files are contained:

  - main.cpp
  - hzhu_npd.h (The NPD for a section of the gaze signal)
  - hzhu_npd_trial.h (The NPD for a trial of gaze signal)
  - hzhu_mat.h (Handling GSL matrix)
  - hzhu_gen.h (Handling OpenMP and auxiliary tasks)

  By default, our algorithm reads in files X.csv and Y.csv, where the X and Y coordinates of the tracked eye-gaze positions are stored in each of the csv files. Data in the X.csv and Y.csv must be sqaure matrices of same dimensions, and each trail of the tracked eye-gaze positions forms the rows of the input file. Our algorithm can perform the NPD on each row of the dataset independently.
  
  The output of the NPD regarding the [i] trial (ith row of the data files) are the following:
  
  - Result_[i]_detail.csv
  
    The detailed detection result, a N by 10 matrix. N is the total number of detected events of the [i] trial. For detailed explanation of the meaning of each element in a row, please check file _hzhu_npd.h_ for function _hzhu_mat hzhu_npd_results(hzhu_mat &detect, hzhu_mat &result_X, hzhu_mat &result_Y)_.
    
  - Result_[i]_detect.csv
  
    The detection result in the form of 0s and 1s, where 0 indicates there is no saccadic event and 1 indicates there is a event
  
  - Result_[i]data_x.csv
  
    Same as the ith row in X.csv
  
  - Result_[i]data_y.csv
    
    Same as the ith row in Y.csv
  
  - Result_[i]noise_x_raw.csv
  
    Estimated variance as a function of n for horizontal movement of the eye
  
  - Result_[i]noise_y_raw.csv
  
    Estimated variance as a function of n for vertical movement of the eye
  
  - Result_[i]result_x.csv
  
    Unprocessed result regarding the horizontal movement of the eye, for more information please check 
  
  - Result_[i]result_y.csv
  - Result-[i]SIGMA_inv_x.csv
  - Result-[i]SIGMA_inv_y.csv
