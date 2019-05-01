# NPD-micro-saccade-detection

What is contained in this repository is the algorithm (C++) that detects (micro)saccades in eye-gaze signal.

  The algorithm is built and tested with visual studio 2017 with intel compiler 2019. The GNU Scientific Library (GSL) version 2.5 and OpenMP are used in our project.

  - For the installation of the GSL, please visit https://www.gnu.org/software/gsl/. (version 2.5 or higher required)
  - For the installation of the intel compiler, please visit https://software.intel.com/en-us/parallel-studio-xe. (If not using the intel compiler, proper configuration would be required before compiling the project due to the use of the header: omp.h)

