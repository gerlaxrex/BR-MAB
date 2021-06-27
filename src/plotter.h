#ifndef PLOTTER_H
#define PLOTTER_H
#include<string>

using namespace std;

/**
 * @brief Class adept at plotting the results of the experiments
 */
class Plotter {
public:
  /**
   * @brief plots the results of the experiments inside iìthe images folder
   */
  static void plot_experiment(int num_exp, string name_exp);
};

#endif
