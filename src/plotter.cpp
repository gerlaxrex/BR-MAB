#include "plotter.h"
#include <string>
#include <stdlib.h>

void Plotter::plot_experiment(int num_exp, string name_exp) {
  //string command = "python scripts/plot_results.py " + to_string(num_exp);
  string command = "python scripts/plot_results.py " + to_string(num_exp) + " " + name_exp;
	const char* c = command.c_str();
	system(c);
}
