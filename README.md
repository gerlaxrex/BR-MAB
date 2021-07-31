[![made-with-c++](https://img.shields.io/badge/Made%20with-C%2B%2B-red?logo=c%2B%2B)](https://www.cplusplus.com//)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-blue?logo=python)](https://www.python.org/)
# BR-MAB
### Code for the BR-MAB framework applicable to non-stationary multi-armed bandits.

📜 **Paper**: <a href="https://2021.ecmlpkdd.org/wp-content/uploads/2021/07/sub_1017.pdf">Exploiting History Data for Nonstationary MAB</a>

#### Purpose

This repository contains the code used for testing and performing the experiments described in the paper _**"Exploiting History Data for Nonstationary Multi-armed Bandits"**_ by _G. Re, F. Chiusano, F. Trovò, D. Carrera, G. Boracchi and M. Restelli_.
There are some "ready-to-use" experiments that can be ran simply by following the instructions written in the last section of this description, but the user could define its own settings (only a basic knowledge of Python is required in order to do it).

#### Code

The core of the code has been taken from <a href="https://github.com/fabiochiusano/SwitchingBandit">here</a> (repository of Fabio Chiusano).
This contains various classes in order to define MABs and some algorithms, both stationary and nonstationary, already used in literature (_UCB1[1], Thompson Sampling[2], CUSUM-UCB[3] ..._).
In addition, here we added:
  - Some missing algorithms, like Bayes-UCB[4] and M-UCB[5].
  - New Change Detection Algorithms (GRL[6], Monitor[5]).
  - Statistical and similarity tests.
  - the framework proposed in our paper, the **BR-MAB**.

#### Instructions (RR)

Reproducible Research (RR) instructions:

  1. Run the script ```'scripts/exp_demo_generator.py'``` in order to generate the experimental file (default is the synthetic setting in the paper).
  2. ```make```
  2. Run the command ```'./bin/main'``` and begin the experiment.
  3. Images of the experiments will be saved in the ```'results'``` folder with a text file of the final regret results.

#### Bibliography


[1] Auer, P., Cesa-Bianchi, N., and Fischer, P. (2002). Finite-time analysis of the multiarmed bandit problem. Machine learning, 47(2-3):235–256. <br>
[2] Thompson, W. R. (1933). On the likelihood that one unknown probability exceeds another in view of the evidence of two samples. <br>
[3] Liu, F., Lee, J., and Shroff, N. B. (2018). A change-detection based framework for piecewise-stationary multi-armed bandit problem. In AAAI.<br>
[4] Kaufmann, Emilie, Olivier Cappé, and Aurélien Garivier. "On Bayesian upper confidence bounds for bandit problems." Artificial intelligence and statistics. PMLR, 2012.<br>
[5] Cao, Y., Wen, Z., Kveton, B., and Xie, Y. (2019). Nearly optimal adaptive procedure with change detection for piecewise-stationary bandit.<br>
[6]Neyman, J. and Pearson, E. S. (1933). On the problem of the most efficient tests of statistical hypotheses. Philosophical Transactions of the Royal Society of London. Series A, Containing Papers of a Mathematical or Physical Character, 231:289–337.<br>
[7] Auer, P., Cesa-Bianchi, N., and Fischer, P. (2002). Finite-time analysis of the multiarmed bandit problem. Machine learning, 47(2-3):235–256.<br>


