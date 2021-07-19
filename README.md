[![made-with-c++](https://img.shields.io/badge/Made%20with-C%2B%2B-red?logo=c%2B%2B)](https://www.cplusplus.com//)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-blue?logo=python)](https://www.python.org/)
# BR-MAB
### Code for the BR-MAB framework applicable to non-stationary multi-armed bandits.

📜 **Paper**: <a href="https://boracchi.faculty.polimi.it/docs/2021_ECML_Nonstationary_MAB.pdf">Exploiting History Data for Nonstationary MAB</a>

#### Purpose

This repository contains the code used for testing and performing the experiments described in the paper _**"Exploiting History Data for Nonstationary Multi-armed Bandits"**_ by _G.Re, F. Chiusano, F. Trovò, D. Carrera, G. Boracchi and M. Restelli_.
There are some "ready-to-use" experiments that can be ran simply by following the instructions written in the last section of this description, but the user could define its own settings (only a basic knowledge of Python is required in order to do it).

#### Code

The core of the code has been taken from <a href="https://github.com/fabiochiusano/SwitchingBandit">here</a> (repository of Fabio Chiusano).
This contains various classes in order to define MABs and some algorithms, both stationary and nonstationary, already used in literature (_UCB1[1], Thompson Sampling[2], CUSUM-UCB[3] ..._).
In addition, here we added:
  - Some missing algorithms, like Bayes-UCB[4] and M-UCB[5].
  - New Change Detection Algorithms (GRL[6], Monitor[5]).
  - Statistical and similarity tests.
  - the framework proposed in our paper, the **BR-MAB**.

#### Instructions

