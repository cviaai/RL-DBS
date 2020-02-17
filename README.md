[![License](https://img.shields.io/github/license/analysiscenter/pydens.svg)](https://www.apache.org/licenses/LICENSE-2.0)
[![Python](https://img.shields.io/badge/python-3.6-blue.svg)](https://python.org)

# RL-suppression-in-oscillatory-ensembles
Reinforcement learning for suppression of collective activity in oscillatory ensembles

The hybrid model relies on two major components: an environment of oscillators and a policy-based reinforcement learning block. This repository features a model-agnostic synchrony control based on proximal policy optimization and two artificial neural networks in an Actor-Critic configuration. 

A class of physically meaningful reward functions enabling the suppression of collective oscillatory mode is proposed. The synchrony suppression is demonstrated for two models of neuronal populations – for the ensembles of globally coupled limit-cycle Bonhoeffer-van der Pol oscillators and for the bursting Hindmarsh–Rose neurons.

<p align="center">
<img src="principle.png" alt>
</p>
<p align="center">
<em>Principle diagram of Reinforcement Learning via PPO Actor-Critic algorithm.</em>
</p>



### Installation as a project repository:

```
git clone https://github.com/cviaai/RL-suppression-in-oscillatory-ensembles.git
```

In this case you need to manually install the dependencies.

### There are few important notes:

Notebook file Baseline shows how model interacts with environment and in this notebook you can find an example of trained model and training.

Environment uses Gym Notation, class that describes all information is 
```
gym_oscillator/envs/osc_env.py
```
A C++ code that creates oscillations is the following:
```
/source_for_build_files/gfg.c
```
## Citing 

If you use this package in your publications or other work please cite the package as follows:

```
Dmitriy Krylov, Dmitry V. Dylov, Michael Rosenblum. 2019.
```

```
@misc{RLsuppr_2019,
  author       = {Krylov D. and Dylov D. and Rosenblum M.},
  title        = {Reinforcement learning for suppression of collective activity in oscillatory ensembles},
  year         = 2019
}
```
