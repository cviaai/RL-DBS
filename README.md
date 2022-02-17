[![License](https://img.shields.io/github/license/analysiscenter/pydens.svg)](https://www.apache.org/licenses/LICENSE-2.0)
[![Python](https://img.shields.io/badge/python-3.6-blue.svg)](https://python.org)

# RL-DBS
# Reinforcement learning for suppression of collective neuronal activity in Deep Brain Stimulation (DBS)

This is a convenient gym environment for developing and comparing interaction of RL agents with several types of synthetic neuronal models of pathological brain activity. The ODEs that simulate synchronized neuronal signaling are wrapped into the the framework as individual environments, allowing simple switching between different physical models, and enabling convenient approbation of various RL approaches and multiple agents. The policy gradient algorithm PPO is shown to provide a robust data-driven control of neronal synchoryzation, agnostic of the neuronal model. 

<p align="center">
<img src="RL-DBS-diagram.png" width="600" alt>
</p>
<p align="center">
<em>Principle diagram of Reinforcement Learning for Deep Brain Stimulation systems.</em>
</p>

We propose a class of physically meaningful reward functions enabling the suppression of collective oscillatory mode. The synchrony suppression is demonstrated for two models of neuronal populations – for the ensembles of globally coupled limit-cycle Bonhoeffer-van der Pol oscillators and for the bursting Hindmarsh–Rose neurons. The suppression workflow proposed here is universal and could be used to create benchmarks among different physical models, to create different control algorithms, and to pave the way towards clinical realization of deep brain stimulation via reinforcement learning. 

<p align="center">
<!---<img src="RL-DBS-demo.png" alt> -->
<img src="2.gif" alt>

</p>
<p align="center">
<em>Demonstration of synchrony suppression via PPO A2C algorithm.</em>
</p>




### Installation as a project repository:

```
git clone https://github.com/cviaai/RL-DBS.git
```

In this case, you need to manually install the dependencies.

### Important notes:

Notebook file Baseline shows how the model interacts with the environment and you can find an example of a trained model in this notebook. Separate training from scratch is also shown.

Environment uses generic Gym notation. A class that describes all relevant information is:
```
gym_oscillator/envs/osc_env.py
```
A C++ code that emulates neuronal oscillations is the following file:
```
/source_for_build_files/gfg.c
```
## Citing 

If you use this package in your publications or in other work, please cite these papers:


@article{krylov-RL-DBS,

title = {Reinforcement Learning Framework for Deep Brain Stimulation Study},
author = {Krylov, Dmitrii and des Combes, Remi and Laroche, Romain and Rosenblum, Michael and Dylov, Dmitry V},
journal = {Proceedings of the Twenty-Ninth International Joint Conference on Artificial Intelligence, {IJCAI-20}},
pages = {2847--2854},
doi = {10.24963/ijcai.2020/394},
url = {https://www.ijcai.org/proceedings/2020/0394},
year = {2020}

}

@article{krylov-RL-Chaos,

    title={Reinforcement Learning for Suppression of Collective Activity in Oscillatory Ensembles},
    author={Dmitrii Krylov and Dmitry V. Dylov and Michael Rosenblum},
    journal = {Chaos: An Interdisciplinary Journal of Nonlinear Science},
    number = {3},
    pages = {33126},
    volume = {30},
    doi = {10.1063/1.5128909},
    url = {https://aip.scitation.org/doi/10.1063/1.5128909}
    year = {2020}
}
