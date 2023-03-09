![PlatoSim logo](docs/figures/LogoPlatoSim.png "PlatoSim logo")

# PlatoSim: The end-to-end PLATO camera simulator
[![Image](https://img.shields.io/badge/license-MIT-blue.svg "")](https://github.com/IvS-KULeuven/PlatoSim3/blob/master/LICENSE.txt)
[![Image](https://img.shields.io/badge/documentation-%E2%9C%93-blue.svg)](https://ivs-kuleuven.github.io/PlatoSim3/)
[![Image](https://img.shields.io/badge/tutorials-%E2%9C%93-blue.svg)](https://github.com/IvS-KULeuven/PlatoSim3/tree/master/docs/tutorials)
[![Image](https://img.shields.io/badge/c++-00599C?style=flat&logo=c%2B%2B&logoColor=white)](https://ivs-kuleuven.github.io/PlatoSim3/_simulation_steps.html)
[![Image](https://img.shields.io/badge/Python-3766AB?style=flat&logo=Python&logoColor=white)](https://ivs-kuleuven.github.io/PlatoSim3/user-overview.html)
[![Image](https://img.shields.io/badge/Jenkins-D24939?style=flat&logo=Jenkins&logoColor=white)](https://ivs-kuleuven.github.io/PlatoSim3/user-jenkins.html)
<!-- [![Image](https://img.shields.io/badge/Docker-2496ED?style=flat-square&logo=Docker&logoColor=white)]() -->
<!-- [![Image](https://img.shields.io/badge/pip%20install-transitleastsquares-blue.svg)]() -->
<!-- [![Image](https://img.shields.io/badge/arXiv-1901.02015-blue.svg)]() -->

## Motivation

To accommodate PLATO's need of versatile simulations prior to mission launch - that at the same time describe accurately the innovative but complex multi-telescope design - we here present Platosim, the end-to-end PLATO camera simulator specifically developed for purpose. PlatoSim allows the user to simulate photometric time series of CCD images and light curves in accordance to the expected observations of PLATO. In the context of the PLATO payload, PlatoSim uses a general formalism of modelling the stellar field and sky background, the short and long-term barycentric pixel displacement of the stellar sources, the cameras and their optics, the CCDs and their electronics, and all main random and systematic noise sources. With its strong predictive powers and diverse applicability, PlatoSim is key simulator for PLATO Mission Consortium.

## Access and contribution

We welcome new users and encourage contributions. Simply contact one of the PlatoSim developers from KU Leuven and we will make sure to give you access. Note that in agreement with the regulations of the PLATO mission, you need a *Non Disclosure Agreement* (NDA) before we can give you access to PlatoSim.

## Installation

PlatoSim has the following installation procedures:

* [conda installation procedure](https://ivs-kuleuven.github.io/PlatoSim3/user-install.html) (recommended)
* [straight from GitHub](https://ivs-kuleuven.github.io/PlatoSim3/dev-fork-clone.html)

## Getting Started

In order to provide a smooth start of your PlatoSim journey, we suggest that you both consult our <br> 
[Documentation Page](https://ivs-kuleuven.github.io/PlatoSim3/) and our [Python Turorials](https://github.com/IvS-KULeuven/PlatoSim3/tree/master/docs/tutorials).

## Reference

Please cite *Jannsen et al. (in prep.)* if use PlatoSim in your research. The ADS BibTex entry for this paper will be available soon.

<!-- The BibTeX entry for the paper is: -->

## Feedback and Issues

If you have any questions or experience any trouble while using PlatoSim, please open a [GitHub Issue](https://github.com/IvS-KULeuven/PlatoSim3/issues). <br> 
In case of an issue -- to help us help you -- we recommend to provide the following information:

* PlatoSim version (bash command: `platosim --version`)
* Appropiate issue *label* (see right-hand menu)
* Concise explanation of the problem
* What you tried to investigate the problem
* Please provide the `inputfile.yaml` and `log.txt` files (in debug mode: `--verbosity 3`)

---

Copyright 2023 - KU Leuven & The PlatoSim team
