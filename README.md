# Rényi security framework against coherent attacks applied to decoy-state QKD

This is a public version of the code used in *Rényi security framework against coherent attacks applied to decoy-state QKD* \[[arXiv](https://arxiv.org/abs/2504.12248)\]. This was built on top of [v2.0.2](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/releases/tag/v2.0.2) of the Open QKD Security package.

The necessary functions for each protocol can be found in the corresponding folders contained in `RenyiProtocols`. Each folder contains a `main` file which generates the data of the corresponding figure. Furthermore, every folder contains a `plot` file which creates the plot. In the table below we summarize which files need to be run for each figure.

| Figure                                                                      | Folder                           | main file                            | plot function                    |
| --------------------------------------------------------------------------- | -------------------------------- | ------------------------------------ | -------------------------------- |
| Fig. 1 - Qubit BB84                                                         | `QubitBB84`                      | `mainRenyiQubitBB84Lossy.`           | `plotQubitBB84.m`                |
| Fig. 2 - active decoy BB84 (2 decoy)                                        | `ActiveDecoyBB84/2 Decoy`        | `mainRenyiDecoyBB84_1Decoy.m`        | `plotDecoyBB84_1Decoy.m`         |
| Fig. 3 - active decoy BB84 (1 decoy)                                        | `ActiveDecoyBB84/1 Decoy`        | `mainRenyiDecoyBB84_2Decoy.m`        | `plotDecoyBB84_2Decoy.m`         |
| Fig. 4 - passive decoy BB84 <br> with intensity imperfections (2 decoy)     | `PassiveDecoyBB84/2 Decoy`       | `mainRenyiDecoyBB84Passive_2Decoy.m` | `plotPassiveDecoyBB84_1Decoy.m`  |
| Fig. 5 - passive decoy BB84 <br> with intensity imperfections (1 decoy)     | `PassiveDecoyBB84/1 Decoy`       | `mainRenyiDecoyBB84Passive_1Decoy.m` | `plotPassiveDecoyBB84_2Decoy.m`  |
| Fig. 6 - passive decoy 4-6 <br> w/o phase imperfections (1 decoy), q=1      | `Decoy46/PerfectDecoy46`         | `mainRenyiPerfectDecoy46.m`          | `plot46Decoy.m`                  |
| Fig. 6 - passive decoy 4-6 <br> with phase imperfections (1 decoy), q=0.99  | `Decoy46/PhaseImperfectDecoy46`  | `mainRenyiPhaseImpDecoy46.m`         | `plot46Decoy.m`                  |

Figure 4 until Figure 6 require MOSEK 10.0.44 (or above). For instalation instructions on the newest version of MOSEK, follow \[[MOSEK](https://www.mosek.com/)].

## Installation instructions
> [!CAUTION]
> This repository is for archival and transparency purposes; we do not guarantee compatibility with other versions of the Open QKD Security package beyond the ones listed above.

### As zip
1. Download the linked version of the code from above and follow all [installation instructions](LINK TO COMMIT).
2. Also follow the additional Mosek install instructions if you want an exact match.
3. Download the latest release on the side bar and unzip in your preferred directory and add this folder to the Matlab path.


### with git
1. Clone this repository and its exact submodules navigate to your desired directory and run,
```
git clone --recurse-submodules https://github.com/Optical-Quantum-Communication-Theory/Renyi-security-framework
```
2. Follow all further [installation instructions](LINK TO COMMIT).
3. Also follow the additional Mosek install instructions if you want an exact match.
