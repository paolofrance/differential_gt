# The differential GT package

This package implements the solution for the two-players Linear-Quadratic differential game-theoretic problem.

This implementation has been used in the works proposed in "Franceschi, Paolo, Nicola Pedrocchi, and Manuel Beschi. [Human–Robot Role Arbitration via Differential Game Theory.](https://ieeexplore.ieee.org/abstract/document/10275780) IEEE Transactions on Automation Science and Engineering (2023).", and in [Modeling and analysis of pHRI with Differential Game Theory](https://arxiv.org/abs/2307.10739).

This repository addressed the solution of the Non-Cooperative case, the so-called Nash Equilibrium, and the solution of the Cooperative case, the Pareto frontier.


Namely, it addresses the NC problem described as 

$$\begin{split}
    & \min_{u_1}J_{1} (z,u_1,u_2) \\
    & \min_{u_2}J_{2} (z,u_1,u_2) \\
    s.t. & \; \dot{z} = Az+B_1 \, u_1 + B_2\, u_2 \\
    & z(t_0) =z_0
  \end{split}$$

and the Cooperative 

$$\begin{split}
    & \min_{u}J_{c}  =  \int_{0}^{\infty}( z^T\,Q_{c}\,z + u^T\,R_{c}\,u )\,dt \\
    & s.t. \; \dot{z} = A_{c}\,z_{c}+B_{c}u_{c} \\
    & z_c(t_0) =z_0
  \end{split}$$

You can find more details in the related publications [here]((https://ieeexplore.ieee.org/abstract/document/10275780)) and [here](https://arxiv.org/abs/2307.10739), and in [Engwerda, Jacob. LQ dynamic optimization and differential games. John Wiley & Sons, 2005.](https://www.wiley.com/en-us/LQ+Dynamic+Optimization+and+Differential+Games-p-9780470015247).


## Citation
If you are using this repository, please consider citing it as 
```
@ARTICLE{10275780,
author={Franceschi, Paolo and Pedrocchi, Nicola and Beschi, Manuel},
journal={IEEE Transactions on Automation Science and Engineering}, 
title={Human–Robot Role Arbitration via Differential Game Theory}, 
year={2023},
volume={},
number={},
pages={1-16},
doi={10.1109/TASE.2023.3320708}}
```

### TODO
1. Better documentation
2. purge from ros



## Contacts
if you have any questions about this repository, are interested in applications, want to contribute or anything else,
please contact me at paolo.franceschi@supsi.ch .

### Acnowledgements
This repository was developed within the [EU Horizon 2020 DrapeBot project](https://www.drapebot.eu/) context. 




