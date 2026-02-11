# Trajectory Optimization for Humanoids via Centroidal Dynamics 

<p align="center">
  <img src="./walking.gif" width="400" alt="Robot walking"/>
</p>

## Project Description
This repository contains the implementation of a trajectory optimization framework for humanoid robots based on **Stiffness-Based Centroidal Dynamics (SBCD)**.
The project builds upon the work by Tazaki et al. (2024) and aims to reproduce and evaluate closed-form centroidal dynamics for **standing** and **walking** tasks.

The core idea is to parametrize contact wrenches through stiffness-like variables, allowing:
- closed-form integration of centroidal dynamics,
- long-horizon trajectory optimization,
- reduced computational complexity compared to full-body models.

The framework is implemented in Python and validated in simulation using the HRP-4 humanoid robot model.



## Repository Structure 
```
AMR-FP1-centroidal_dyn-main/
├── meshes/                            # .dae files defining robot visualization
├── urdf/                              # .urdf files defining ground and robot  
├── report/                            # LaTeX report files
├── outputs/                           # .txt files containing TO solutions 
├── plots/                             # .png files containing TO solutions
├── foot_trajectory_generatory.py      # Foot trajectory reference generator
├── footstep_planner.py                # Contact sequence and step planning
├── graphics.py                        # Visualization plots 
├── main.py                            # Main  
├── sbcdyn.py                          # Stiffness Based Centroidal Dynamics Formulas Definition  
├── simulation.py                      # Visualization of Robot TO Solution  
├── trajectory_generation.py           # Reference Trajectory Generator
└── utils.py                           # Utility functions 
```





## Proposed Method
The proposed method is based on **Stiffness-Based Centroidal Dynamics**, a reduced-order model that describes the robot motion through the evolution of its **Center of Mass (CoM)** and **angular momentum**.

At the centroidal level, the robot dynamics are expressed in terms of the total external wrench acting on the system. Contact interactions are modeled through a spring-like parametrization, where each contact is associated with a stiffness parameter, a CMP-like virtual pivot point, and an optional pure moment term.

Assuming piecewise-constant contact parameters over each contact phase, the resulting dynamics admit **closed-form analytical solutions** for the CoM position, velocity, and angular momentum. These solutions naturally lead to a discrete-time formulation that is well suited for trajectory optimization. The base-link orientation is handled separately through quaternion-based integration, ensuring a consistent representation of rotational motion.



## Mathematical Model (SBCD)
The centroidal dynamics describe the evolution of the CoM position and angular momentum under the effect of gravity and external contact wrenches. At this level, the translational and rotational dynamics are given by:

$$
m\ddot{\mathbf{p}} = \mathbf{f} - m\mathbf{g}
$$

$$
\dot{\mathbf{L}} = \boldsymbol{\eta}
$$

where $\mathbf{p}$ denotes the CoM position, $\mathbf{L}$ the angular momentum about the CoM, $\mathbf{f}$ and $\boldsymbol{\eta}$ the total external force and moment, and $m$ the robot mass.

Each contact wrench is modeled through a stiffness-based parametrization. In this formulation, contact forces behave as virtual springs pulling the CoM toward a pivot point, while contact moments are scaled consistently by the same stiffness parameter:

$$
\mathbf{f}_l = m\lambda_l^2 \left( \mathbf{p} - (\mathbf{p}_l + \mathbf{r}_l) \right), \quad
\boldsymbol{\eta}_l = m\lambda_l^2 \hat{\boldsymbol{\eta}}_l
$$

By aggregating all active contacts, the centroidal dynamics can be rewritten in a compact form, leading to the **Stiffness-Based Centroidal Dynamics (SBCD)** equations. Under the assumption of piecewise-constant parameters, these equations admit closed-form solutions over each contact phase.

In particular, the CoM trajectory over a time interval $[t_k, t_{k+1}]$ is given by:

In particular, the CoM trajectory over a time interval $[t_k, t_{k+1}]$ is given by:

$$
\mathbf{p}(t)
=
(\bar{\mathbf{p}}_k + \bar{\mathbf{r}}_k)
+ C_k(\Delta t)\left(\mathbf{p}_k - (\bar{\mathbf{p}}_k + \bar{\mathbf{r}}_k)\right)
+ \frac{S_k(\Delta t)}{\bar{\lambda}_k}\mathbf{v}_k
$$

with

$$
C_k(\Delta t) = \cosh(\bar{\lambda}_k \Delta t),
\qquad
S_k(\Delta t) = \sinh(\bar{\lambda}_k \Delta t).
$$




## Trajectory Optimization
The trajectory planning problem is formulated as a **finite-horizon optimal control problem**. The system state includes the CoM position and velocity, the base orientation and angular momentum, the positions and orientations of the end-effectors, as well as the time variable used to parametrize phase durations. Control inputs consist of phase durations, end-effector linear and angular velocities, contact stiffness parameters, and CMP offsets.

The objective function combines multiple contributions. A task-related term encourages tracking of reference trajectories for the CoM, feet, and base motion. Limit-related terms enforce physical feasibility through inequality constraints such as friction cones, center-of-pressure bounds, stiffness limits, and phase duration bounds. Additionally, contact-dependent costs ensure consistency between contact states, velocities, and stiffness values.

In compact form, the optimization problem can be written as:

$$
\min_{\mathbf{x}, \mathbf{u}}
\sum_k
\left(
\mathcal{L}_{\text{task},k}
+ \mathcal{L}_{\text{limit},k}
+ \mathcal{L}_{\text{contact},k}
\right)
\quad \text{s.t.} \quad
\mathbf{x}_{k+1} = f(\mathbf{x}_k, \mathbf{u}_k)
$$


The problem is solved using **CasADi** and the **IPOPT** nonlinear programming solver.



<!-- ## Tasks Implemented
Two representative tasks are implemented to validate the proposed framework.

**Still Task.**  
Both feet remain in contact with the ground while the robot maintains a stable posture. The objective is to compensate gravity and preserve balance, allowing the evaluation of force distribution and centroidal stability.

**Walking Task.**  
The robot performs a walking motion characterized by alternating single-support and double-support phases. A predefined contact sequence and footstep plan are used, while CoM and foot trajectories are optimized to generate forward locomotion with lateral sway and vertical motion. -->


## Tasks Implemented & Results 
Two representative tasks are implemented to validate the proposed framework.

| Task    | Description |
|--------|-------------|
| Still | Static balance with both feet in contact, used to validate force distribution and stability. |
| Walking | Dynamic locomotion with alternating support phases and optimized CoM and foot trajectories. |


Simulation results show that the robot successfully maintains balance during the still task, with stable CoM motion and bounded angular momentum. During walking, the optimized trajectories reproduce realistic CoM evolution, foot swing phases, and physically consistent contact forces. In all cases, friction and moment constraints are satisfied, and the resulting motions are dynamically feasible while maintaining a relatively low computational cost.

Plots and numerical results are stored in the `plots/` and `outputs/` folders.


## Installation
You need a Python installation and some dependencis. If using PIP, you can run the following
```
pip install dartpy casadi scipy matplotlib osqp
```
You need dartpy 0.2 : If pip does not allow you to install this version on your system, you might want to use conda

To run the simulation, choose one of the following tasks 
- **still**, to simulate a fixed robot
- **walking**, to simulate a walking robot
```
python main.py [task] --make_video
```
if you want to save the video (more time)
or
```
python main.py [task] 
```

Then, to see the robot simulation in Dartpy, run
```
python simulation.py
```


## Acknowledgments
* The original paper:

```bib
@article{tazaki2024trajectory,
  title={Trajectory generation for legged robots based on a closed-form solution of centroidal dynamics},
  author={Tazaki, Yuichi},
  journal={IEEE Robotics and Automation Letters},
  year={2024},
  publisher={IEEE}
}
```

* Other Related Works : 

```bib
@article{scianca2020mpc,
  title={MPC for humanoid gait generation: Stability and feasibility},
  author={Scianca, Nicola and De Simone, Daniele and Lanari, Leonardo and Oriolo, Giuseppe},
  journal={IEEE Transactions on Robotics},
  volume={36},
  number={4},
  pages={1171--1188},
  year={2020},
  publisher={IEEE}
}
```

```bib
@article{cipriano2023humanoid,
  title={Humanoid motion generation in a world of stairs},
  author={Cipriano, Michele and Ferrari, Paolo and Scianca, Nicola and Lanari, Leonardo and Oriolo, Giuseppe},
  journal={Robotics and Autonomous Systems},
  volume={168},
  pages={104495},
  year={2023},
  publisher={Elsevier}
}
```

