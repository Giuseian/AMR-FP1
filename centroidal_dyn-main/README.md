# AMR-FP1

## Project Description




## Repository Structure 
```
AMR-FP1-centroidal_dyn-main/
├── meshes/                            # .dae files defining robot visualization
├── urdf/                              # .urdf files defining ground and robot  
├── report/                            # LaTeX report files
├── outputs/                           # .txt files containing TO solutions 
├── plots/                             # .png files containing TO solutions
├── foot_trajectory_generatory.py      #
├── footstep_planner.py                # 
├── graphics.py                        # Visualization plots 
├── main.py                            # Run the code 
├── sbcdyn.py                          # Stiffness Based Centroidal Dynamics Formulas Definition  
├── simulation.py                      # Visualization of Robot TO Solution  
├── trajectory_generation.py           # Reference Trajectory Generator
└── utils.py                           # Utils file  
```




## Proposed Method


## Results 


# Installation
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
@ARTICLE{10669176,
  author={Tazaki, Yuichi},
  journal={IEEE Robotics and Automation Letters}, 
  title={Trajectory Generation for Legged Robots Based on a Closed-Form Solution of Centroidal Dynamics}, 
  year={2024},
  volume={9},
  number={11},
  pages={9239-9246},
  keywords={Trajectory;Closed-form solutions;Mathematical models;Dynamics;Legged locomotion;Reduced order systems;Closed-form solutions;Reduced order systems;Robot motion;Centroidal dynamics;closed-form solution;legged robots;trajectory generation},
  doi={10.1109/LRA.2024.3455944}}
```

* Other Related Works : 

```bib
@ARTICLE{8955951,
  author={Scianca, Nicola and De Simone, Daniele and Lanari, Leonardo and Oriolo, Giuseppe},
  journal={IEEE Transactions on Robotics}, 
  title={MPC for Humanoid Gait Generation: Stability and Feasibility}, 
  year={2020},
  volume={36},
  number={4},
  pages={1171-1188},
  keywords={Trajectory;Humanoid robots;Stability criteria;Timing;Predictive control;Gait generation;humanoid robots;internal stability;legged locomotion;predictive control;recursive feasibility},
  doi={10.1109/TRO.2019.2958483}}
```

```bib
@article{CIPRIANO2023104495,
title = {Humanoid motion generation in a world of stairs},
journal = {Robotics and Autonomous Systems},
volume = {168},
pages = {104495},
year = {2023},
issn = {0921-8890},
doi = {https://doi.org/10.1016/j.robot.2023.104495},
url = {https://www.sciencedirect.com/science/article/pii/S0921889023001343},
author = {Michele Cipriano and Paolo Ferrari and Nicola Scianca and Leonardo Lanari and Giuseppe Oriolo},
keywords = {Humanoid robot, Footstep Planning, Gait Generation, MPC, Uneven ground, Sensor-based},
abstract = {Consider the problem of generating humanoid motions in an environment consisting of horizontal patches located at different heights (world of stairs). To this end, the paper proposes an integrated scheme which combines footstep planning and gait generation. In particular, footsteps are produced by a randomized algorithm that guarantees both feasibility and quality of the plan according to a chosen criterion; whereas for 3D gait generation we devise an ad hoc extension of the Intrinsically Stable MPC scheme. In its basic form, the proposed scheme addresses the off-line case (known environments), but a sensor-based adaptation is developed for the on-line case (unknown environments) based on an anytime version of the footstep planner. In order to validate the proposed approach, we present simulations in CoppeliaSim for the HRP-4 humanoid robot navigating scenarios of different complexity, both in the on-line and off-line case.}
}
```

