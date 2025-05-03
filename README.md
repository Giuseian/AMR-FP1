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

