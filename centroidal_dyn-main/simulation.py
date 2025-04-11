import numpy as np
import dartpy as dart
import copy
from utils import *
import os
import ismpc
import footstep_planner
import inverse_dynamics as id
import filter
import foot_trajectory_generator as ftg
from logger import Logger

class Hrp4Controller(dart.gui.osg.RealTimeWorldNode):
    def __init__(self, world, hrp4):
        super(Hrp4Controller, self).__init__(world)
        self.world = world
        self.hrp4 = hrp4
        self.time = 0
        self.params = {
            'g': 9.81,
            'h': 0.72,
            'foot_size': 0.1,
            'step_height': 0.02,
            'ss_duration': 70,
            'ds_duration': 30,
            'world_time_step': world.getTimeStep(),
            'first_swing': 'rfoot',
            'µ': 0.5,
            'N': 100,
            'dof': self.hrp4.getNumDofs(),
        }
        self.params['eta'] = np.sqrt(self.params['g'] / self.params['h'])

        # robot links
        self.lsole = hrp4.getBodyNode('l_sole')
        self.rsole = hrp4.getBodyNode('r_sole')
        self.torso = hrp4.getBodyNode('torso')
        self.base  = hrp4.getBodyNode('body')

        for i in range(hrp4.getNumJoints()):
            joint = hrp4.getJoint(i)
            dim = joint.getNumDofs()

            # set floating base to passive, everything else to torque
            if   dim == 6: joint.setActuatorType(dart.dynamics.ActuatorType.PASSIVE)
            elif dim == 1: joint.setActuatorType(dart.dynamics.ActuatorType.FORCE)

        # set initial configuration
        initial_configuration = {'CHEST_P': 0., 'CHEST_Y': 0., 'NECK_P': 0., 'NECK_Y': 0., \
                                 'R_HIP_Y': 0., 'R_HIP_R': -3., 'R_HIP_P': -25., 'R_KNEE_P': 50., 'R_ANKLE_P': -25., 'R_ANKLE_R':  3., \
                                 'L_HIP_Y': 0., 'L_HIP_R':  3., 'L_HIP_P': -25., 'L_KNEE_P': 50., 'L_ANKLE_P': -25., 'L_ANKLE_R': -3., \
                                 'R_SHOULDER_P': 4., 'R_SHOULDER_R': -8., 'R_SHOULDER_Y': 0., 'R_ELBOW_P': -25., \
                                 'L_SHOULDER_P': 4., 'L_SHOULDER_R':  8., 'L_SHOULDER_Y': 0., 'L_ELBOW_P': -25.}

        for joint_name, value in initial_configuration.items():
            self.hrp4.setPosition(self.hrp4.getDof(joint_name).getIndexInSkeleton(), value * np.pi / 180.)

        # position the robot on the ground
        lsole_pos = self.lsole.getTransform(withRespectTo=dart.dynamics.Frame.World(), inCoordinatesOf=dart.dynamics.Frame.World()).translation()
        rsole_pos = self.rsole.getTransform(withRespectTo=dart.dynamics.Frame.World(), inCoordinatesOf=dart.dynamics.Frame.World()).translation()
        self.hrp4.setPosition(3, - (lsole_pos[0] + rsole_pos[0]) / 2.)
        self.hrp4.setPosition(4, - (lsole_pos[1] + rsole_pos[1]) / 2.)
        self.hrp4.setPosition(5, - (lsole_pos[2] + rsole_pos[2]) / 2.)

        # initialize state
        self.initial = self.retrieve_state()
        self.contact = 'lfoot' if self.params['first_swing'] == 'rfoot' else 'rfoot' # there is a dummy footstep
        self.desired = copy.deepcopy(self.initial)

        # selection matrix for redundant dofs
        redundant_dofs = [ \
            "NECK_Y", "NECK_P", \
            "R_SHOULDER_P", "R_SHOULDER_R", "R_SHOULDER_Y", "R_ELBOW_P", \
            "L_SHOULDER_P", "L_SHOULDER_R", "L_SHOULDER_Y", "L_ELBOW_P"]
        
        # initialize inverse dynamics
        self.id = id.InverseDynamics(self.hrp4, redundant_dofs)

        # initialize footstep planner
        reference = [(0.1, 0., 0.2)] * 5 + [(0.1, 0., -0.1)] * 10 + [(0.1, 0., 0.)] * 10
        self.footstep_planner = footstep_planner.FootstepPlanner(
            reference,
            self.initial['lfoot']['pos'],
            self.initial['rfoot']['pos'],
            self.params
            )

        # initialize MPC controller
        self.mpc = ismpc.Ismpc(
            self.initial, 
            self.footstep_planner, 
            self.params
            )

        # initialize foot trajectory generator
        self.foot_trajectory_generator = ftg.FootTrajectoryGenerator(
            self.initial, 
            self.footstep_planner, 
            self.params
            )


        # initialize kalman filter
        A = np.identity(3) + self.params['world_time_step'] * self.mpc.A_lip
        B = self.params['world_time_step'] * self.mpc.B_lip
        d = np.zeros(9)
        d[7] = - self.params['world_time_step'] * self.params['g']
        H = np.identity(3)
        Q = block_diag(1., 1., 1.)
        R = block_diag(1e1, 1e2, 1e4)
        P = np.identity(3)
        x = np.array([self.initial['com']['pos'][0], self.initial['com']['vel'][0], self.initial['zmp']['pos'][0], \
                      self.initial['com']['pos'][1], self.initial['com']['vel'][1], self.initial['zmp']['pos'][1], \
                      self.initial['com']['pos'][2], self.initial['com']['vel'][2], self.initial['zmp']['pos'][2]])
        self.kf = filter.KalmanFilter(block_diag(A, A, A), \
                                      block_diag(B, B, B), \
                                      d, \
                                      block_diag(H, H, H), \
                                      block_diag(Q, Q, Q), \
                                      block_diag(R, R, R), \
                                      block_diag(P, P, P), \
                                      x)

        # initialize logger and plots
        self.logger = Logger(self.initial)
        self.logger.initialize_plot(frequency=10)
        
    def customPreStep(self):
        
        # create current and desired states
        self.current = self.retrieve_state()

        # here our code

        phases_duration = np.loadtxt('./outputs/phase_durations.txt', delimiter=',')
        desired_positions, desired_velocities = [], []
        contacts = ["ds", "rfoot"] # change this wrt sigma sequence

        
        pos = np.loadtxt("./outputs/com_pos.txt", delimiter=',').tolist()
        vel = np.loadtxt("./outputs/com_vel.txt", delimiter=',').tolist()
        
        for i in range(len(pos)): 
            desired_positions.append(pos[i])
            desired_velocities.append(vel[i])
        
        phase1 = phases_duration[0]*100
        fin_phase2 = (phases_duration[0]+phases_duration[1])*100 
        if self.time >= 0 and self.time < phase1:
            contact = contacts[0]
        elif self.time >= phase1 and self.time < fin_phase2:
            contact = contacts[1]
                
        if self.time < fin_phase2 - 2:
            print(f"time: {self.time}, contact: {contact}") 
            
            desired_pos = desired_positions[self.time]
            
            desired_vel = desired_velocities[self.time]                
            self.desired['com']['pos'] = desired_pos
            self.desired['com']['vel'] = desired_vel
            if contact != "ds":
                acc_x = desired_vel[0] / (self.time) 
                
            else:
                acc_x = 0.0

            self.desired['com']['acc'] = [acc_x, 0.0, 0.0] 
            print(f"COM pos: {desired_pos}")
            print(f"COM vel: {desired_vel}")
            print(f"COM x acc: {acc_x}")
               
            commands = self.id.get_joint_torques(self.desired, self.current, contact)    
            for i in range(self.params['dof'] - 6):
                self.hrp4.setCommand(i + 6, commands[i])

            self.time += 1
        else:
            exit()
        
        # get foot trajectories
        feet_trajectories = self.foot_trajectory_generator.generate_feet_trajectories_at_time(self.time)
        for foot in ['lfoot', 'rfoot']:
            for key in ['pos', 'vel', 'acc']:
                self.desired[foot][key] = feet_trajectories[foot][key]

        # set torso and base references to the average of the feet
        for link in ['torso', 'base']:
            for key in ['pos', 'vel', 'acc']:
                self.desired[link][key] = (self.desired['lfoot'][key][:3] + self.desired['rfoot'][key][:3]) / 2.
        
        
        # log and plot
        self.logger.log_data(self.current, self.desired)
        #self.logger.update_plot(self.time)

        
    def retrieve_state(self):
        # com and torso pose (orientation and position)
        com_position = self.hrp4.getCOM()
        torso_orientation = get_rotvec(self.hrp4.getBodyNode('torso').getTransform(withRespectTo=dart.dynamics.Frame.World(), inCoordinatesOf=dart.dynamics.Frame.World()).rotation())
        base_orientation  = get_rotvec(self.hrp4.getBodyNode('body' ).getTransform(withRespectTo=dart.dynamics.Frame.World(), inCoordinatesOf=dart.dynamics.Frame.World()).rotation())

        # feet poses (orientation and position)
        l_foot_transform = self.lsole.getTransform(withRespectTo=dart.dynamics.Frame.World(), inCoordinatesOf=dart.dynamics.Frame.World())
        l_foot_orientation = get_rotvec(l_foot_transform.rotation())
        l_foot_position = l_foot_transform.translation()
        left_foot_pose = np.hstack((l_foot_orientation, l_foot_position))
        r_foot_transform = self.rsole.getTransform(withRespectTo=dart.dynamics.Frame.World(), inCoordinatesOf=dart.dynamics.Frame.World())
        r_foot_orientation = get_rotvec(r_foot_transform.rotation())
        r_foot_position = r_foot_transform.translation()
        right_foot_pose = np.hstack((r_foot_orientation, r_foot_position))

        # velocities
        com_velocity = self.hrp4.getCOMLinearVelocity(relativeTo=dart.dynamics.Frame.World(), inCoordinatesOf=dart.dynamics.Frame.World())
        torso_angular_velocity = self.hrp4.getBodyNode('torso').getAngularVelocity(relativeTo=dart.dynamics.Frame.World(), inCoordinatesOf=dart.dynamics.Frame.World())
        base_angular_velocity = self.hrp4.getBodyNode('body').getAngularVelocity(relativeTo=dart.dynamics.Frame.World(), inCoordinatesOf=dart.dynamics.Frame.World())
        l_foot_spatial_velocity = self.lsole.getSpatialVelocity(relativeTo=dart.dynamics.Frame.World(), inCoordinatesOf=dart.dynamics.Frame.World())
        r_foot_spatial_velocity = self.rsole.getSpatialVelocity(relativeTo=dart.dynamics.Frame.World(), inCoordinatesOf=dart.dynamics.Frame.World())

        # compute total contact force
        force = np.zeros(3)
        for contact in world.getLastCollisionResult().getContacts():
            force += contact.force

        # compute zmp
        zmp = np.zeros(3)
        zmp[2] = com_position[2] - force[2] / (self.hrp4.getMass() * self.params['g'] / self.params['h'])
        for contact in world.getLastCollisionResult().getContacts():
            if contact.force[2] <= 0.1: continue
            zmp[0] += (contact.point[0] * contact.force[2] / force[2] + (zmp[2] - contact.point[2]) * contact.force[0] / force[2])
            zmp[1] += (contact.point[1] * contact.force[2] / force[2] + (zmp[2] - contact.point[2]) * contact.force[1] / force[2])

        if force[2] <= 0.1: # threshold for when we lose contact
            zmp = np.array([0., 0., 0.]) # FIXME: this should return previous measurement
        else:
            # sometimes we get contact points that dont make sense, so we clip the ZMP close to the robot
            midpoint = (l_foot_position + l_foot_position) / 2.
            zmp[0] = np.clip(zmp[0], midpoint[0] - 0.3, midpoint[0] + 0.3)
            zmp[1] = np.clip(zmp[1], midpoint[1] - 0.3, midpoint[1] + 0.3)
            zmp[2] = np.clip(zmp[2], midpoint[2] - 0.3, midpoint[2] + 0.3)

        # create state dict
        return {
            'lfoot': {'pos': left_foot_pose,
                      'vel': l_foot_spatial_velocity,
                      'acc': np.zeros(6)},
            'rfoot': {'pos': right_foot_pose,
                      'vel': r_foot_spatial_velocity,
                      'acc': np.zeros(6)},
            'com'  : {'pos': com_position,
                      'vel': com_velocity,
                      'acc': np.zeros(3)},
            'torso': {'pos': torso_orientation,
                      'vel': torso_angular_velocity,
                      'acc': np.zeros(3)},
            'base' : {'pos': base_orientation,
                      'vel': base_angular_velocity,
                      'acc': np.zeros(3)},
            'joint': {'pos': self.hrp4.getPositions(),
                      'vel': self.hrp4.getVelocities(),
                      'acc': np.zeros(self.params['dof'])},
            'zmp'  : {'pos': zmp,
                      'vel': np.zeros(3),
                      'acc': np.zeros(3)}
        }

import os
import dartpy as dart
import numpy as np

if __name__ == "__main__":
    
    # Initialize world
    world = dart.simulation.World()
    
    # URDF loader
    urdfParser = dart.utils.DartLoader()
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Paths to URDF files
    hrp4_path = os.path.join(current_dir, "urdf", "hrp4.urdf")
    ground_path = os.path.join(current_dir, "urdf", "ground.urdf")
    
    # Check if files exist
    assert os.path.exists(hrp4_path), f"URDF file not found: {hrp4_path}"
    assert os.path.exists(ground_path), f"URDF file not found: {ground_path}"
    
    # Parse URDF files
    hrp4 = urdfParser.parseSkeleton(hrp4_path)
    if hrp4 is None:
        raise RuntimeError(f"Failed to parse HRP4 URDF: {hrp4_path}")
    
    ground = urdfParser.parseSkeleton(ground_path)
    if ground is None:
        raise RuntimeError(f"Failed to parse Ground URDF: {ground_path}")
    
    # Add skeletons to the world
    world.addSkeleton(hrp4)
    world.addSkeleton(ground)
    world.setGravity([0, 0, -9.81])
    world.setTimeStep(0.01)
    
    # Set default inertia for bodies with zero mass
    default_inertia = dart.dynamics.Inertia(1e-8, np.zeros(3), 1e-10 * np.identity(3))
    for body in hrp4.getBodyNodes():
        if body.getMass() == 0.0:
            print(f"Setting default mass and inertia for body: {body.getName()}")
            body.setMass(1e-8)
            body.setInertia(default_inertia)
    
    # Initialize HRP4 controller
    try:
        node = Hrp4Controller(world, hrp4)
    except NameError:
        raise RuntimeError("Hrp4Controller is not defined. Ensure it is implemented and imported.")
    
    # Create viewer
    viewer = dart.gui.osg.Viewer()
    node.setTargetRealTimeFactor(10)  # Speed up visualization by 10x
    viewer.addWorldNode(node)
    
    # Viewer settings
    viewer.setUpViewInWindow(0, 0, 1280, 720)
    viewer.setCameraHomePosition([5., -1., 1.5],
                                 [1., 0., 0.5],
                                 [0., 0., 1.])
    
    # Run the viewer
    viewer.run()


if __name__ == "__main__":
    world = dart.simulation.World()

    urdfParser = dart.utils.DartLoader()
    current_dir = os.path.dirname(os.path.abspath(__file__))
    hrp4 = urdfParser.parseSkeleton(os.path.join(current_dir, "urdf", "hrp4.urdf"))
    ground = urdfParser.parseSkeleton(os.path.join(current_dir, "urdf", "ground.urdf"))
    world.addSkeleton(hrp4)
    world.addSkeleton(ground)
    world.setGravity([0, 0, -9.81])
    world.setTimeStep(0.01)

    # set default inertia
    default_inertia = dart.dynamics.Inertia(1e-8, np.zeros(3), 1e-10 * np.identity(3))
    for body in hrp4.getBodyNodes():
        if body.getMass() == 0.0:
            body.setMass(1e-8)
            body.setInertia(default_inertia)
   
    node = Hrp4Controller(world, hrp4)

    # create world node and add it to viewer
    viewer = dart.gui.osg.Viewer()
    node.setTargetRealTimeFactor(10) # speed up the visualization by 10x
    viewer.addWorldNode(node)

    #viewer.setUpViewInWindow(0, 0, 1920, 1080)
    viewer.setUpViewInWindow(0, 0, 1280, 720)
    #viewer.setUpViewInWindow(0, 0, 640, 480)
    viewer.setCameraHomePosition([5., -1., 1.5],
                                 [1.,  0., 0.5],
                                 [0.,  0., 1. ])
    viewer.run()
    
