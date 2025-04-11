import numpy as np

class FootTrajectoryGenerator:
  def __init__(self, initial, footstep_planner, params):
    self.delta = params['world_time_step']
    self.step_height = params['step_height']
    self.initial = initial
    self.footstep_planner = footstep_planner
    #self.plan = self.footstep_planner.plan
    
    #self.footstep_planner.plan
    
    ref_foot = "./outputs/ref_foot.txt"  
    data = np.loadtxt(ref_foot, delimiter=",")

    # Extract the required values: this is only valid for [0,0,0], [0,-1,0] contact sequence
    # Extract required values
    # [0, 0, 0], [0, -1, 0]
    
    lfoot1 = data[-3:, 0]
    rfoot1 = data[:3,1]
    lfoot2 = data[-3:, 2]
    print(f"first feet (R): {lfoot1}")
    print(f"second feet (L): {rfoot1}")
    print(f"third feet (R): {lfoot2}")
    supports = [lfoot1, rfoot1, lfoot2]
    supports_name = ["lfoot", "rfoot", "lfoot"]
    phases_durations = "./outputs/phase_durations.txt"
    duration = np.loadtxt(phases_durations, delimiter=",")
    ang = np.array([0,0,0])
    
    self.plan = [
                    {
                    'pos'        : lfoot1,
                    'ang'        : ang,
                    'ss_duration': 0,
                    'ds_duration': phases_durations[0], 
                    'foot_id'    : "lfoot"
                    },
                    {
                    'pos'        : rfoot1,
                    'ang'        : ang,
                    'ss_duration': phases_durations[1],
                    'ds_duration': phases_durations[1], 
                    'foot_id'    : "rfoot"
                    },
                    {
                    'pos'        : lfoot2,
                    'ang'        : ang,
                    'ss_duration': phases_durations[1],
                    'ds_duration': phases_durations[1], 
                    'foot_id'    : "lfoot"
                    },
                    {
                    'pos'        : rfoot1,
                    'ang'        : ang,
                    'ss_duration': phases_durations[1],
                    'ds_duration': phases_durations[1], 
                    'foot_id'    : "rfoot"
                    }
    ]
    
  def generate_feet_trajectories_at_time(self, time):
    step_index = self.footstep_planner.get_step_index_at_time(time)
    time_in_step = time - self.footstep_planner.get_start_time(step_index)
    phase = self.footstep_planner.get_phase_at_time(time)
    support_foot = self.footstep_planner.plan[step_index]['foot_id']
    swing_foot = 'lfoot' if support_foot == 'rfoot' else 'rfoot'
    single_support_duration = self.footstep_planner.plan[step_index]['ss_duration']

    # if first step, return initial foot poses with zero velocities and accelerations
    if step_index == 0:
        zero_vel = np.zeros(6)
        zero_acc = np.zeros(6)
        return {
            'lfoot': {
                'pos': self.initial['lfoot']['pos'],
                'vel': zero_vel,
                'acc': zero_acc
            },
            'rfoot': {
                'pos': self.initial['rfoot']['pos'],
                'vel': zero_vel,
                'acc': zero_acc
            }
        }

    # if double support, return planned foot poses with zero velocities and accelerations
    if phase == 'ds':
        support_pose = np.hstack((
            self.plan[step_index]['ang'],
            self.plan[step_index]['pos']
        ))
        swing_pose = np.hstack((
            self.plan[step_index + 1]['ang'],
            self.plan[step_index + 1]['pos']
        ))
        zero_vel = np.zeros(6)
        zero_acc = np.zeros(6)
        return {
            support_foot: {
                'pos': support_pose,
                'vel': zero_vel,
                'acc': zero_acc
            },
            swing_foot: {
                'pos': swing_pose,
                'vel': zero_vel,
                'acc': zero_acc
            }
        }
    
    # get positions and angles for cubic interpolation
    start_pos  = self.plan[step_index - 1]['pos']
    target_pos = self.plan[step_index + 1]['pos']
    start_ang  = self.plan[step_index - 1]['ang']
    target_ang = self.plan[step_index + 1]['ang']

    # time variables
    t = time_in_step
    T = single_support_duration

    # cubic polynomial for position and angle
    A = - 2 / T**3
    B =   3 / T**2
    swing_pos     = start_pos + (target_pos - start_pos) * (    A * t**3 +     B * t**2)
    swing_vel     =             (target_pos - start_pos) * (3 * A * t**2 + 2 * B * t   ) / self.delta
    swing_acc     =             (target_pos - start_pos) * (6 * A * t    + 2 * B       ) / self.delta**2
    swing_ang_pos = start_ang + (target_ang - start_ang) * (    A * t**3 +     B * t**2)
    swing_ang_vel =             (target_ang - start_ang) * (3 * A * t**2 + 2 * B * t   ) / self.delta
    swing_ang_acc =             (target_ang - start_ang) * (6 * A * t    + 2 * B       ) / self.delta**2

    # quartic polynomial for vertical position
    A =   16 * self.step_height / T**4
    B = - 32 * self.step_height / T**3
    C =   16 * self.step_height / T**2
    swing_pos[2] =       A * t**4 +     B * t**3 +     C * t**2
    swing_vel[2] = ( 4 * A * t**3 + 3 * B * t**2 + 2 * C * t   ) / self.delta
    swing_acc[2] = (12 * A * t**2 + 6 * B * t    + 2 * C       ) / self.delta**2

    #print("swing_pos: ", swing_pos)
    #print("swing_vel: ", swing_vel)
    #print("swing_acc: ", swing_acc)
    # support foot remains stationary
    support_pos = self.plan[step_index]['pos']
    support_ang = self.plan[step_index]['ang']
    zero_vel = np.zeros(3)
    zero_acc = np.zeros(3)

    # assemble pose, velocity, and acceleration for each foot
    support_data = {
        'pos': np.hstack((support_ang, support_pos)),
        'vel': np.hstack((np.zeros(3), zero_vel)),
        'acc': np.hstack((np.zeros(3), zero_acc))
    }

    swing_data = {
        'pos': np.hstack((swing_ang_pos, swing_pos)),
        'vel': np.hstack((swing_ang_vel, swing_vel)),
        'acc': np.hstack((swing_ang_acc, swing_acc))
    }
    return {
        support_foot: support_data,
        swing_foot: swing_data
    }