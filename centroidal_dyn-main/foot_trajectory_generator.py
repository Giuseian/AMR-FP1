import numpy as np

class FootTrajectoryGenerator:
  def __init__(self, initial, footstep_planner, params):
    self.delta = params['world_time_step']
    self.step_height = params['step_height']
    self.initial = initial
    self.footstep_planner = footstep_planner
    
    ref_foot = "./outputs/ref_foot_walking.txt"  
    data = np.loadtxt(ref_foot, delimiter=",")
    
    lfoot0 = data[-3:, 0]   
    rfoot1 = data[:3, 1]    
    lfoot1 = data[-3:, 3]      
    rfoot2 = data[:3, 4]      
    lfoot2 = data[-3:, 6]   
    rfoot3 = data[:3, 8]
    lfoot3 = data[-3:, 10]
    rfoot4 = data[:3, 12]
    lfoot4 = data[-3:, 14]
    rfoot5 = data[:3, 16]
    lfoot5 = data[-3:, 18]
    rfoot6 = data[:3, 20]
    lfoot6 = data[-3:, 22]
    rfoot7 = data[:3, 24]
    foot = [lfoot0, rfoot1, lfoot1, rfoot2, lfoot2, rfoot3, lfoot3, rfoot4, lfoot4, rfoot5, lfoot5, rfoot6, lfoot6, rfoot7]
    
    phases_durations = "./outputs/phase_durations_walking.txt"
    duration = np.loadtxt(phases_durations, delimiter=",")
    p1 = duration[0]*100
    p2 = duration[1]*100
    p3 = duration[2]*100
    p4 = duration[3]*100
    p5 = duration[4]*100
    p6 = duration[5]*100
    p7 = duration[6]*100
    p8 = duration[7]*100
    p9 = duration[8]*100
    p10 = duration[9]*100
    p11 = duration[10]*100
    p12 = duration[11]*100
    p13 = duration[12]*100
    p14 = duration[13]*100
    p15 = duration[14]*100
    p16 = duration[15]*100
    p17 = duration[16]*100
    p18 = duration[17]*100
    p19 = duration[18]*100
    p20 = duration[19]*100
    p21 = duration[20]*100
    p22 = duration[21]*100
    p23 = duration[22]*100
    p24 = duration[23]*100
    ang = np.array([0.0, 0.0, 0.0])
    
    self.plan = [
        { 'pos': lfoot0, 'ang': ang, 'ss_duration': 0,   'ds_duration': p1, 'foot_id': 'lfoot' },   # 0,0

        { 'pos': rfoot1, 'ang': ang, 'ss_duration': p2,  'ds_duration': p3, 'foot_id': 'rfoot' },   # 0,-1

        { 'pos': lfoot1, 'ang': ang, 'ss_duration': p3,  'ds_duration': p4, 'foot_id': 'lfoot' },    # 0, 0

        { 'pos': rfoot2, 'ang': ang, 'ss_duration': p5,  'ds_duration': p6, 'foot_id': 'rfoot' },   # -1, 0

        { 'pos': lfoot2, 'ang': ang, 'ss_duration': p6,  'ds_duration': p7, 'foot_id': 'lfoot' },   # 0, 0

        { 'pos': rfoot3, 'ang': ang, 'ss_duration': p8,  'ds_duration': p9, 'foot_id': 'rfoot' },   # 0, -1
        
        { 'pos': lfoot3, 'ang': ang, 'ss_duration': p10,  'ds_duration': p11, 'foot_id': 'lfoot' },   

        { 'pos': rfoot4, 'ang': ang, 'ss_duration': p12,  'ds_duration': p13, 'foot_id': 'rfoot' },   

        { 'pos': lfoot4, 'ang': ang, 'ss_duration': p13,  'ds_duration': p14, 'foot_id': 'lfoot' },   

        { 'pos': rfoot5, 'ang': ang, 'ss_duration': p15,  'ds_duration': p16, 'foot_id': 'rfoot' },   

        { 'pos': lfoot5, 'ang': ang, 'ss_duration': p17,  'ds_duration': p18, 'foot_id': 'lfoot' },   

        { 'pos': rfoot6, 'ang': ang, 'ss_duration': p19,  'ds_duration': p20, 'foot_id': 'rfoot' },   

        { 'pos': lfoot6, 'ang': ang, 'ss_duration': p21,  'ds_duration': p22, 'foot_id': 'lfoot' },   

        { 'pos': rfoot7, 'ang': ang, 'ss_duration': p23,  'ds_duration': p24, 'foot_id': 'rfoot' },   


    ]
    self.footstep_planner.plan = self.plan
    
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