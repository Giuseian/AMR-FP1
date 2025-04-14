import casadi as ca
from scipy.spatial.transform import Rotation as R
import numpy as np

def rotation_vector_difference(rotvec_a, rotvec_b):
    R_a = R.from_rotvec(rotvec_a)
    R_b = R.from_rotvec(rotvec_b)
    R_diff = R_b.inv() * R_a
    return R_diff.as_rotvec()

def pose_difference(pose_a, pose_b):
    pos_diff = pose_a[:3] - pose_b[:3]
    rot_diff = rotation_vector_difference(pose_a[3:], pose_b[3:])
    return np.hstack((pos_diff, rot_diff))

# converts a rotation matrix to a rotation vector
def get_rotvec(rot_matrix):
    rotation = R.from_matrix(rot_matrix)
    return rotation.as_rotvec()

def block_diag(*arrays):
    arrays = [np.atleast_2d(a) if np.isscalar(a) else np.atleast_2d(a) for a in arrays]

    rows = sum(arr.shape[0] for arr in arrays)
    cols = sum(arr.shape[1] for arr in arrays)
    block_matrix = np.zeros((rows, cols), dtype=arrays[0].dtype)

    current_row = 0
    current_col = 0

    for arr in arrays:
        r, c = arr.shape
        block_matrix[current_row:current_row + r, current_col:current_col + c] = arr
        current_row += r
        current_col += c

    return block_matrix

# solves a constrained QP with casadi
class QPSolver:
    def __init__(self, n_vars, n_eq_constraints=0, n_ineq_constraints=0):
        self.n_vars = n_vars
        self.n_eq_constraints = n_eq_constraints
        self.n_ineq_constraints = n_ineq_constraints

        self.opti = ca.Opti('conic')
        self.x = self.opti.variable(self.n_vars)

        # objective function: (1/2) * x.T @ H @ x + F.T @ x
        self.F_ = self.opti.parameter(self.n_vars)
        self.H_ = self.opti.parameter(self.n_vars, self.n_vars)
        objective = 0.5 * self.x.T @ self.H_ @ self.x + self.F_.T @ self.x
        self.opti.minimize(objective)

        # equality constraints: A_eq * x == b_eq
        self.A_eq_ = self.opti.parameter(self.n_eq_constraints, self.n_vars)
        self.b_eq_ = self.opti.parameter(self.n_eq_constraints)
        if self.n_eq_constraints > 0:
            self.opti.subject_to(self.A_eq_ @ self.x == self.b_eq_)

        # inequality constraints: A_ineq * x <= b_ineq
        if self.n_ineq_constraints > 0:
            self.A_ineq_ = self.opti.parameter(self.n_ineq_constraints, self.n_vars)
            self.b_ineq_ = self.opti.parameter(self.n_ineq_constraints)
            self.opti.subject_to(self.A_ineq_ @ self.x <= self.b_ineq_)
        else:
            self.A_ineq_ = None
            self.b_ineq_ = None

        # solver options
        p_opts = {'expand': True}
        s_opts = {'max_iter': 1000, 'verbose': False}
        self.opti.solver('osqp', p_opts, s_opts)

    def set_values(self, H, F, A_eq=None, b_eq=None, A_ineq=None, b_ineq=None):
        self.opti.set_value(self.H_, H)
        self.opti.set_value(self.F_, F)
        if self.n_eq_constraints > 0 and A_eq is not None and b_eq is not None:
            self.opti.set_value(self.A_eq_, A_eq)
            self.opti.set_value(self.b_eq_, b_eq)
        if self.n_ineq_constraints > 0 and A_ineq is not None and b_ineq is not None:
            self.opti.set_value(self.A_ineq_, A_ineq)
            self.opti.set_value(self.b_ineq_, b_ineq)

    def solve(self):
        try:
            solution = self.opti.solve()
            x_sol = solution.value(self.x)
        except RuntimeError as e:
            print("QP Solver failed:", e)
            x_sol = np.zeros(self.n_vars)
        return x_sol

def parabolic_trajectory_3D(p_start, p_end, h_max, num_points):
    """
    p_start, p_end: 3D start/end foot positions, e.g. np.array([x, y, z])
    h_max: peak height above the start/end z (assuming start.z = end.z = 0)
    num_points: how many samples we want from start to end
    Returns: a (num_points, 3) array of 3D points
    """
    traj = np.zeros((num_points, 3))
    t_vals = np.linspace(0, 1, num_points)
    
    # We assume the foot starts at p_start[2] and ends at p_end[2].
    # If your terrain is flat, you can set them to 0. 
    # We'll define "flat_z" = average of start/end for reference.
    z_start = p_start[2]
    z_end   = p_end[2]
    
    for i, t in enumerate(t_vals):
        # Linear XY (horizontal) interpolation
        xy = p_start[:2] + t * (p_end[:2] - p_start[:2])
        
        # Parabolic Z that goes from z_start to z_end 
        # but lifts up by h_max in the middle
        # We'll do a standard parabola shape and then shift up if needed:
        z_parab = 4.0 * h_max * t * (1 - t)
        # Add the baseline from start to end
        # so that if z_start == z_end == 0, it’s just the standard parabola.
        z_lin   = z_start + t * (z_end - z_start)
        z_total = z_lin + z_parab
        
        traj[i, 0] = xy[0]
        traj[i, 1] = xy[1]
        traj[i, 2] = z_total
    
    return traj