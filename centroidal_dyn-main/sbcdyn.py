import casadi as cs # type: ignore
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from matplotlib.animation import FuncAnimation # type: ignore
import matplotlib.animation as animation # type: ignore
from mpl_toolkits.mplot3d import Axes3D # type: ignore
from scipy.interpolate import interp1d # type: ignore

class StiffnessBasedCentroidalDynamics:

    def __init__(self, n_e, n_div, surfaces, N, sigma, opti):
        
        self.n_e = n_e
        self.n_div = n_div
        self.surfaces = surfaces
        self.N = N
        self.sigma = sigma
        self.opti = opti
        
        self.epsilon = 1e-6
        self.g = [0,0,9.81]
        self.mu = 5
        self.mu_z = 0.6
        self.tau_min = 0.4
        self.tau_max = 5
        self.m = 10
        self.LAMBDA_max = 1000
        self.feet_length = 0.1
        
        # for cost
        self.SIGMA_L_k = opti.parameter(n_e, N+1)  # (n_e, N)  
        self.ETA_NORMAL_L = opti.parameter(3, len(surfaces))  # (3, surfaces)  # (3,3)
        self.O_L = opti.parameter(3, len(surfaces))   # (3, surfaces)          # (3,3)
        self.w_compl = opti.parameter(1)  # Complementarity weight
        self.W_x_k = opti.parameter(28,28)  # State weight matrix
        self.W_u_k = opti.parameter(27,27)  # Input weight matrix
        self.x_0 = opti.parameter(28)

        weights_x = [1.0]*7 + [0.0001]*6 + [1]*15  # Higher weight for the first 18 components
        W_x_k = cs.diag(weights_x)  # Create a diagonal matrix from the weights
        weights_u = [1] + [0.0001] * 12 + [1]*2 + [1]*12# Higher weight for the first 18 components
        W_u_k = cs.diag(weights_u)  # Create a diagonal matrix from the weights
        
        # Set the weight matrix in the optimizer
        opti.set_value(self.W_x_k, W_x_k)
        opti.set_value(self.W_u_k, W_u_k)
        opti.set_value(self.ETA_NORMAL_L, cs.DM([0, 0, 1]))
        opti.set_value(self.O_L, cs.DM([0, 0, 0]))
        opti.set_value(self.SIGMA_L_k, sigma)
        opti.set_value(self.w_compl, 1000)
        
    def init_state(self, X):
        n_e = self.n_e
        
        X_init = np.zeros(X.shape[0])
        start_pos = np.array([-6.17867734e-04, 4.43297775e-04, 7.23981584e-01])      # Initial CoM position
        start_orient = np.array([1,0,0,0])
        start_vel = np.zeros(3)
        start_feet_pos = np.array([[1.03109240e-17, -1.01638576e-01, -1.38777878e-17], [1.03109240e-17, 1.01638576e-01, -1.38777878e-17]])  # Initial feet positions (distance 0.5 from CoM)
        start_feet_orient = np.array([[1, 0, 0, 0], [1, 0, 0, 0]])  # Neutral orientations
        
        X_init[0:3] = start_pos
        X_init[3:7] = start_orient
        X_init[7:10] = start_vel
        X_init[14:14+3*n_e] = start_feet_pos.flatten()
        X_init[14 + 3*n_e:14 + 3*n_e + 4*n_e] = start_feet_orient.flatten()
        return X_init
    
    # defining f_l (+func) and eta_l (+func)
    def define_contact_wrench(self, p_k, LAMBDA_L_k, P_L_k, R_L_k, ETA_HAT_L_k):
        # Taking f_l and eta_l from eq.3.1 and eq.3.2
        # Eq. 3.1 : Translational component of the contact wrench
        P_k = cs.repmat(p_k, 1, self.n_e) # (3, n_e)
        F_L = self.m * cs.mtimes((P_k - (P_L_k + R_L_k)), cs.diag(LAMBDA_L_k**2))  # (3, n_e)
        
        # Eq. 3.2 : Angular component of the contact wrench
        ETA_L = self.m * cs.mtimes(ETA_HAT_L_k, cs.diag(LAMBDA_L_k**2))   # (3, n_e)
        return F_L, ETA_L

    def quaternion(self, w, tau_k_prime):
        w_norm = cs.norm_2(w)
        theta = 0.5 * tau_k_prime * w_norm  # Half-angle for rotation
        q_v = cs.if_else(w_norm > 1e-6, cs.sin(theta) * w / w_norm, cs.DM.zeros(3))
        q_w = cs.cos(theta)
        return cs.vertcat(q_w, q_v)  # Quaternion [w, x, y, z]
    
    def quaternion_product(self,q1, q2):
        w1, x1, y1, z1 = q1[0], q1[1], q1[2], q1[3]
        w2, x2, y2, z2 = q2[0], q2[1], q2[2], q2[3]
    
        w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
        x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
        y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2
        z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2
    
        return cs.vertcat(w, x, y, z)
        
    def update_q_l_k(self, W_L_k, tau_k, q_k):
        q_quaternions = []
        _, col = W_L_k.size()    # it should be n_e
    
        for i in range(col):
            q_temp = self.quaternion(W_L_k[:,i], tau_k)  # (4,1)
            q_prod = self.quaternion_product(q_temp, q_k)   # (4,1)
            q_quaternions.append(q_prod)
    
        q_k_plus_one = cs.horzcat(*q_quaternions)
        
        return q_k_plus_one
    
    def get_x_comp(self, x_k):
        p_k = x_k[0:3]  # (3,1), corresponds to p_k (position)
        q_k = x_k[3:7]  # (4,1), corresponds to q_k (orientation as quaternion)
        v_k = x_k[7:10]  # (3,1), corresponds to v_k (velocity)
        L_k = x_k[10:13]  # (3,1), corresponds to L_k (angular momentum)
        t_k = x_k[13:14]  # (1,1), current instant of time (instant of time inside a phase)
        P_L_k_flat = x_k[14:14+3*self.n_e]  # (3*n_e, 1) positions of end-effectors
        Q_L_k_flat = x_k[14+3*self.n_e:14+3*self.n_e+4*self.n_e]  # (4*n_e, 1) orientations (quaternions) of end-effectors
        P_L_k = cs.reshape(P_L_k_flat, 3, self.n_e)
        Q_L_k = cs.reshape(Q_L_k_flat, 4, self.n_e)
        return p_k, q_k, v_k, L_k, t_k, P_L_k, Q_L_k

    def get_u_comp(self, u_k):
        tau_k = u_k[0]
        V_L_k_flat = u_k[1:1+3*self.n_e]
        W_L_k_flat = u_k[1+3*self.n_e:1+6*self.n_e]
        LAMBDA_L_k = u_k[1+6*self.n_e:1+7*self.n_e]
        R_L_k_flat = u_k[1+7*self.n_e:1+10*self.n_e]
        ETA_HAT_L_k_flat = u_k[1+10*self.n_e:1+13*self.n_e]
        V_L_k = cs.reshape(V_L_k_flat, 3, self.n_e)
        W_L_k = cs.reshape(W_L_k_flat, 3, self.n_e)
        R_L_k = cs.reshape(R_L_k_flat, 3, self.n_e)
        ETA_HAT_L_k = cs.reshape(ETA_HAT_L_k_flat, 3, self.n_e)
        return tau_k, V_L_k, W_L_k, LAMBDA_L_k, R_L_k, ETA_HAT_L_k
    
    # returns x_k+1
    def update_dynamics(self, x_k, u_k):
        # extracting components from x_k
        p_k, q_k, v_k, L_k, t_k, P_L_k, Q_L_k = self.get_x_comp(x_k)
        
        # extracting components from u_k
        tau_k, V_L_k, W_L_k, LAMBDA_L_k, R_L_k, ETA_HAT_L_k = self.get_u_comp(u_k)
        
        # Eq. 5.1 : lambda_bar
        lambda_bar = cs.sqrt(cs.sumsqr(LAMBDA_L_k) + self.epsilon**2)  # (1,1)
        
        # Eq. 5.2 : p_bar
        p_bar = ((P_L_k @ (LAMBDA_L_k**2)) + self.g) / lambda_bar**2  # (3,1)
        
        # Eq. 5.3 : r_bar
        r_bar = (R_L_k @ (LAMBDA_L_k**2)) / lambda_bar**2  # (3,1)
        
        # Eq. 5.4 : eta_bar
        eta_bar = (
            lambda_bar**2 * cs.cross(p_bar, r_bar) +
            (ETA_HAT_L_k - cs.cross(P_L_k, R_L_k, 1)) @ (LAMBDA_L_k**2)
        )  # (3,1)

        C_k = cs.cosh(lambda_bar * tau_k)    # C_k(t-t_k), scalar
        S_k = cs.sinh(lambda_bar * tau_k)    # S_k(t-t_k), scalar

        p_k_plus_one = (
            p_bar + r_bar +
            C_k * (p_k - (p_bar + r_bar)) +
            (S_k / lambda_bar) * v_k
        )
        
        v_k_plus_one = (
            lambda_bar * S_k * (p_k - (p_bar + r_bar)) +
            C_k * v_k
        )

        L_k_plus_one = (
            L_k +
            self.m * (
                cs.cross((v_k_plus_one - v_k), r_bar) +
                tau_k * eta_bar
            )
        )

        t_k_plus_one = t_k + tau_k

        # update q_k
        w, x, y, z = q_k[0], q_k[1], q_k[2], q_k[3]
        
        R_t = cs.vertcat(
            cs.horzcat(1 - 2 * (y**2 + z**2), 2 * (x * y - w * z), 2 * (x * z + w * y)),
            cs.horzcat(2 * (x * y + w * z), 1 - 2 * (x**2 + z**2), 2 * (y * z - w * x)),
            cs.horzcat(2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x**2 + y**2))
        )
        
        tau_k_prime = tau_k / self.n_div
        q_new = q_k
    
        for i in range(self.n_div):
            t_i_prime = t_k + tau_k_prime * i  
            C_prime = cs.cosh(lambda_bar * (t_i_prime - t_k))    # C_k(t-t_k), scalar
            S_prime = cs.sinh(lambda_bar * (t_i_prime - t_k))
            v_t_prime = ( lambda_bar * S_prime * (p_k - (p_bar + r_bar)) + C_prime * v_k )
            L_t_prime = (L_k + self.m * (cs.cross((v_t_prime - v_k), r_bar) + (t_i_prime - t_k) * eta_bar) )
        
            # Reference linear displacement
            L_ref_t = cs.DM.eye(3)  # Assuming L_ref = 0
        
            # Inverse of reference inertia matrix
            I_ref_inv = cs.inv(cs.diag(cs.DM([0.25, 0.25, 0.15])))  # Constant matrix
        
            # Compute w_t
            w_t = R_t @ I_ref_inv @ (R_t.T @ L_t_prime - L_ref_t)

            t_i_prime = t_k + tau_k_prime * i
            omega =  w_t * tau_k_prime
            q_temp = self.quaternion(w, tau_k_prime)
            q_new = self.quaternion_product(q_temp, q_new)

        q_k_plus_one = q_new
        
        p_l_k_plus_one = P_L_k + V_L_k * tau_k
        p_l_k_plus_one_flat = cs.reshape(p_l_k_plus_one, (p_l_k_plus_one.shape[0]*p_l_k_plus_one.shape[1], 1))
        
        q_l_k_plus_one = self.update_q_l_k(W_L_k, tau_k, q_k)
        q_l_k_plus_one_flat = cs.reshape(q_l_k_plus_one, (q_l_k_plus_one.shape[0]*q_l_k_plus_one.shape[1], 1))
        
        x_k_plus_one = cs.vertcat(
                p_k_plus_one,              # Position update
                q_k_plus_one,              # Quaternion update
                v_k_plus_one,              # Velocity update
                L_k_plus_one,              # Angular momentum update
                t_k_plus_one,              # Time update
                p_l_k_plus_one_flat,       # End position update
                q_l_k_plus_one_flat        # Flattened orientation matrix update
            )
        
        return x_k_plus_one
    
    
    # GET INTRA DYNAMICS 
    def get_x_t(self, x_k, u_k, t):
        # extracting components from x_k
        p_k, q_k, v_k, L_k, t_k, P_L_k, Q_L_k = self.get_x_comp(x_k)
        
        # extracting components from u_k
        tau_k, V_L_k, W_L_k, LAMBDA_L_k, R_L_k, ETA_HAT_L_k = self.get_u_comp(u_k)
        
        # Eq. 5.1 : lambda_bar
        lambda_bar = cs.sqrt(cs.sumsqr(LAMBDA_L_k) + self.epsilon**2)  # (1,1)
        
        # Eq. 5.2 : p_bar
        p_bar = ((P_L_k @ (LAMBDA_L_k**2)) + self.g) / lambda_bar**2  # (3,1)
        
        # Eq. 5.3 : r_bar
        r_bar = (R_L_k @ (LAMBDA_L_k**2)) / lambda_bar**2  # (3,1)
        
        # Eq. 5.4 : eta_bar
        eta_bar = (
            lambda_bar**2 * cs.cross(p_bar, r_bar) +
            (ETA_HAT_L_k - cs.cross(P_L_k, R_L_k, 1)) @ (LAMBDA_L_k**2)
        )  # (3,1)

        C_k = cs.cosh(lambda_bar * t )    # C_k(t-t_k), scalar
        S_k = cs.sinh(lambda_bar * t )    # S_k(t-t_k), scalar

        p_k_plus_one = (
            p_bar + r_bar +
            C_k * (p_k - (p_bar + r_bar)) +
            (S_k / lambda_bar) * v_k
        )
        
        v_k_plus_one = (
            lambda_bar * S_k * (p_k - (p_bar + r_bar)) +
            C_k * v_k
        )

        L_k_plus_one = (
            L_k +
            self.m * (
                cs.cross((v_k_plus_one - v_k), r_bar) +
                tau_k * eta_bar
            )
        )

        t_k_plus_one = t_k + t 

        # update q_k
        w, x, y, z = q_k[0], q_k[1], q_k[2], q_k[3]
        
        R_t = cs.vertcat(
            cs.horzcat(1 - 2 * (y**2 + z**2), 2 * (x * y - w * z), 2 * (x * z + w * y)),
            cs.horzcat(2 * (x * y + w * z), 1 - 2 * (x**2 + z**2), 2 * (y * z - w * x)),
            cs.horzcat(2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x**2 + y**2))
        )
        
        tau_k_prime = t / self.n_div
        q_new = q_k
    
        for i in range(self.n_div):
            t_i_prime = t_k + tau_k_prime * i  
            C_prime = cs.cosh(lambda_bar * (t_i_prime - t_k))    # C_k(t-t_k), scalar
            S_prime = cs.sinh(lambda_bar * (t_i_prime - t_k))
            v_t_prime = ( lambda_bar * S_prime * (p_k - (p_bar + r_bar)) + C_prime * v_k )
            L_t_prime = (L_k + self.m * (cs.cross((v_t_prime - v_k), r_bar) + (t_i_prime - t_k) * eta_bar) )
        
            # Reference linear displacement
            L_ref_t = cs.DM.eye(3)  # Assuming L_ref = 0
        
            # Inverse of reference inertia matrix
            I_ref_inv = cs.inv(cs.diag(cs.DM([2, 3, 4])))  # Constant matrix
        
            # Compute w_t
            w_t = R_t @ I_ref_inv @ (R_t.T @ L_t_prime - L_ref_t)

            t_i_prime = t_k + tau_k_prime * i
            omega =  w_t * tau_k_prime
            q_temp = self.quaternion(w, tau_k_prime)
            q_new = self.quaternion_product(q_temp, q_new)

        q_k_plus_one = q_new
        
        p_l_k_plus_one = P_L_k + V_L_k * t
        p_l_k_plus_one_flat = cs.reshape(p_l_k_plus_one, (p_l_k_plus_one.shape[0]*p_l_k_plus_one.shape[1], 1))
        
        q_l_k_plus_one = self.update_q_l_k(W_L_k, t, q_k)
        q_l_k_plus_one_flat = cs.reshape(q_l_k_plus_one, (q_l_k_plus_one.shape[0]*q_l_k_plus_one.shape[1], 1))
        
        x_k_plus_one = cs.vertcat(
                p_k_plus_one,              # Position update
                q_k_plus_one,              # Quaternion update
                v_k_plus_one,              # Velocity update
                L_k_plus_one,              # Angular momentum update
                t_k_plus_one,              # Time update
                p_l_k_plus_one_flat,       # End position update
                q_l_k_plus_one_flat        # Flattened orientation matrix update
            )
        
        return x_k_plus_one


    def constraint_lambda(self, p_k_ref, P_L_k_ref, LAMBDA_L_k_ref, g):
        diff = cs.mtimes(cs.repmat(p_k_ref, 1, self.n_e) - P_L_k_ref, LAMBDA_L_k_ref**2) - g
        lambda_obj = 0.5 * cs.mtimes(diff.T, diff)
        self.opti.minimize(lambda_obj)

        self.opti.solver("ipopt")
        solution = self.opti.solve()        
        LAMBDA_L_sol = solution.value(LAMBDA_L_k_ref)
        
        return LAMBDA_L_sol

    def constraint_box(self, p_k, P_L_k):

        self.opti.subject_to(cs.sumsqr(P_L_k[0:3] - p_k) <= 1) # euclidian distance have square root, but it cause optimization problems
        self.opti.subject_to(cs.sumsqr(P_L_k[3:6] - p_k) <= 1) # but we have an equivalent formulation by setting <= maxdistance^2 (1**2 = 1vin this case)
        
        return 
    
    def constraint_slip_condition(self, F_L, k):
        SIGMA_L = np.array(self.opti.value(self.SIGMA_L_k)).astype("int")
    
        for i in range(self.n_e):
            
            if SIGMA_L[i,k] == 0:                
                f_l_x, f_l_y, f_l_z = F_L[0, i], F_L[1, i], F_L[2, i]
                
                # Eq. 20: Non-slip condition
                self.opti.subject_to(cs.fabs(self.mu * f_l_z) >= cs.sqrt(f_l_x**2 + f_l_y**2))
        
        return 
    
    def constraint_vector(self, q_k, p_k, LAMBDA_L_k, P_L_k, R_L_k, ETA_HAT_L_k, tau_k):
        F_L, ETA_L = self.define_contact_wrench(p_k, LAMBDA_L_k, P_L_k, R_L_k, ETA_HAT_L_k)
        mu = self.mu
        mu_z = self.mu_z 
        feet_length = self.feet_length
        g_vector = []
        
        g_vector.append(tau_k - self.tau_min)
        g_vector.append(self.tau_max - tau_k)
        g_vector.append(LAMBDA_L_k)
        g_vector.append(self.LAMBDA_max - LAMBDA_L_k)
        
        for i in range(self.n_e):
                        
            f_l_x, f_l_y, f_l_z = F_L[0, i], F_L[1, i], F_L[2, i]
            eta_l_x, eta_l_y, eta_l_z = ETA_L[0, i], ETA_L[1, i], ETA_L[2, i]
            
            # Eq. 20: Non-slip condition
            g_vector.append(mu * f_l_z - cs.sqrt(f_l_x**2 + f_l_y**2))
            
            # Eq. 21abc: Contact moment for eta_l_x y and z
            """
            g_vector.append(eta_l_x /  f_l_z + feet_length)
            g_vector.append(- eta_l_x /  f_l_z - feet_length)
            g_vector.append(eta_l_y / f_l_z - feet_length)
            g_vector.append(- eta_l_y / f_l_z + feet_length)
            g_vector.append(eta_l_z / f_l_z + mu_z)
            g_vector.append(-eta_l_z   /  f_l_z + mu_z)
            """
        
        # Convert g_vector to a CasADi MX vector
        g_constraint = cs.vertcat(*g_vector)  # (14, 1) = (n_g, 1), where n_g is the total number of inequality constraints
        return g_constraint

    def limit_cost(self, q_k, p_k, LAMBDA_L_k, P_L_k, R_L_k, ETA_HAT_L_k, tau_k):
        # Log-Barrier Function
        m = self.m
        
        epsilon = self.epsilon
        g_constraint = self.constraint_vector(q_k, p_k, LAMBDA_L_k, P_L_k, R_L_k, ETA_HAT_L_k, tau_k)
        
        L_limit = -cs.sum1(cs.log(cs.fmax(epsilon, g_constraint)))  # Log-barrier term
        
        return L_limit

    def contact_dependent_cost(self, P_L_k, V_L_k, W_L_k, LAMBDA_L_k, SIGMA_L_k):
        # Initialize the cost
        L_compl_k = 0
        
        for l in range(self.n_e):  # Loop over ends
            # First term: Distance to contact face
            contact_cost = 0
            for i in self.surfaces:    # SIGMA_L_k (n_e, N)
                delta_contact = cs.if_else(SIGMA_L_k[l] == i, 1, 0)
                contact_cost += delta_contact * (cs.dot(self.ETA_NORMAL_L[:, i], P_L_k[:, l] - self.O_L[:, i]) ** 2)
        
            # Second term: End velocity cost
            delta_velocity = cs.if_else(SIGMA_L_k[l] != -1, 1, 0)  # Not off-contact
            velocity_cost = delta_velocity * (cs.norm_2(V_L_k[:, l]) ** 2 + cs.norm_2(W_L_k[:, l]) ** 2)
        
            # Third term: Stiffness cost
            delta_off_contact = cs.if_else(SIGMA_L_k[l] == -1, 1, 0)  # Off-contact
            stiffness_cost = delta_off_contact * (LAMBDA_L_k[l] ** 2)
        
            # Combine terms
            L_compl_k += contact_cost + stiffness_cost + velocity_cost
        
        # Scale by w_compl^2
        L_compl_k *= self.w_compl ** 2
        return L_compl_k
    
    # defining task-related cost function
    def task_related_cost(self, x_k, u_k, x_k_ref, u_k_ref):
        W_x_k = self.W_x_k  # State weight matrix
        W_u_k = self.W_u_k  # Input weight matrix
        
        # Eq. 16 : Task-related cost function
        L_task_k = 0.5 * cs.mtimes([(x_k - x_k_ref).T, W_x_k, (x_k- x_k_ref)]) + \
                   0.5 * cs.mtimes([(u_k - u_k_ref).T, W_u_k, (u_k - u_k_ref)])
        
        return L_task_k

    def impose_constraints(self, X, U, X_ref, U_ref):
        
        self.opti.subject_to()
        N_intervals = U.shape[1]
        cost = 0
        SIGMA_L = self.SIGMA_L_k
        
        for k in range(N_intervals): # number of phases 
            x_k = X[:, k]
            u_k = U[:, k]
            x_k_ref = X_ref[:, k]
            u_k_ref = U_ref[:, k]

            p_k, q_k, v_k, L_k, t_k, P_L_k, Q_L_k = self.get_x_comp(x_k)
            
            tau_k, V_L_k, W_L_k, LAMBDA_L_k, R_L_k, ETA_HAT_L_k = self.get_u_comp(u_k)
            
            p_k_ref, q_k_ref, v_k_ref, L_k_ref, t_k_ref, P_L_k_ref, Q_L_k_ref = self.get_x_comp(x_k_ref)
        
            tau_k_ref, V_L_k_ref, W_L_k_ref, LAMBDA_L_k_ref, R_L_k_ref, ETA_HAT_L_k_ref = self.get_u_comp(u_k_ref) 

            SIGMA_L_k = SIGMA_L[:, k]   # !!!
            task_cost_k = self.task_related_cost(x_k, u_k, x_k_ref, u_k_ref)
            contact_dependent_cost_k = self.contact_dependent_cost(P_L_k, V_L_k, W_L_k, LAMBDA_L_k, SIGMA_L_k)
            #limit_cost_k = self.limit_cost(q_k, p_k, LAMBDA_L_k, P_L_k, R_L_k, ETA_HAT_L_k, tau_k)
            
            cost += task_cost_k + contact_dependent_cost_k #+ limit_cost_k
            x_next = self.update_dynamics(x_k, u_k)
            #lambda_sol = self.constraint_lambda(p_k_ref, P_L_k_ref, LAMBDA_L_k_ref, self.g)
            #self.constraint_box(p_k, P_L_k)
            self.opti.subject_to(X[:, k + 1] == x_next)
            #self.opti.subject_to(tau_k > 0)
            
        last_x = X[:,N_intervals]
        last_x_ref = X_ref[:, N_intervals]
        
        cost += 0.5 * cs.mtimes([(last_x - last_x_ref).T, self.W_x_k, (last_x- last_x_ref)])
        x_0 = self.init_state(X)
        self.opti.subject_to(X[:, 0] == x_0)  # Enforce the initial state
            
        self.opti.minimize(cost)

        return

    def solve(self, X,U, X_ref, U_ref):
        X_sol = None
        U_sol = None
        self.impose_constraints(X, U, X_ref, U_ref)
        self.opti.solver("ipopt", {
                            "ipopt.max_iter": 100000,               # Allow more iterations
                            "ipopt.tol": 1e-4,                      # Overall convergence tolerance (default: 1e-8)
                            "ipopt.constr_viol_tol": 1e-4,          # Constraint violation tolerance (default: 1e-6)
                            "ipopt.compl_inf_tol": 1e-4,            # Complementarity tolerance (default: 1e-6)
                            "ipopt.acceptable_tol": 1e-3,           # Acceptable solution tolerance (less strict)
                            "ipopt.acceptable_constr_viol_tol": 1e-3,  # Acceptable constraint violation
                    })

        solution = self.opti.solve()
        X_sol = solution.value(X)  # Optimal states
        U_sol = solution.value(U)  # Optimal controls
        return X_sol, U_sol
    
    def time_domain_solution(self, X, U):
        N = X.shape[1] - 1
        phases_duration = []
        full_list = []
        for n in range(N):
            intra_list = []
            x_n = X[:,n].reshape(-1, 1)
            t_n_plus_one = X[13:14,n+1].item()
            t_n = x_n[13:14].item()
            duration = int((t_n_plus_one - t_n)// 0.01)
            phases_duration.append(int(duration))
            u_n = U[:,n]
            for _ in range(duration):
                intra_list.append(np.array(x_n))
                x_n = self.get_x_t(x_n, u_n, 0.01)
            full_list += intra_list
        return np.array(full_list).squeeze(2), phases_duration

