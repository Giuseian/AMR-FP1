import numpy as np

def ref_trajectory_generation(n_e, N, ref_type, sigma):
    
    X_ref = None
    U_ref = None
    
    start_pos = np.array([-6.17867734e-04, 4.43297775e-04, 7.23981584e-01])      # Initial CoM position
    start_orient = np.array([1,0,0,0])
    start_vel = np.zeros(3)
    start_feet_pos = np.array([[1.03109240e-17, -1.01638576e-01, -1.38777878e-17], [1.03109240e-17, 1.01638576e-01, -1.38777878e-17]])  # Initial feet positions (distance 0.5 from CoM)
    start_feet_orient = np.array([[1, 0, 0, 0], [1, 0, 0, 0]])  # Neutral orientations
        
        
    if ref_type == "still": # immobile
        X_ref = np.zeros((28, N+1))
        U_ref = np.zeros((27, N)) + 0.000001
        time_k = 0
        phase_duration = 1 # tau
        for t in range(N+1):
            X_ref[0:3, t] = start_pos                         # Initial CoM position
            X_ref[3:7, t] = start_orient
            X_ref[7:10, t] = start_vel
            X_ref[13:14, t] = time_k
            time_k += phase_duration
            X_ref[14:14+3*n_e, t] = start_feet_pos.flatten()  # Initial feet positions
            X_ref[14 + 3*n_e:14 + 3*n_e + 4*n_e, t] = start_feet_orient.flatten()
            
        for t in range(N):
            U_ref[1 + 6*n_e : 1 + 7*n_e, t] = [np.sqrt(9.81),np.sqrt(9.81)]
            U_ref[0, t] = phase_duration
        
        return X_ref, U_ref

    if ref_type == "walking":
        X_ref = np.zeros((28, N+1))
        U_ref = np.zeros((27, N)) + 0.000001
        time_k = 0
        start_vel = [0, -0.07, 0]
        X_ref[0:3, 0] = start_pos                         # Initial CoM position
        X_ref[3:7, 0] = start_orient
        X_ref[7:10, 0] = start_vel
        X_ref[13:14, 0] = time_k
        X_ref[14:14+3*n_e, 0] = start_feet_pos.flatten()  # Initial feet positions
        X_ref[14 + 3*n_e:14 + 3*n_e + 4*n_e, 0] = start_feet_orient.flatten()

        foot_vel = np.zeros(6)

        sig_idx = 0
        #print("sigma" , sigma)
        for t in range(N):
            right_contact = sigma[sig_idx]
            left_contact = sigma[sig_idx+1]
            
            lambda_right, lambda_left = np.sqrt(9.81),np.sqrt(9.81)
            phase_duration = 1  # ds
            if right_contact == -1 and left_contact == 0:
                phase_duration = 1
                lambda_right = 0
                lambda_left = np.sqrt(9.81/0.5)
                
            if left_contact == -1 and right_contact == 0:
                phase_duration = 1
                lambda_left = 0
                lambda_right = np.sqrt(9.81/0.5)

            U_ref[1 + 6*n_e : 1 + 7*n_e, t] = [lambda_right, lambda_left]
            U_ref[0, t] = phase_duration
            sig_idx += 2

        sig_idx = 2
        for t in range(1, N+1):
            right_contact = sigma[sig_idx]  
            left_contact = sigma[sig_idx+1]
            c_v_x = 0.1
            if t == 2:
                c_v_x /= 2
            
            # update foot position and velocity
            right_vel_x, left_vel_x, right_vel_z, left_vel_z = 0.0, 0.0, 0.0, 0.0
            phase_duration = 1    # ss
            
            if right_contact == -1 and sigma[sig_idx-2] == 0:   # rf 0 -> -1
                phase_duration = 1  # ds
                c_v_x = 0 
                c_v_y = 0.07
            
            if left_contact == -1 and sigma[sig_idx-1] == 0:    # lf 0 -> -1
                phase_duration = 1  # ds
                c_v_x = 0
                c_v_y = -0.07
                
            if right_contact == 0 and sigma[sig_idx-2] == -1:   # rf -1 -> 0
                right_vel_x = c_v_x*2
                c_v_y = -0.07

            if left_contact == 0 and sigma[sig_idx-1] == -1:    # lf -1 -> 0
                left_vel_x = c_v_x*2
                c_v_y = 0.07

            # update com position 
            com_vel = np.array([c_v_x, c_v_y, 0.0])
            X_ref[7:10, t-1] = com_vel   # constant com velocity
            com_ds = com_vel * U_ref[0, t-1]
            X_ref[0:3, t] = X_ref[0:3, t-1] + com_ds
            
            # update time
            time_k += phase_duration
            X_ref[13:14, t] = time_k
            
            feet_vel = np.array([right_vel_x, 0.0, right_vel_z, left_vel_x, 0.0, left_vel_z])
            U_ref[1:1+3*n_e, t-1] = feet_vel
            right_ds = feet_vel[:3]*U_ref[0, t-1]
            left_ds = feet_vel[3:6]*U_ref[0, t-1]
            X_ref[14:17, t] = X_ref[14:17, t-1] + right_ds
            X_ref[17:20, t] = X_ref[17:20, t-1] + left_ds
            X_ref[3:7, t] = start_orient
            X_ref[14 + 3*n_e:14 + 3*n_e + 4*n_e, t] = start_feet_orient.flatten()
            sig_idx += 2
            
        return X_ref, U_ref

    if ref_type == "running":
        X_ref = np.zeros((28, N+1))
        U_ref = np.zeros((27, N)) + 0.000001
        time_k = 0
        start_vel = [0, -0.07, 0]
        X_ref[0:3, 0] = start_pos                         # Initial CoM position
        X_ref[3:7, 0] = start_orient
        X_ref[7:10, 0] = start_vel
        X_ref[13:14, 0] = time_k
        X_ref[14:14+3*n_e, 0] = start_feet_pos.flatten()  # Initial feet positions
        X_ref[14 + 3*n_e:14 + 3*n_e + 4*n_e, 0] = start_feet_orient.flatten()

        sig_idx = 0
        for t in range(N):
            right_contact = sigma[sig_idx]
            left_contact = sigma[sig_idx+1]
            
            lambda_right, lambda_left = np.sqrt(9.81),np.sqrt(9.81)
            phase_duration = 0.2  # ds
            if right_contact == -1 and left_contact == 0:
                phase_duration = 0.4
                lambda_right = 0
                lambda_left = np.sqrt(9.81/0.5)
                
            if left_contact == -1 and right_contact == 0:
                phase_duration = 0.4
                lambda_left = 0
                lambda_right = np.sqrt(9.81/0.5)

            U_ref[1 + 6*n_e : 1 + 7*n_e, t] = [lambda_right, lambda_left]
            U_ref[0, t] = phase_duration
            sig_idx += 2
        
        sig_idx = 2
        for t in range(1, N+1):
            right_contact = sigma[sig_idx]  
            left_contact = sigma[sig_idx+1]
            c_v_x = 0.2
            if t == 2:
                c_v_x /= 2
            
            # update foot position and velocity
            right_vel_x, left_vel_x, right_vel_z, left_vel_z = 0.0, 0.0, 0.0, 0.0
            phase_duration = 0.4    # ss
            
            if right_contact == -1 and sigma[sig_idx-2] == 0:   # rf 0 -> -1
                phase_duration = 0.2  # ds
                c_v_x = 0 
                c_v_y = 0.08
            
            if left_contact == -1 and sigma[sig_idx-1] == 0:    # lf 0 -> -1
                phase_duration = 0.2  # ds
                c_v_x = 0
                c_v_y = -0.08
                
            if right_contact == 0 and sigma[sig_idx-2] == -1:   # rf -1 -> 0
                right_vel_x = c_v_x*2
                c_v_y = -0.08

            if left_contact == 0 and sigma[sig_idx-1] == -1:    # lf -1 -> 0
                left_vel_x = c_v_x*2
                c_v_y = 0.08

            # update com position 
            com_vel = np.array([c_v_x, c_v_y, 0.0])
            X_ref[7:10, t-1] = com_vel   # constant com velocity
            com_ds = com_vel * U_ref[0, t-1]
            X_ref[0:3, t] = X_ref[0:3, t-1] + com_ds
            
            # update time
            time_k += phase_duration
            X_ref[13:14, t] = time_k
            
            feet_vel = np.array([right_vel_x, 0.0, right_vel_z, left_vel_x, 0.0, left_vel_z])
            right_ds = feet_vel[:3]*U_ref[0, t-1]
            left_ds = feet_vel[3:6]*U_ref[0, t-1]
            X_ref[14:17, t] = X_ref[14:17, t-1] + right_ds
            X_ref[17:20, t] = X_ref[17:20, t-1] + left_ds
            X_ref[3:7, t] = start_orient
            X_ref[14 + 3*n_e:14 + 3*n_e + 4*n_e, t] = start_feet_orient.flatten()
            sig_idx += 2
            
        
        return X_ref, U_ref