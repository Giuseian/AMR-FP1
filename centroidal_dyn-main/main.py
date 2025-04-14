import casadi as cs # type: ignore
from graphics import *
from trajectory_generation import *
from sbcdyn import * 

if __name__ == "__main__":
    sigma = {"still":   cs.DM([[0],
                            [0]]),
            "walking" : cs.DM([[0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0],
                           [0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0] ])
            #"walking" : cs.DM([[0, 0, 0, -1, 0, 0, 0],
            #               [0, -1, 0, 0, 0, -1, 0] ])
            
            #"walking" : cs.DM([[0, 0, 0],
            #                   [0, -1, 0]])
            }

    ref = "walking"
    n_e, n_div, surfaces, N, opti = 2, 10, [0], sigma[ref].shape[1]-1, cs.Opti()
    X = opti.variable(14+7*n_e, N+1)
    U = opti.variable(1+13*n_e, N)
    X_ref = opti.parameter(14+7*n_e, N+1)
    U_ref = opti.parameter(1+13*n_e, N)

    print(f"Solving for {ref} reference...\n")
    X_init_ref, U_init_ref = ref_trajectory_generation(n_e, N, ref, sigma[ref])

    X_init = np.full(X.shape, 0.7)
    U_init = np.full(U.shape, 0.5)
    opti.set_initial(X, X_init)
    opti.set_initial(U, U_init)
    opti.set_value(X_ref, X_init_ref)
    opti.set_value(U_ref, U_init_ref)

    dynamic_model = StiffnessBasedCentroidalDynamics(n_e, n_div, surfaces, N, sigma[ref],opti)
    X_sol, U_sol = dynamic_model.solve(X, U, X_ref, U_ref)
    full_array, phases_duration = dynamic_model.time_domain_solution(X_sol, U_sol)
    
    ref_foot_traj = X_init_ref[14:20]
    com_pos = full_array[:, 0:3]
    com_vel = full_array[:, 7:10]
    time_durations = X_sol[13:14]
    phases_durations = U_sol[0:1]

    np.savetxt("./outputs/ref_foot.txt", ref_foot_traj, delimiter=",", fmt="%.8f") 
    np.savetxt("./outputs/com_pos.txt", com_pos, delimiter=",", fmt="%.8f")  
    np.savetxt("./outputs/com_vel.txt", com_vel, delimiter=",", fmt="%.8f")  
    np.savetxt("./outputs/time_durations.txt", time_durations, delimiter=",", fmt="%.8f") 
    np.savetxt("./outputs/phase_durations.txt", phases_durations, delimiter=",", fmt="%.8f") 
    print("\n LOGS saved under outputs folders...\n")
    plot_components(full_array, phases_duration)
    evolution_contact_forces(X=X_sol, X_ref=X_init_ref, U=U_sol, U_ref=U_init_ref, dm=dynamic_model, ref_type=ref, opti=opti)
    evolution_plots(X=X_sol, U=U_sol, X_ref=X_init_ref, U_ref=U_init_ref, dm=dynamic_model, ref_type=ref)
    checking_sum_forces(X=X_sol, X_ref=X_init_ref, U = U_sol, U_ref = U_init_ref, dm=dynamic_model, ref_type=ref, opti=opti)
    print("\n Plots saved under plots folders...\n")
    animate_trajectories_td(full_array, n_e=2, save_path="./videos/walking_in_time.gif", ref_type=ref, label="Actual")
    animate_trajectories(X_sol, n_e=2, save_path="./videos/walking_in_contacts.gif", ref_type=ref, label="Actual")
    print("\n Videos saved under videos folder ...\nExit!")
    