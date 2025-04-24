import casadi as cs # type: ignore
from graphics import *
from trajectory_generation import *
from sbcdyn import * 
import argparse

if __name__ == "__main__":
    

    parser = argparse.ArgumentParser(description="Trajectory optimization and visualization.")
    parser.add_argument("task", choices=["walking", "running", "still"], help="Task to perform")
    parser.add_argument("--make_video", action="store_true", help="Generate videos if this flag is set.")

    args = parser.parse_args()
    ref = args.task

    sigma = {"still":   cs.DM([[0, 0, 0, 0],
                            [0, 0, 0, 0]]),
            "walking" : cs.DM([[0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0],
                           [0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0] ]),
            "running": cs.DM([
                [0, 0, -1, -1, -1, 0],
                [0, -1, -1, 0, -1, -1]
            ])
            }

    n_e, n_div, surfaces, N, opti = 2, 10, [0], sigma[ref].shape[1]-1, cs.Opti()
    X = opti.variable(14+7*n_e, N+1)
    U = opti.variable(1+13*n_e, N)
    X_ref = opti.parameter(14+7*n_e, N+1)
    U_ref = opti.parameter(1+13*n_e, N)
    print(f"Solving for {ref} reference...\n")
    X_init_ref, U_init_ref = ref_trajectory_generation(n_e, N, ref, sigma[ref])

    X_init = np.full(X.shape, 0.0001)
    U_init = np.full(U.shape, 0.001)
    opti.set_initial(X, X_init)
    opti.set_initial(U, U_init)
    opti.set_value(X_ref, X_init_ref)
    opti.set_value(U_ref, U_init_ref)

    dynamic_model = StiffnessBasedCentroidalDynamics(n_e, n_div, surfaces, N, sigma[ref], opti, ref)
    X_sol, U_sol = dynamic_model.solve(X, U, X_ref, U_ref)
    full_array, phases_duration = dynamic_model.time_domain_solution(X_sol, U_sol)
        
    ref_foot_traj = X_init_ref[14:20]
    com_pos = full_array[:, 0:3]
    com_vel = full_array[:, 7:10]
    time_durations = X_sol[13:14]
    phases_durations = U_sol[0:1]

    np.savetxt(f"./outputs/ref_foot_{ref}.txt", ref_foot_traj, delimiter=",", fmt="%.8f") 
    np.savetxt(f"./outputs/com_pos_{ref}.txt", com_pos, delimiter=",", fmt="%.8f")  
    np.savetxt(f"./outputs/com_vel_{ref}.txt", com_vel, delimiter=",", fmt="%.8f")  
    np.savetxt(f"./outputs/time_durations_{ref}.txt", time_durations, delimiter=",", fmt="%.8f") 
    np.savetxt(f"./outputs/phase_durations_{ref}.txt", phases_durations, delimiter=",", fmt="%.8f") 
    print("\nLOGS saved under outputs folders...\n")
    plot_components(full_array, phases_duration, ref)
    evolution_contact_forces(X=X_sol, X_ref=X_init_ref, U=U_sol, U_ref=U_init_ref, dm=dynamic_model, ref_type=ref, opti=opti)
    evolution_plots(X=X_sol, U=U_sol, X_ref=X_init_ref, U_ref=U_init_ref, dm=dynamic_model, ref_type=ref)
    checking_sum_forces(X=X_sol, X_ref=X_init_ref, U = U_sol, U_ref = U_init_ref, dm=dynamic_model, ref_type=ref, opti=opti)
    print("\nPlots saved under plots folders...\n")
    # Conditionally generate videos
    if args.make_video:
        animate_trajectories_td(full_array, ref_type=ref, label="Actual")
        animate_trajectories(X_sol, ref_type=ref, label="Actual")
        print("\nVideos saved under videos folder ...")
    else:
        print("\nVideo creation skipped (use --make_video to enable).\n")

    print("Exit!")
    