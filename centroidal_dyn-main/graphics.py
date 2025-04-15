import numpy as np # type: ignore
import matplotlib.pyplot as plt
import casadi as cs # type: ignore
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

def print_solutions(solutions):
    x = ["p_x", "p_y", "p_z", "q_1", "q_2", "q_3", "q_4", "v_x", "v_y", "v_z", "L_x", "L_y", "L_z",
     "t", "p_right_x", "p_right_y", "p_right_z", "p_left_x", "p_left_y", "p_left_z", "Q_right_1",
    "Q_right_2", "Q_right_3", "Q_right_4", "Q_left_1", "Q_left_2", "Q_left_3", "Q_left_4"]

    u = ["tau", "v_right_x", "v_right_y", "v_right_z", "v_left_x", "v_left_y","v_left_z",
    "w_right_x", "w_right_y", "w_right_z", "w_left_x", "w_left_y", "w_left_z", "lambda_right",
    "lambda_left", "r_right_x", "r_right_y", "r_right_z", "r_left_x", "r_left_y", "r_left_z",
    "eta_right_x", "eta_right_y", "eta_right_z", "eta_left_x", "eta_left_y", "eta_left_z"]

    for sol in solutions:
        X = sol[0]  # NumPy array (28 x N)
        U = sol[1]  # NumPy array (27 x N)
        N = X.shape[1]  # Number of time steps (columns in X)

        print("\nSolution:")
        for i in range(N):  # Loop over columns (time steps)

            print(f"\n--- Step {i} ---")
            print(f"{'State':<20} | {'Value':<10} || {'Input':<20} | {'Value':<10}")
            print("-" * 65)

            for j in range(X.shape[0]):  # Loop over max of states/inputs
                state_str = f"{x[j]:<20} | {X[j, i]:<10.4f}"
                if N-1 > 1:
                    if i<U.shape[1] and j < 27:
                        input_str = f"{u[j]:<20} | {U[j, i]:<10.4f}" if j < U.shape[0] else ""
                    else:
                        input_str = ""
                else:
                    if i<U.shape[1] and j < 27:
                        input_str = f"{u[j]:<20} | {U[j]:<10.4f}"
                    else:
                        input_str = ""
                print(state_str + " || " + input_str)

            print("-" * 65)  # Line separator after each step
    return


def _plt(phases_duration, components, labels, title):
    time_axis = np.cumsum([0] + phases_duration)[:-1]  # Start at 0
    plt.figure(figsize=(6, 4))
    for i in range(len(components)):
        plt.plot(range(len(components[i])), components[i], label=labels[i])

    phase_start = 0
    for i, duration in enumerate(phases_duration):
        phase_end = phase_start + duration
        plt.axvspan(phase_start, phase_end, color="gray", alpha=0.2 if i % 2 == 0 else 0.1)
        phase_start = phase_end

    plt.xlabel("Time Steps")
    plt.ylabel("Values")
    plt.title(title)
    plt.legend()
    plt.grid(True)
    #plt.show()
    save_path = os.path.join("./plots/", f"{title}.png")  # Save as PNG (change extension if needed)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")  # High-quality save
    return

def plot_components(full_array, phases_duration):
    # com trajectory
    p_x = full_array[:, 0]
    p_y = full_array[:, 1]
    p_z = full_array[:, 2]
    _plt(phases_duration, [p_x], ["p_x"], "CoM x Trajectory")
    _plt(phases_duration, [p_y], ["p_y"], "CoM y Trajectory")
    _plt(phases_duration, [p_z], ["p_z"], "CoM z Trajectory")

    # com velocity
    v_x = full_array[:, 7]
    v_y = full_array[:, 8]
    v_z = full_array[:, 9]
    _plt(phases_duration, [v_x, v_y, v_z], ["v_x", "v_y", "v_z"], "CoM Velocity")

    # angular momentum
    L_x = full_array[:,10]
    L_y = full_array[:,11]
    L_z = full_array[:,12]
    _plt(phases_duration, [L_x, L_y, L_z], ["L_x", "L_y", "L_z"], "Angular Momentum")

    # Feet positions
    rf_x = full_array[:,14]
    rf_y = full_array[:,15]
    rf_z = full_array[:,16]
    lf_x = full_array[:,17]
    lf_y = full_array[:,18]
    lf_z = full_array[:,19]
    _plt(phases_duration, [rf_z, lf_z], ["rf_z", "lf_z"], "Feet_z")

    return

def evolution_contact_forces(X, X_ref, U, U_ref, dm, ref_type, opti):

    N_intervals = X.shape[1]  # Number of intervals
    time = np.arange(N_intervals)  # Create a time vector

    translational_fs, rotational_etas = [], []
    ref_translational_fs, ref_rotational_etas = [], []

    SIGMA_L_k = np.array(opti.value(dm.SIGMA_L_k)).astype("int")

    # Loop to extract the relevant components
    for k in range(N_intervals):
        x_k = X[:, k]
        x_ref_k = X_ref[:, k]
        if N_intervals > 1 and k<U.shape[1]:
            u_ref_k = U_ref[:,k]
            u_k = U[:, k]
        elif N_intervals > 1 and k==U.shape[1]:
            u_ref_k = U_ref[:,k-1]
            u_k = U[:, k-1]
        elif N_intervals == 1:
            u_ref_k = U_ref[:]
            u_k = U[:]

        # Decompose state and control
        p_k, q_k, v_k, L_k, t_k, P_L_k, Q_L_k = dm.get_x_comp(x_k)
        _, _, _, LAMBDA_L_k, R_L_k, ETA_HAT_L_k = dm.get_u_comp(u_k)
        ref_pos, _, ref_vel, ref_for, _, ref_PLk, _ = dm.get_x_comp(x_ref_k)
        _, _, _, ref_LAMBDALk, ref_RLk, ref_ETAHATLk = dm.get_u_comp(u_ref_k)

        F_L_k, ETA_L_k = dm.define_contact_wrench(p_k, LAMBDA_L_k, P_L_k, R_L_k, ETA_HAT_L_k)
        F_L_k_ref, ETA_L_k_ref = dm.define_contact_wrench(ref_pos, ref_LAMBDALk, ref_PLk, ref_RLk, ref_ETAHATLk)

        translational_fs.append(F_L_k)
        rotational_etas.append(ETA_L_k)
        ref_translational_fs.append(F_L_k_ref)
        ref_rotational_etas.append(ref_ETAHATLk)


    # Convert lists to numpy arrays
    translational_fs = np.array(translational_fs)
    rotational_etas = np.array(rotational_etas)

    ref_translational_fs = np.array(ref_translational_fs)
    ref_rotational_etas = np.array(ref_rotational_etas)


    fig_left, axs_left = plt.subplots(2, 3, figsize=(15, 10))  # Left foot plots
    fig_right, axs_right = plt.subplots(2, 3, figsize=(15, 10))  # Right foot plots

    labels = ['x', 'y', 'z']

    properties = {
        "Translational Forces F": (translational_fs, ref_translational_fs),
        "Rotational Forces ETA": (rotational_etas, ref_rotational_etas),
    }

    time = np.arange(N_intervals)  # Fix the shape issue

    for i, (title, (values, ref_values)) in enumerate(properties.items()):
        for j in range(3):  # Loop through x, y, z components

            # Right Foot (index 0)
            axs_left[i, j].plot(time, values[:, j, 0], label=f"{title} ({labels[j]}) - Right Foot", color="b")

            if N_intervals > 1:
                axs_left[i, j].set_xticks(time, labels = SIGMA_L_k[0, :], fontsize = 10)
            else:
                axs_left[i, j].set_xticks(time, labels = np.array(SIGMA_L_k), fontsize = 10)
            #axs_left[i, j].set_xticklabels(opti.value(dm.SIGMA_L_k[0, :]))

            if ref_values is not None:
                axs_left[i, j].plot(time, ref_values[:, j, 0], linestyle="dashed", label=f"Ref {title} ({labels[j]}) - Left Foot", color="r")


            axs_left[i, j].set_ylabel(f"{title} ({labels[j]})")
            axs_left[i, j].set_title(f"{title} - {labels[j]} (Right Foot)")
            axs_left[i, j].legend()
            axs_left[i, j].grid()

            # Left Foot (index 1)
            axs_right[i, j].plot(time, values[:, j, 1], label=f"{title} ({labels[j]}) - Right Foot", color="g")
            if N_intervals > 1:
                axs_right[i, j].set_xticks(time, labels = SIGMA_L_k[1, :], fontsize = 10)
            else:
                axs_left[i, j].set_xticks(time, labels = np.array(SIGMA_L_k), fontsize = 10)
            if ref_values is not None:
                axs_right[i, j].plot(time, ref_values[:, j, 1], linestyle="dashed", label=f"Ref {title} ({labels[j]}) - Left Foot", color="orange")

            axs_right[i, j].set_ylabel(f"{title} ({labels[j]})")
            axs_right[i, j].set_title(f"{title} - {labels[j]} (Left Foot)")
            axs_right[i, j].legend()
            axs_right[i, j].grid()

    plt.tight_layout()
    save_path = os.path.join("./plots/", "contact_forces.png")  # Save as PNG (change extension if needed)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")  # High-quality save
    
    return

def evolution_plots(X, U, X_ref, U_ref, dm, ref_type):
    
    N_intervals = X.shape[1] # Number of intervals
    N = N_intervals-1
    time = np.arange(N_intervals)  # Create a time vector

    # Extract components for plotting
    positions, velocities, momentum, foot_positions_L, foot_velocities_L = [], [], [], [], []
    ref_position, ref_velocity, ref_momentum, ref_P_L_k, ref_foot_velocities_L = [], [], [], [], []
    
    # Loop to extract the relevant components
    for k in range(N_intervals):
        x_k = X[:, k]
        x_ref_k = X_ref[:, k]
        if N > 1 and N<U.shape[1]:
            u_k = U[:,k]
            u_k_ref = U_ref[:,k]
        elif N == 1:
            u_k = U[:]
            u_k_ref = U_ref[:]
        else:
            u_k = U[:,k-1]
            u_k_ref = U_ref[:,k-1]
        
        # Decompose state and control
        p_k, q_k, v_k, L_k, t_k, P_L_k, Q_L_k = dm.get_x_comp(x_k)
        _, V_L_k, _, _, _, _ = dm.get_u_comp(u_k)
            
        ref_pos, _, ref_vel, ref_for, _, ref_PLk, _ = dm.get_x_comp(x_ref_k)
        _, ref_V_L_k, _, _, _, _ = dm.get_u_comp(u_k_ref)
        
        
        positions.append(p_k)
        velocities.append(v_k)
        momentum.append(L_k)
        foot_positions_L.append(P_L_k)
        foot_velocities_L.append(V_L_k)
        
        ref_position.append(ref_pos)
        ref_velocity.append(ref_vel)
        ref_momentum.append(ref_for)
        ref_P_L_k.append(ref_PLk)
        ref_foot_velocities_L.append(ref_V_L_k)
        
    r = 5
    # Convert lists to numpy arrays
    positions = np.round(np.array(positions),r)
    velocities = np.round(np.array(velocities),r)
    momentum = np.round(np.array(momentum),r)
    foot_positions_L = np.round(np.array(foot_positions_L),r)  # Shape (N_intervals, 3, 2)
    foot_velocities_L = np.round(np.array(foot_velocities_L),r)
    ref_position = np.round(np.array(ref_position),r)  
    ref_velocity = np.round(np.array(ref_velocity),r)
    ref_momentum = np.round(np.array(ref_momentum),r)
    ref_P_L_k = np.round(np.array(ref_P_L_k),r)
    ref_foot_velocities_L = np.round(np.array(ref_foot_velocities_L),r)
    
    fig, axs = plt.subplots(5, 3, figsize=(15, 16), sharex=True)

    labels = ['x', 'y', 'z']
    properties = {
        "CoM Position": (positions, ref_position, axs[0]),
        "CoM Velocity": (velocities, ref_velocity, axs[1]),
        "Angular Momentum": (momentum, ref_momentum, axs[2])
    }

    for i, (title, (values, ref_values, ax_row)) in enumerate(properties.items()):
        for j in range(3):  # Loop through x, y, z
            ax_row[j].plot(time, values[:, j], label=f"{title}_{labels[j]}")
            if ref_values is not None:
                ax_row[j].plot(time, ref_values[:, j], color='r', linestyle='--', label=f"ref_{title}_{labels[j]}")
            ax_row[j].set_ylabel(f"{title} ({labels[j]})")
            ax_row[j].set_title(f"{title} - {labels[j]}")
            ax_row[j].legend()
            ax_row[j].grid()

    # Foot Positions (P_L_k) for Left and Right Foot separately
    for j in range(3):  # Loop through x, y, z
        axs[3, j].plot(time, foot_positions_L[:, j, 0], label=f"Foot 1 {labels[j]}")
        axs[3, j].plot(time, foot_positions_L[:, j, 1], label=f"Foot 2 {labels[j]}")
        axs[3, j].plot(time, ref_P_L_k[:, j, 0], color='r', linestyle='--', label=f"ref_Foot 1 {labels[j]}")
        axs[3, j].plot(time, ref_P_L_k[:, j, 1], color='g', linestyle='--', label=f"ref_Foot 2 {labels[j]}")
        axs[3, j].set_ylabel(f"Foot Position {labels[j]} (m)")
        axs[3, j].set_title(f"Foot Position - {labels[j]}")
        axs[3, j].legend()
        axs[3, j].grid()

        axs[4, j].plot(time, foot_velocities_L[:, j, 0], label=f"Foot 1 {labels[j]}")
        axs[4, j].plot(time, foot_velocities_L[:, j, 1], label=f"Foot 2 {labels[j]}")
        axs[4, j].plot(time, ref_foot_velocities_L[:, j, 0], color='r', linestyle='--', label=f"ref_Foot 1 {labels[j]}")
        axs[4, j].plot(time, ref_foot_velocities_L[:, j, 1], color='g', linestyle='--', label=f"ref_Foot 2 {labels[j]}")
        axs[4, j].set_ylabel(f"Foot Velocity {labels[j]} (m)")
        axs[4, j].set_title(f"Foot Velocity - {labels[j]}")
        axs[4, j].legend()
        axs[4, j].grid()
        
        
    axs[-1, 0].set_xlabel("Time (steps)")
    axs[-1, 1].set_xlabel("Time (steps)")
    axs[-1, 2].set_xlabel("Time (steps)")
    plt.tight_layout()
    save_path = os.path.join("./plots/", "contact_x.png")  # Save as PNG (change extension if needed)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")  # High-quality save
    
    return

def checking_sum_forces(X, X_ref, U, U_ref, dm, ref_type, opti):
    N_intervals = X.shape[1] - 1  # Number of intervals
    SIGMA_L = dm.SIGMA_L_k

    gravity_force = cs.DM(dm.m) * cs.DM(dm.g)

    # Open a file for writing (you can change 'output.txt' to your preferred file name)
    with open('./outputs/forces_balance.txt', 'w') as file:
        # Write gravity force at the beginning
        file.write(f"gravity_force: {gravity_force}\n")

        # Loop to extract the relevant components
        for k in range(N_intervals):
            
            # Decompose state and control
            x_k = X[:, k]
            x_ref_k = X_ref[:, k]
            if N_intervals > 1:
                u_k = U[:, k]
                u_ref_k = U_ref[:, k]
                SIGMA_L_k = SIGMA_L[:, k]
            else:
                u_k = U[:]
                u_ref_k = U_ref[:]
                SIGMA_L_k = SIGMA_L

            p_k, q_k, v_k, L_k, t_k, P_L_k, Q_L_k = dm.get_x_comp(x_k)
            _, _, _, LAMBDA_L_k, R_L_k, ETA_HAT_L_k = dm.get_u_comp(u_k)

            F_L_k, ETA_L_k = dm.define_contact_wrench(p_k, LAMBDA_L_k, P_L_k, R_L_k, ETA_HAT_L_k)

            # Calculate forces for right and left feet
            right_foot_x = F_L_k[0, 0]
            left_foot_x = F_L_k[0, 1]
            right_foot_y = F_L_k[1, 0]
            left_foot_y = F_L_k[1, 1]
            right_foot_z = F_L_k[2, 0]
            left_foot_z = F_L_k[2, 1]

            # Sum forces (just a placeholder for your logic)
            sum_forces = cs.sum2(F_L_k)
            SIGMA_L_k_value = opti.value(SIGMA_L_k)

            # Write the data to the file
            file.write(f"{k}-th interval\n")
            file.write(f"right foot x: {right_foot_x} and left foot x: {left_foot_x}\n")
            file.write(f"right foot y: {right_foot_y} and left foot y: {left_foot_y}\n")
            file.write(f"right foot z: {right_foot_z} and left foot z: {left_foot_z}\n")
            file.write(f"SIGMA_L_k: {SIGMA_L_k_value} and sum_forces: {sum_forces}\n")
            file.write("\n")  # Add a blank line for readability

    
    return



def animate_trajectories_td(X, n_e, save_path="trajectory_animation_td.gif", ref_type="", label=""):
    
    
    # Extract positions (assuming X has shape [529, 28])
    com_positions = X[:, 0:3]  # CoM (Center of Mass) positions
    left_foot_positions = X[:, 14:17]  # Left foot positions
    right_foot_positions = X[:, 17:20]  # Right foot positions

    # Set up the figure and 3D axis
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Initialize trajectory lines
    com_line, = ax.plot([], [], [], label="CoM", color="blue", linestyle="-")
    right_foot_line, = ax.plot([], [], [], label="Right Foot", color="red", linestyle="-")
    left_foot_line, = ax.plot([], [], [], label="Left Foot", color="green", linestyle="-")

    # Add markers at the starting positions
    ax.scatter(*com_positions[0], color="blue", s=100, label="CoM Start", edgecolor="black", marker="o", alpha=0.7)
    ax.scatter(*(right_foot_positions[0] + [0.01, 0, 0]), color="red", s=100, label="Right Foot Start", edgecolor="black", marker="s", alpha=0.7)
    ax.scatter(*(left_foot_positions[0] - [0.01, 0, 0]), color="green", s=100, label="Left Foot Start", edgecolor="black", marker="^", alpha=0.7)

    
    # Labels and title
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_zlabel("Z (m)")
    ax.set_title(f"{label} {ref_type}")
    ax.legend()

    # Animation update function
    def update(frame):
        com_line.set_data(com_positions[:frame+1, 0], com_positions[:frame+1, 1])
        com_line.set_3d_properties(com_positions[:frame+1, 2])

        right_foot_line.set_data(right_foot_positions[:frame+1, 0], right_foot_positions[:frame+1, 1])
        right_foot_line.set_3d_properties(right_foot_positions[:frame+1, 2])

        left_foot_line.set_data(left_foot_positions[:frame+1, 0], left_foot_positions[:frame+1, 1])
        left_foot_line.set_3d_properties(left_foot_positions[:frame+1, 2])

        return com_line, right_foot_line, left_foot_line

    # Create animation
    anim = FuncAnimation(fig, update, frames=X.shape[0], interval=200, blit=False)

    # Save animation
    anim.save(save_path, writer='pillow', fps=30)

    return


def animate_trajectories(X, n_e, save_path="trajectory_animation.gif", ref_type="", label=""):
    # Extract positions directly from X
    com_positions = X[0:3, :]  # CoM positions over time
    left_foot_positions = X[14:17, :]  # Left foot positions over time
    right_foot_positions = X[17:20, :]  # Right foot positions over time

    # Set up the figure and 3D axis
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Initialize lines for CoM and feet
    com_line, = ax.plot([], [], [], label="CoM", color="blue", linestyle="-")
    right_foot_line, = ax.plot([], [], [], label="Right Foot", color="red", linestyle="-")
    left_foot_line, = ax.plot([], [], [], label="Left Foot", color="green", linestyle="-")

    # Add markers at the starting positions with slight offsets
    ax.scatter(com_positions[0, 0], com_positions[1, 0], com_positions[2, 0],
               color="blue", s=100, label="CoM Start",
               edgecolor="black", alpha=0.7, marker="o")  # Circle

    ax.scatter(right_foot_positions[0, 0] + 0.01, right_foot_positions[1, 0], right_foot_positions[2, 0],
               color="red", s=100, label="Right Foot Start",
               edgecolor="black", alpha=0.7, marker="s")  # Square

    ax.scatter(left_foot_positions[0, 0] - 0.01, left_foot_positions[1, 0], left_foot_positions[2, 0],
               color="green", s=100, label="Left Foot Start",
               edgecolor="black", alpha=0.7, marker="^")  # Triangle

    # Set axis limits (adjust as needed)
    ax.set_zlim([0, np.max(X[2, :]) + 0.5])
    ax.set_xlim([0, np.max(X[0, :]) + 0.5])

    # Set labels, title, and legend
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_zlabel("Z (m)")
    ax.set_title(f"{label} {ref_type}")
    ax.legend()

    # Animation update function
    def update(frame):
        # Update lines with positions from X
        com_line.set_data(com_positions[0, :frame+1], com_positions[1, :frame+1])
        com_line.set_3d_properties(com_positions[2, :frame+1])

        right_foot_line.set_data(right_foot_positions[0, :frame+1], right_foot_positions[1, :frame+1])
        right_foot_line.set_3d_properties(right_foot_positions[2, :frame+1])

        left_foot_line.set_data(left_foot_positions[0, :frame+1], left_foot_positions[1, :frame+1])
        left_foot_line.set_3d_properties(left_foot_positions[2, :frame+1])

        return com_line, right_foot_line, left_foot_line

    # Create the animation
    anim = FuncAnimation(fig, update, frames=X.shape[1]-1, interval=200, blit=False)

    # Save the animation as a video
    anim.save(save_path, writer="pillow", fps=10)
    
    return