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
    label2unit_dict = {
        "p_x": "x (m)",
        "p_y": "y (m)",
        "p_z": "z (m)",
        "v_x": "v (m/s)",
        "v_y": "v (m/s)",
        "v_z": "v (m/s)",
        "L_x": "L_x (kg*m^2/s)",
        "L_y": "L_y (kg*m^2/s)",
        "L_z": "L_z (kg*m^2/s)",
        "rf_z": "z (m)",
        "lf_z": "z (m)",
    }
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
    #add units to y-axis label
    plt.ylabel(f"{label2unit_dict.get(labels[0], '')}")
    plt.legend()
    plt.grid(True)
    save_path = os.path.join("./plots/", f"{title}.png")  # Save as PNG (change extension if needed)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")  # High-quality save
    return

def plot_components(full_array, phases_duration, ref):
    # com trajectory
    p_x = full_array[:, 0]
    p_y = full_array[:, 1]
    p_z = full_array[:, 2]
    _plt(phases_duration, [p_x], ["p_x"], f"CoM x Trajectory {ref}")
    _plt(phases_duration, [p_y], ["p_y"], f"CoM y Trajectory {ref}")
    _plt(phases_duration, [p_z], ["p_z"], f"CoM z Trajectory {ref}")

    # com velocity
    v_x = full_array[:, 7]
    v_y = full_array[:, 8]
    v_z = full_array[:, 9]
    _plt(phases_duration, [v_x, v_y, v_z], ["v_x", "v_y", "v_z"], f"CoM Velocity {ref}")

    # angular momentum
    L_x = full_array[:,10]
    L_y = full_array[:,11]
    L_z = full_array[:,12]
    _plt(phases_duration, [L_x, L_y, L_z], ["L_x", "L_y", "L_z"], f"Angular Momentum {ref}")

    # Feet positions
    rf_x = full_array[:,14]
    rf_y = full_array[:,15]
    rf_z = full_array[:,16]
    lf_x = full_array[:,17]
    lf_y = full_array[:,18]
    lf_z = full_array[:,19]
    _plt(phases_duration, [rf_z, lf_z], ["rf_z", "lf_z"], f"Feet_z {ref}")

    return

def evolution_contact_forces(X, X_ref, U, U_ref, dm, ref_type, opti):

    N_intervals = X.shape[1]  # Number of intervals
    time = np.arange(N_intervals-1)  # Create a time vector

    translational_fs, rotational_etas = [], []
    ref_translational_fs, ref_rotational_etas = [], []

    SIGMA_L_k = np.array(opti.value(dm.SIGMA_L_k)).astype("int")[:, :-1]

    # Loop to extract the relevant components
    for k in range(N_intervals-1):
        x_k = X[:, k]
        x_ref_k = X_ref[:, k]
        p_k, q_k, v_k, L_k, t_k, P_L_k, Q_L_k = dm.get_x_comp(x_k)
        ref_pos, _, ref_vel, ref_for, _, ref_PLk, _ = dm.get_x_comp(x_ref_k)
        u_k = U[:, k]
        u_ref_k = U_ref[:, k]
        _, _, _, LAMBDA_L_k, R_L_k, ETA_HAT_L_k = dm.get_u_comp(u_k)
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
        "Translational Forces F": (translational_fs, None),
        "Rotational Forces ETA": (rotational_etas, None),
    }

    for i, (title, (values, ref_values)) in enumerate(properties.items()):
        for j in range(3):  # Loop through x, y, z components

            # Right Foot (index 0)
            axs_right[i, j].plot(time, values[:, j, 0], label=f"{title} ({labels[j]}) - Right Foot", color="b")

            if N_intervals > 1:
                axs_right[i, j].set_xticks(time, labels = SIGMA_L_k[0, :], fontsize = 10)
            else:
                axs_right[i, j].set_xticks(time, labels = np.array(SIGMA_L_k), fontsize = 10)
            #axs_left[i, j].set_xticklabels(opti.value(dm.SIGMA_L_k[0, :]))

            if ref_values is not None:
                axs_right[i, j].plot(time, ref_values[:, j, 0], linestyle="dashed", label=f"Ref {title} ({labels[j]}) - Right Foot", color="r")

            #set y label (F or ETA and unit)
            if title == "Translational Forces F":
                axs_right[i, j].set_ylabel(f"F_{labels[j]} (N)")
            elif title == "Rotational Forces ETA":
                axs_right[i, j].set_ylabel(f"ETA_{labels[j]} (N*m)")
            axs_right[i, j].set_title(f"{title} - {labels[j]} (Right Foot)")
            axs_right[i, j].legend()
            axs_right[i, j].grid()

            # Left Foot (index 1)
            axs_left[i, j].plot(time, values[:, j, 1], label=f"{title} ({labels[j]}) - Left Foot", color="g")
            if N_intervals > 1:
                axs_left[i, j].set_xticks(time, labels = SIGMA_L_k[1, :], fontsize = 10)
            else:
                axs_left[i, j].set_xticks(time, labels = np.array(SIGMA_L_k), fontsize = 10)
            if ref_values is not None:
                axs_left[i, j].plot(time, ref_values[:, j, 1], linestyle="dashed", label=f"Ref {title} ({labels[j]}) - Left Foot", color="orange")

            #set y label (F or ETA and unit)
            if title == "Translational Forces F":
                axs_left[i, j].set_ylabel(f"F_{labels[j]} (N)")
            elif title == "Rotational Forces ETA":
                axs_left[i, j].set_ylabel(f"ETA_{labels[j]} (N*m)")
            axs_left[i, j].set_title(f"{title} - {labels[j]} (Left Foot)")
            axs_left[i, j].legend()
            axs_left[i, j].grid()

    plt.tight_layout()
    save_path_left = os.path.join("./plots/", f"contact_forces_{ref_type}_left.png")  # Save as PNG (change extension if needed)
    save_path_right = os.path.join("./plots/", f"contact_forces_{ref_type}_right.png")  # Save as PNG (change extension if needed)
    
    fig_left.savefig(save_path_left, dpi = 300, bbox_inches="tight")
    fig_right.savefig(save_path_right, dpi = 300, bbox_inches="tight")

    #plt.savefig(save_path, dpi=300, bbox_inches="tight")  # High-quality save
    
    return

def evolution_plots(X, U, X_ref, U_ref, dm, ref_type):
    
    N_intervals = X.shape[1] # Number of intervals
    time = np.arange(N_intervals)  # Create a time vector

    # Extract components for plotting
    positions, velocities, momentum, foot_positions_L, foot_velocities_L = [], [], [], [], []
    ref_position, ref_velocity, ref_momentum, ref_P_L_k, ref_foot_velocities_L = [], [], [], [], []
    
    # Loop to extract the relevant components
    for k in range(N_intervals):
        x_k = X[:, k]
        x_ref_k = X_ref[:, k]
        p_k, q_k, v_k, L_k, t_k, P_L_k, Q_L_k = dm.get_x_comp(x_k)   
        ref_pos, _, ref_vel, ref_for, _, ref_PLk, _ = dm.get_x_comp(x_ref_k)
        if k < U.shape[1]:
            u_k = U[:,k]
            _, V_L_k, _, _, _, _ = dm.get_u_comp(u_k)
            u_k_ref = U_ref[:,k]
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
            if title == "Angular Momentum":
                ax_row[j].set_ylabel(f"L_{labels[j]} (kg*m^2/s)")
            elif title == "CoM Position":
                ax_row[j].set_ylabel(f"{labels[j]} (m)")
            elif title == "CoM Velocity":
                ax_row[j].set_ylabel(f"v_{labels[j]} (m/s)")
            ax_row[j].set_title(f"{title} - {labels[j]}")
            ax_row[j].legend()
            ax_row[j].grid()

    # Foot Positions (P_L_k) for Left and Right Foot separately
    for j in range(3):  # Loop through x, y, z
        axs[3, j].plot(time, foot_positions_L[:, j, 0], label=f"Foot 1 {labels[j]}")
        axs[3, j].plot(time, foot_positions_L[:, j, 1], label=f"Foot 2 {labels[j]}")
        axs[3, j].plot(time, ref_P_L_k[:, j, 0], color='r', linestyle='--', label=f"ref_Foot 1 {labels[j]}")
        axs[3, j].plot(time, ref_P_L_k[:, j, 1], color='g', linestyle='--', label=f"ref_Foot 2 {labels[j]}")
        axs[3, j].set_ylabel(f"{labels[j]} (m)")
        axs[3, j].set_title(f"Foot Position - {labels[j]}")
        axs[3, j].legend()
        axs[3, j].grid()

        axs[4, j].plot(time, foot_velocities_L[:, j, 0], label=f"Foot 1 {labels[j]}")
        axs[4, j].plot(time, foot_velocities_L[:, j, 1], label=f"Foot 2 {labels[j]}")
        axs[4, j].plot(time, ref_foot_velocities_L[:, j, 0], color='r', linestyle='--', label=f"ref_Foot 1 {labels[j]}")
        axs[4, j].plot(time, ref_foot_velocities_L[:, j, 1], color='g', linestyle='--', label=f"ref_Foot 2 {labels[j]}")
        axs[4, j].set_ylabel(f"v_{labels[j]} (m)")
        axs[4, j].set_title(f"Foot Velocity - {labels[j]}")
        axs[4, j].legend()
        axs[4, j].grid()
        
        
    axs[-1, 0].set_xlabel("Time (steps)")
    axs[-1, 1].set_xlabel("Time (steps)")
    axs[-1, 2].set_xlabel("Time (steps)")
    plt.tight_layout()
    save_path = os.path.join("./plots/", f"contact_x_{ref_type}.png")  # Save as PNG (change extension if needed)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")  # High-quality save
    
    return

def checking_sum_forces(X, X_ref, U, U_ref, dm, ref_type, opti):
    N_intervals = X.shape[1] - 1  # Number of intervals
    SIGMA_L = dm.SIGMA_L_k

    gravity_force = cs.DM(dm.m) * cs.DM(dm.g)

    # Open a file for writing (you can change 'output.txt' to your preferred file name)
    with open(f'./outputs/forces_balance_{ref_type}.txt', 'w') as file:
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

def animate_com_plots(X, X_ref, phases_durations, timestep):
    """
    X: np.ndarray, shape (T, 3)
    X_ref: np.ndarray, shape (3, N)  # Note: shape adjusted to match ref_com = X_ref[i, :]
    phases_durations: np.ndarray, shape (N,) or (1, N), durations in seconds
    timestep: float, duration of one frame in seconds (e.g., 0.01 for 100 Hz)
    """
    coords = ['x', 'y', 'z']

    # Flatten durations if needed and compute when to show each ref point
    durations = np.array(phases_durations).flatten()
    ref_display_frames = np.cumsum(durations) / timestep  # convert seconds to frame indices
    ref_display_frames = ref_display_frames.astype(int)

    for i, coord in enumerate(coords):
        com = X[:, i]                # shape (T,)
        ref_com = X_ref[i, :]        # shape (N,)
        frames = len(com)

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.set_xlim(0, frames)
        min_val = min(np.min(com), np.min(ref_com))
        max_val = max(np.max(com), np.max(ref_com))
        ax.set_ylim(min_val - 0.1 * abs(min_val), max_val + 0.1 * abs(max_val))

        ax.set_xlabel('Time')
        ax.set_ylabel(f'CoM {coord}')
        ax.set_title(f'CoM vs Ref - {coord}')

        com_line, = ax.plot([], [], label='CoM', color='blue')
        ref_scatter = ax.scatter([], [], label='Ref CoM (Phases)', color='red')
        ax.legend()

        def init():
            com_line.set_data([], [])
            ref_scatter.set_offsets(np.empty((0, 2)))
            return com_line, ref_scatter

        def update(frame):
            t = np.arange(frame + 1)
            com_line.set_data(t, com[:frame + 1])

            # Display ref points whose phase starts at or before current frame
            indices_to_show = np.where(ref_display_frames <= frame)[0]
            if len(indices_to_show) > 0:
                scatter_x = ref_display_frames[indices_to_show]
                scatter_y = ref_com[indices_to_show]
                offsets = np.column_stack((scatter_x, scatter_y))
                ref_scatter.set_offsets(offsets)
            else:
                ref_scatter.set_offsets(np.empty((0, 2)))

            return com_line, ref_scatter

        ani = FuncAnimation(fig, update, frames=frames, init_func=init,
                            blit=True, interval=timestep * 1000)

        ani.save(f"./videos/animated_plot_com_{coord}.gif", writer='pillow')
        # plt.close(fig)  # Keep open if debugging

    return

def animate_contact_2d(X, ref_type, label):
    save_path = f"./videos/trajectory_2d_points_{ref_type}.gif"

    # Extract 2D (X, Y) positions
    com_xy = X[:, 0:2]
    left_foot_xy = X[:, 14:16]
    right_foot_xy = X[:, 17:19]

    # Set up the 2D plot
    fig, ax = plt.subplots(figsize=(8, 6))

    com_line, = ax.plot([], [], color='blue', label='CoM Trajectory')
    left_dot, = ax.plot([], [], 'go', label='Left Foot', markersize=8)
    right_dot, = ax.plot([], [], 'ro', label='Right Foot', markersize=8)

    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_title(f'{label} - {ref_type}')
    ax.legend()
    ax.set_xlim(np.min(X[:, [0,14,17]]) - 0.1, np.max(X[:, [0,14,17]]) + 0.1)
    ax.set_ylim(np.min(X[:, [1,15,18]]) - 0.1, np.max(X[:, [1,15,18]]) + 0.1)
    ax.grid(True)

    # Initialize CoM trace
    com_x, com_y = [], []

    def update(frame):
        # Update CoM line
        com_x.append(com_xy[frame, 0])
        com_y.append(com_xy[frame, 1])
        com_line.set_data(com_x, com_y)

        # Update foot positions (as moving dots)
        left_dot.set_data(left_foot_xy[frame, 0], left_foot_xy[frame, 1])
        right_dot.set_data(right_foot_xy[frame, 0], right_foot_xy[frame, 1])

        return com_line, left_dot, right_dot

    # Create animation
    anim = FuncAnimation(fig, update, frames=X.shape[0], interval=100, blit=True)

    # Save the animation
    anim.save(save_path, writer='pillow', fps=30)
    return 


def animate_trajectories_td(X, ref_type, label):
    save_path=f"./videos/trajectory_animation_td_{ref_type}.gif"
    
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


def animate_trajectories(X, ref_type, label=""):
    save_path=f"./videos/trajectory_animation_{ref_type}.gif"
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