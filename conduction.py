import numpy as np
import matplotlib.pyplot as plt

def solve_exp4a():
    print("\n2D Steady State Heat Conduction Solver (Fixed Boundaries)")
    
    # Inputs
    try:
        L = float(input("Enter plate length L (m): "))
        n_div = int(input("Enter number of divisions (n_div): "))
        t_left = float(input("Enter Left boundary temperature (K): "))
        t_right = float(input("Enter Right boundary temperature (K): "))
        t_bottom = float(input("Enter Bottom boundary temperature (K): "))
        t_top = float(input("Enter Top boundary temperature (K): "))
    except ValueError:
        print("Invalid input.")
        return

    # Parameters
    n_nodes = n_div + 1
    total_nodes = n_nodes**2
    dx = L / n_div
    dy = dx # Assuming square grid
    
    print(f"\nComputed Parameters:")
    print(f"dx = dy = {dx:g}")
    print(f"Grid size = {n_nodes} x {n_nodes}")
    print(f"Total equations = {total_nodes}\n")

    # Build Matrix A and Vector B
    A = np.zeros((total_nodes, total_nodes))
    B = np.zeros(total_nodes)

    # Global index helper
    def get_idx(i, j): return i * n_nodes + j

    # Apply Governing Equation and Boundaries
    for i in range(n_nodes):
        for j in range(n_nodes):
            idx = get_idx(i, j)
            # Find bound mapping using inline conditions, evaluate to interior (None) if not an edge
            if i == 0 and j == 0:
                bound = (t_bottom + t_left) / 2.0
            elif i == 0 and j == n_nodes - 1:
                bound = (t_bottom + t_right) / 2.0
            elif i == n_nodes - 1 and j == 0:
                bound = (t_top + t_left) / 2.0
            elif i == n_nodes - 1 and j == n_nodes - 1:
                bound = (t_top + t_right) / 2.0
            elif i == 0:
                bound = t_bottom
            elif i == n_nodes - 1:
                bound = t_top
            elif j == 0:
                bound = t_left
            elif j == n_nodes - 1:
                bound = t_right
            else:
                bound = None
            
            if bound is not None:
                A[idx, idx], B[idx] = 1, bound
            else:
                A[idx, idx], B[idx] = -4, 0
                for di, dj in [(1,0), (-1,0), (0,1), (0,-1)]:
                    A[idx, get_idx(i+di, j+dj)] = 1

    # Print Equations Table
    print("Generated Equation Table")
    print(f"{'NODE':<12} | {'EQUATION':<65} | {'RHS':<12}")
    for i in range(n_nodes):
        for j in range(n_nodes):
            node_str = f"({j},{i})"
            if i == 0 and j == 0:
                eq, rhs = f"1*T({j},{i})", (t_bottom + t_left) / 2.0
            elif i == 0 and j == n_nodes - 1:
                eq, rhs = f"1*T({j},{i})", (t_bottom + t_right) / 2.0
            elif i == n_nodes - 1 and j == 0:
                eq, rhs = f"1*T({j},{i})", (t_top + t_left) / 2.0
            elif i == n_nodes - 1 and j == n_nodes - 1:
                eq, rhs = f"1*T({j},{i})", (t_top + t_right) / 2.0
            elif i == 0:
                eq, rhs = f"1*T({j},{i})", t_bottom
            elif i == n_nodes - 1:
                eq, rhs = f"1*T({j},{i})", t_top
            elif j == 0:
                eq, rhs = f"1*T({j},{i})", t_left
            elif j == n_nodes - 1:
                eq, rhs = f"1*T({j},{i})", t_right
            else:
                eq = f"-4*T({j},{i}) + T({j+1},{i}) + T({j-1},{i}) + T({j},{i+1}) + T({j},{i-1})"
                rhs = 0
            
            print(f"{node_str:<12} | {eq:<65} | {rhs:<12.2f}")

    # Solve
    T_flat = np.linalg.solve(A, B)
    T = T_flat.reshape((n_nodes, n_nodes))

    # Output specific results (Internal nodes only)
    print("\nFinal Results (Internal Nodal Temperatures):")
    for i in range(1, n_nodes - 1):
        for j in range(1, n_nodes - 1):
            print(f"Node ({j},{i}): {T[i,j]:.2f} K")

    # Plotting
    x_vals = np.linspace(0, L, n_nodes)
    y_vals = np.linspace(0, L, n_nodes)

    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 12
    plt.figure(figsize=(10, 8))
    plt.contourf(x_vals, y_vals, T, cmap='hot', levels=20)
    plt.colorbar(label='Temperature (K)')
    
    # Add grid lines to match nodes
    plt.xticks(x_vals)
    plt.yticks(y_vals)
    plt.grid(color='black', linestyle='-', linewidth=1, alpha=0.5)

    for i in range(n_nodes):
        for j in range(n_nodes):
            ha = 'left' if j == 0 else 'right' if j == n_nodes - 1 else 'center'
            va = 'bottom' if i == 0 else 'top' if i == n_nodes - 1 else 'center'
            xo, yo = (5 if j == 0 else -5 if j == n_nodes - 1 else 0), (5 if i == 0 else -5 if i == n_nodes - 1 else 0)
            plt.annotate(f"({j},{i})\n{T[i,j]:.0f} K", xy=(x_vals[j], y_vals[i]), xytext=(xo, yo),
                         textcoords='offset points', ha=ha, va=va, fontsize=9,
                         bbox=dict(facecolor='white', alpha=0.7, edgecolor='black', boxstyle='round,pad=0.3'))

    plt.title('2D Temperature Distribution (Fixed Boundaries)', weight='bold')
    plt.xlabel('Position x (m)', weight='bold')
    plt.ylabel('Position y (m)', weight='bold')
    plt.show()

if __name__ == "__main__":
    solve_exp4a()
    