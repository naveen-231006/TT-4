"""
Design and Analysis of a Tensegrity-Based Tumbleweed Robot
with Controlled Locomotion - Final Simulation v6

Material: Dyneema SK75 cables (UHMWPE fiber)
Geometry: Icosahedral tensegrity (12 nodes, 30 cables)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('default')
plt.rcParams.update({
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'axes.edgecolor': 'black',
    'axes.labelcolor': 'black',
    'text.color': 'black',
    'xtick.color': 'black',
    'ytick.color': 'black',
    'grid.color': '#cccccc',
    'savefig.facecolor': 'white',
})
from mpl_toolkits.mplot3d import Axes3D
import time as timer

# ============================================================
# MATERIAL PROPERTIES - Dyneema SK75 (UHMWPE)
# ============================================================
MATERIAL_NAME = "Dyneema SK75 (UHMWPE)"
E_cable = 100e9          # Young's modulus [Pa]
sigma_allow = 900e6      # Allowable stress [Pa] (UTS ~3.4 GPa, SF~3.8)
d_cable = 0.004          # Cable diameter [m]
A_cable = np.pi * (d_cable/2)**2  # Cross-section area [m^2]
rho_cable = 970           # Density [kg/m^3]

# ============================================================
# ROBOT PARAMETERS
# ============================================================
L_edge = 0.30             # Edge length [m]
m_payload = 0.50          # Central payload mass [kg]
m_node_extra = 0.02       # Extra mass per node (joint hardware) [kg]

# Cable mass
total_cable_length = 30 * L_edge  # 30 cables
m_cables = rho_cable * A_cable * total_cable_length
m_total = m_payload + 12 * m_node_extra + m_cables
m_node = m_total / 12     # Lumped mass per node [kg]

# ============================================================
# STIFFNESS SCALING
# ============================================================
k_phys = E_cable * A_cable / L_edge   # Physical stiffness [N/m]
k_sim = 5000.0                         # Simulation stiffness [N/m]
k_ratio = k_phys / k_sim              # Stress scaling factor

# ============================================================
# PRESTRESS
# ============================================================
PRESTRAIN = 0.00135  # 0.135%

# ============================================================
# SIMULATION PARAMETERS
# ============================================================
g = 9.81
dt = 0.0015             # Time step [s]
T_total = 60.0           # Total simulation time [s]
T_settle = 5.0           # Settling time before wind [s]
N_steps = int(T_total / dt)

# Damping (higher to prevent excessive bouncing)
zeta = 0.15
c_damp = 2 * zeta * np.sqrt(k_sim * m_node)

# Ground contact (penalty method - NO clamping)
k_ground = 50000.0   # Ground spring [N/m]
c_ground = 200.0      # Ground dashpot [N*s/m] (high to prevent bounce)
mu_s = 0.6            # Static friction
mu_k = 0.4            # Kinetic friction
v_max = 5.0           # Velocity cap [m/s] for stability

# Wind / locomotion force
F_wind_total = 0.5    # Total wind force [N] (realistic for small robot)
wind_dir = np.array([1.0, 0.3, 0.0])
wind_dir /= np.linalg.norm(wind_dir)

# ============================================================
# ICOSAHEDRON GEOMETRY
# ============================================================
phi = (1 + np.sqrt(5)) / 2
t_v = L_edge / 2
s_v = L_edge * phi / 2

vertices = np.array([
    [ 0,  t_v, -s_v], [ 0, -t_v, -s_v],
    [ 0,  t_v,  s_v], [ 0, -t_v,  s_v],
    [ t_v, -s_v,  0], [-t_v, -s_v,  0],
    [ t_v,  s_v,  0], [-t_v,  s_v,  0],
    [-s_v,  0, -t_v], [-s_v,  0,  t_v],
    [ s_v,  0, -t_v], [ s_v,  0,  t_v],
], dtype=np.float64)

# Center at origin
vertices -= vertices.mean(axis=0)

# Scale to exact L_edge
dists_all = []
for i in range(12):
    for j in range(i+1, 12):
        dists_all.append(np.linalg.norm(vertices[i] - vertices[j]))
min_dist = min(dists_all)
vertices *= L_edge / min_dist

# Find cables (edges of icosahedron = pairs at distance L_edge)
cables = []
for i in range(12):
    for j in range(i+1, 12):
        if np.linalg.norm(vertices[i] - vertices[j]) < L_edge * 1.05:
            cables.append((i, j))
cables = np.array(cables)
N_cables = len(cables)

# Compute rest lengths with prestrain
L_natural = np.array([np.linalg.norm(vertices[c[0]] - vertices[c[1]]) for c in cables])
L0 = L_natural * (1.0 - PRESTRAIN)

# Lift structure so lowest node is at z = 0.01 (slight gap)
vertices[:, 2] -= vertices[:, 2].min()
vertices[:, 2] += 0.01

# ============================================================
# PRINT SETUP INFO
# ============================================================
print("=" * 65)
print("  TENSEGRITY TUMBLEWEED ROBOT SIMULATION")
print("=" * 65)
print(f"\n--- Material ---")
print(f"  Material:        {MATERIAL_NAME}")
print(f"  Young's Modulus: {E_cable/1e9:.0f} GPa")
print(f"  Cable Diameter:  {d_cable*1000:.1f} mm")
print(f"  Cable Area:      {A_cable*1e6:.4f} mm^2")
print(f"  Allowable Stress:{sigma_allow/1e6:.0f} MPa")
print(f"\n--- Robot ---")
print(f"  Edge Length:     {L_edge*100:.0f} cm")
print(f"  Total Mass:      {m_total:.3f} kg")
print(f"  Mass/Node:       {m_node:.4f} kg")
print(f"  Nodes:           12, Cables: {N_cables}")
print(f"\n--- Stiffness ---")
print(f"  k_phys:          {k_phys:.0f} N/m")
print(f"  k_sim:           {k_sim:.0f} N/m")
print(f"  k_ratio:         {k_ratio:.1f}")
print(f"  Prestrain:       {PRESTRAIN*100:.3f}%")
F_pre_per_cable = k_sim * PRESTRAIN * L_edge
print(f"  F_pre (sim):     {F_pre_per_cable:.3f} N/cable")
sigma_pre = (F_pre_per_cable * k_ratio) / A_cable / 1e6
print(f"  sigma_prestress: {sigma_pre:.1f} MPa")
print(f"\n--- Simulation ---")
print(f"  dt:              {dt*1000:.2f} ms")
print(f"  T_total:         {T_total:.0f} s")
print(f"  N_steps:         {N_steps}")
print(f"  Wind force:      {F_wind_total:.1f} N (after {T_settle:.0f}s)")
print("=" * 65)

# ============================================================
# INITIALIZE STATE
# ============================================================
pos = vertices.copy()
vel = np.zeros((12, 3), dtype=np.float64)

# Pre-allocate recording arrays (sample every 100 steps)
sample_interval = 100
N_samples = N_steps // sample_interval + 1
rec_time = np.zeros(N_samples)
rec_com = np.zeros((N_samples, 3))
rec_KE = np.zeros(N_samples)
rec_PE_grav = np.zeros(N_samples)
rec_PE_elastic = np.zeros(N_samples)
rec_stress_max = np.zeros(N_samples)
rec_stress_mean = np.zeros(N_samples)
rec_min_z = np.zeros(N_samples)
sample_idx = 0

# Cable index arrays for vectorized computation
ci = cables[:, 0]
cj = cables[:, 1]

# ============================================================
# SIMULATION LOOP
# ============================================================
print("\nRunning simulation...")
t_start = timer.time()

for step in range(N_steps):
    t = step * dt

    # --- Cable forces (vectorized) ---
    d_vec = pos[cj] - pos[ci]                    # (N_cables, 3)
    lengths = np.linalg.norm(d_vec, axis=1)       # (N_cables,)
    lengths = np.maximum(lengths, 1e-10)          # avoid /0
    extensions = lengths - L0                      # (N_cables,)
    
    # Tension only (cables can't push)
    F_mag = np.where(extensions > 0, k_sim * extensions, 0.0)
    
    # Damping along cable axis
    d_unit = d_vec / lengths[:, None]             # (N_cables, 3)
    rel_vel = vel[cj] - vel[ci]                   # (N_cables, 3)
    vel_along = np.sum(rel_vel * d_unit, axis=1)  # (N_cables,)
    F_damp = np.where(extensions > 0, c_damp * vel_along, 0.0)
    
    F_total_cable = (F_mag + F_damp)              # (N_cables,)
    F_vec = F_total_cable[:, None] * d_unit       # (N_cables, 3)
    
    # Accumulate forces on nodes
    forces = np.zeros((12, 3), dtype=np.float64)
    np.add.at(forces, ci, F_vec)
    np.add.at(forces, cj, -F_vec)
    
    # --- Gravity ---
    forces[:, 2] -= m_node * g
    
    # --- Ground contact (penalty, per node) ---
    below = pos[:, 2] < 0
    if np.any(below):
        penetration = -pos[below, 2]              # positive when below ground
        # Normal force (spring + dashpot)
        F_normal = k_ground * penetration - c_ground * vel[below, 2]
        F_normal = np.maximum(F_normal, 0.0)      # only push up
        forces[below, 2] += F_normal
        
        # Friction force (Coulomb)
        v_tang = vel[below, :2].copy()            # tangential velocity (x,y)
        v_tang_mag = np.linalg.norm(v_tang, axis=1)
        moving = v_tang_mag > 1e-6
        if np.any(moving):
            F_fric_max = mu_k * F_normal[moving]
            v_dir = v_tang[moving] / v_tang_mag[moving, None]
            # Apply friction opposing motion
            idx_below = np.where(below)[0]
            idx_moving = idx_below[moving]
            forces[idx_moving, :2] -= F_fric_max[:, None] * v_dir
    
    # --- Wind force (after settling) ---
    if t > T_settle:
        F_per_node = F_wind_total / 12.0
        forces[:, 0] += F_per_node * wind_dir[0]
        forces[:, 1] += F_per_node * wind_dir[1]
    
    # --- Symplectic Euler integration ---
    acc = forces / m_node
    vel += acc * dt
    
    # Velocity cap for stability
    v_mag = np.linalg.norm(vel, axis=1)
    too_fast = v_mag > v_max
    if np.any(too_fast):
        vel[too_fast] *= (v_max / v_mag[too_fast, None])
    
    pos += vel * dt
    
    # --- Record data ---
    if step % sample_interval == 0 and sample_idx < N_samples:
        com = pos.mean(axis=0)
        KE = 0.5 * m_node * np.sum(vel**2)
        PE_grav = m_node * g * np.sum(pos[:, 2])
        PE_elastic = 0.5 * k_sim * np.sum(np.where(extensions > 0, extensions**2, 0.0))
        
        # Stress: correct quasi-static approach
        # Prestress is analytical: sigma_pre = E * prestrain (static, exact)
        # Dynamic increment: force in cable is stiffness-independent in
        # quasi-static equilibrium, so sigma_dyn = (F_sim - F_pre_sim) / A
        F_pre_sim = k_sim * PRESTRAIN * L_natural  # sim prestress force per cable
        F_dynamic = np.maximum(F_mag - F_pre_sim, 0.0)  # dynamic force increment
        sigma_pre_val = E_cable * PRESTRAIN / 1e6  # MPa (constant = 135 MPa)
        sigma_dyn = F_dynamic / A_cable / 1e6      # MPa
        stress_phys = np.where(F_mag > 0, sigma_pre_val + sigma_dyn, 0.0)  # MPa
        
        rec_time[sample_idx] = t
        rec_com[sample_idx] = com
        rec_KE[sample_idx] = KE
        rec_PE_grav[sample_idx] = PE_grav
        rec_PE_elastic[sample_idx] = PE_elastic
        rec_stress_max[sample_idx] = stress_phys.max()
        rec_stress_mean[sample_idx] = stress_phys[stress_phys > 0].mean() if np.any(stress_phys > 0) else 0
        rec_min_z[sample_idx] = pos[:, 2].min()
        sample_idx += 1
    
    # --- Progress ---
    if step % (N_steps // 10) == 0:
        com = pos.mean(axis=0)
        print(f"  t={t:6.2f}s  CoM=({com[0]:+.3f}, {com[1]:+.3f}, {com[2]:+.3f})  "
              f"z_min={pos[:,2].min():+.4f}")

t_elapsed = timer.time() - t_start
rec_time = rec_time[:sample_idx]
rec_com = rec_com[:sample_idx]
rec_KE = rec_KE[:sample_idx]
rec_PE_grav = rec_PE_grav[:sample_idx]
rec_PE_elastic = rec_PE_elastic[:sample_idx]
rec_stress_max = rec_stress_max[:sample_idx]
rec_stress_mean = rec_stress_mean[:sample_idx]
rec_min_z = rec_min_z[:sample_idx]

print(f"\nSimulation completed in {t_elapsed:.1f}s ({N_steps/t_elapsed:.0f} steps/s)")

# ============================================================
# FINAL ANALYSIS
# ============================================================
print("\n" + "=" * 65)
print("  RESULTS SUMMARY")
print("=" * 65)

com_final = rec_com[-1]
com_init = rec_com[0]
displacement = np.linalg.norm(com_final[:2] - com_init[:2])
max_stress = rec_stress_max.max()
SF = sigma_allow / 1e6 / max_stress if max_stress > 0 else float('inf')

print(f"\n--- Material Used ---")
print(f"  Cable Material:    {MATERIAL_NAME}")
print(f"  Fiber Type:        Ultra-High-Molecular-Weight Polyethylene (UHMWPE)")
print(f"  Young's Modulus:   {E_cable/1e9:.0f} GPa")
print(f"  Tensile Strength:  ~3,400 MPa (fiber)")
print(f"  Allowable Stress:  {sigma_allow/1e6:.0f} MPa")
print(f"  Density:           {rho_cable} kg/m^3")
print(f"  Cable Diameter:    {d_cable*1000:.1f} mm")
print(f"  Cable Area:        {A_cable*1e6:.4f} mm^2")

print(f"\n--- Structural Performance ---")
print(f"  Prestress (phys):  {sigma_pre:.1f} MPa")
print(f"  Max Dynamic Stress:{max_stress:.1f} MPa")
print(f"  Safety Factor:     {SF:.2f}")

print(f"\n--- Locomotion ---")
print(f"  CoM Initial:       ({com_init[0]:.3f}, {com_init[1]:.3f}, {com_init[2]:.3f}) m")
print(f"  CoM Final:         ({com_final[0]:.3f}, {com_final[1]:.3f}, {com_final[2]:.3f}) m")
print(f"  Lateral Displace:  {displacement:.3f} m")
print(f"  Avg Speed:         {displacement/(T_total-T_settle):.4f} m/s")
print(f"  Wind Force:        {F_wind_total:.1f} N")

print(f"\n--- Energy ---")
print(f"  Final KE:          {rec_KE[-1]:.4f} J")
print(f"  Final PE_grav:     {rec_PE_grav[-1]:.4f} J")
print(f"  Final PE_elastic:  {rec_PE_elastic[-1]:.4f} J")
print("=" * 65)

# ============================================================
# PLOTS
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.patch.set_facecolor('white')
for ax_row in axes:
    for ax in ax_row:
        ax.set_facecolor('white')
fig.suptitle('Tensegrity Tumbleweed Robot — Simulation Results\n'
             f'Material: {MATERIAL_NAME} | L={L_edge*100:.0f}cm | m={m_total:.2f}kg',
             fontsize=13, fontweight='bold')

# 1. Trajectory (top-down)
ax = axes[0, 0]
ax.plot(rec_com[:, 0], rec_com[:, 1], 'b-', linewidth=1, alpha=0.7)
ax.plot(rec_com[0, 0], rec_com[0, 1], 'go', markersize=10, label='Start')
ax.plot(rec_com[-1, 0], rec_com[-1, 1], 'r*', markersize=12, label='End')
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_title('CoM Trajectory (Top View)')
ax.legend()
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)

# 2. CoM height
ax = axes[0, 1]
ax.plot(rec_time, rec_com[:, 2], 'b-', linewidth=1)
ax.plot(rec_time, rec_min_z, 'r--', linewidth=0.8, alpha=0.6, label='Min node z')
ax.axhline(y=0, color='k', linewidth=0.5, linestyle='--')
ax.axvline(x=T_settle, color='orange', linewidth=1, linestyle='--', label=f'Wind ON (t={T_settle}s)')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Height [m]')
ax.set_title('CoM Height & Min Node Z')
ax.legend()
ax.grid(True, alpha=0.3)

# 3. Energy
ax = axes[1, 0]
ax.plot(rec_time, rec_KE, label='Kinetic', linewidth=1)
ax.plot(rec_time, rec_PE_grav, label='PE Gravity', linewidth=1)
ax.plot(rec_time, rec_PE_elastic, label='PE Elastic', linewidth=1)
E_total = rec_KE + rec_PE_grav + rec_PE_elastic
ax.plot(rec_time, E_total, 'k--', label='Total', linewidth=1, alpha=0.6)
ax.axvline(x=T_settle, color='orange', linewidth=1, linestyle='--')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Energy [J]')
ax.set_title('Energy Components')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# 4. Cable stress
ax = axes[1, 1]
ax.plot(rec_time, rec_stress_max, 'r-', linewidth=1, label='Max Cable Stress')
ax.plot(rec_time, rec_stress_mean, 'b-', linewidth=0.8, alpha=0.7, label='Mean Cable Stress')
ax.axhline(y=sigma_allow/1e6, color='k', linewidth=1.5, linestyle='--',
           label=f'Allowable ({sigma_allow/1e6:.0f} MPa)')
ax.axhline(y=sigma_pre, color='g', linewidth=1, linestyle=':', label=f'Prestress ({sigma_pre:.1f} MPa)')
ax.axvline(x=T_settle, color='orange', linewidth=1, linestyle='--')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Stress [MPa]')
ax.set_title('Cable Stress Timeline')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('tensegrity_results.png', dpi=150, bbox_inches='tight', facecolor='white')
print("\nPlots saved to: tensegrity_results.png")

# ============================================================
# 3D SNAPSHOT
# ============================================================
fig3d = plt.figure(figsize=(10, 8))
fig3d.patch.set_facecolor('white')
ax3d = fig3d.add_subplot(111, projection='3d')
ax3d.set_facecolor('white')
ax3d.xaxis.pane.fill = True
ax3d.yaxis.pane.fill = True
ax3d.zaxis.pane.fill = True
ax3d.xaxis.pane.set_facecolor('white')
ax3d.yaxis.pane.set_facecolor('white')
ax3d.zaxis.pane.set_facecolor('white')

# Draw cables
for c in cables:
    p = pos[c]
    ax3d.plot3D(p[:, 0], p[:, 1], p[:, 2], 'b-', linewidth=1.5, alpha=0.7)

# Draw nodes
ax3d.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c='red', s=60, zorder=5)

# Draw ground plane
xlim = [pos[:, 0].min() - 0.1, pos[:, 0].max() + 0.1]
ylim = [pos[:, 1].min() - 0.1, pos[:, 1].max() + 0.1]
xx, yy = np.meshgrid(np.linspace(xlim[0], xlim[1], 2),
                      np.linspace(ylim[0], ylim[1], 2))
ax3d.plot_surface(xx, yy, np.zeros_like(xx), alpha=0.15, color='green')

ax3d.set_xlabel('X [m]')
ax3d.set_ylabel('Y [m]')
ax3d.set_zlabel('Z [m]')
ax3d.set_title(f'Final Configuration (t={T_total}s)\n{MATERIAL_NAME}')

plt.savefig('tensegrity_3d.png', dpi=150, bbox_inches='tight', facecolor='white')
print("3D snapshot saved to: tensegrity_3d.png")

print("\nSimulation complete.")
