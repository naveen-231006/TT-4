import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patheffects as pe
import numpy as np

# ── Figure setup ──
fig, ax = plt.subplots(figsize=(14, 10))
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
ax.set_xlim(0, 14)
ax.set_ylim(0, 11)
ax.axis('off')
fig.patch.set_facecolor('white')

# ── Professional Color Palette ──
C_BG_LOOP   = '#F0F2F5'
C_PHYSICS   = '#E8F5E9'   # soft green
C_PHYSICS_E = '#388E3C'
C_ESTIM     = '#FFF8E1'   # soft amber
C_ESTIM_E   = '#F57F17'
C_CTRL      = '#E3F2FD'   # soft blue
C_CTRL_E    = '#1565C0'
C_ACT       = '#FCE4EC'   # soft pink
C_ACT_E     = '#C62828'
C_IO        = '#ECEFF1'   # grey
C_IO_E      = '#455A64'
C_ARROW     = '#1A1A1A'
C_LABEL     = '#000000'

# ── Helper: draw a rounded box ──
def draw_block(cx, cy, w, h, label, sub, fc, ec, fontsize=12):
    x, y = cx - w/2, cy - h/2
    box = patches.FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.15",
        facecolor=fc, edgecolor=ec, linewidth=2.5,
        zorder=3
    )
    ax.add_patch(box)
    # Main label — bold, black
    ax.text(cx, cy + 0.15, label, ha='center', va='center',
            fontsize=fontsize, fontweight='bold', color='#000000', zorder=4)
    # Sub label (file name) — dark, bold-italic for readability at small sizes
    ax.text(cx, cy - 0.25, sub, ha='center', va='center',
            fontsize=10, fontstyle='italic', fontweight='bold', color='#1A1A1A', zorder=4)
    return cx, cy, w, h

# ── Helper: draw arrow with label placed OUTSIDE the arrow path ──
def draw_arrow(x1, y1, x2, y2, label="", label_side="above", offset=0.18, rad=0.0, dashed=False):
    style = patches.ArrowStyle('->', head_length=8, head_width=5)
    ls = '--' if dashed else '-'
    arrow = patches.FancyArrowPatch(
        (x1, y1), (x2, y2),
        connectionstyle=f"arc3,rad={rad}",
        arrowstyle=style, color=C_ARROW, lw=2.2,
        mutation_scale=14, zorder=2, linestyle=ls
    )
    ax.add_patch(arrow)
    if label:
        mx = (x1 + x2) / 2
        my = (y1 + y2) / 2
        # Place label above or below the midpoint
        if label_side == "above":
            ly = my + offset
        elif label_side == "below":
            ly = my - offset
        elif label_side == "left":
            mx -= offset + 0.6
            ly = my
        elif label_side == "right":
            mx += offset + 0.6
            ly = my
        else:
            ly = my + offset
        ax.text(mx, ly, label, ha='center', va='center',
                fontsize=9.5, fontweight='bold', color='#000000', zorder=5,
                bbox=dict(facecolor='white', edgecolor='#757575',
                          boxstyle='round,pad=0.25', alpha=0.95, linewidth=1.0))

# ══════════════════════════════════════════════════
#  LAYOUT — 4 rows, carefully spaced
# ══════════════════════════════════════════════════

# Row positions (Y)
ROW1 = 9.5   # Init / Vis (outside loop)
ROW2 = 7.2   # Physics row
ROW3 = 4.8   # Estimation + Control row
ROW4 = 2.4   # Actuation row

# Column positions (X)
COL1 = 3.0
COL2 = 7.0
COL3 = 11.0

BW = 2.8  # block width
BH = 1.1  # block height

# ── Background: Loop boundary ──
loop_box = patches.FancyBboxPatch(
    (0.8, 1.3), 12.4, 7.4,
    boxstyle="round,pad=0.3",
    facecolor=C_BG_LOOP, edgecolor='#90A4AE',
    linewidth=2.5, linestyle=(0, (8, 4)), zorder=1
)
ax.add_patch(loop_box)
ax.text(1.2, 8.5, "Simulation Loop", fontsize=14,
        fontweight='bold', color='#37474F', zorder=2,
        fontstyle='italic')

# ── ROW 1: Initialization & Visualization ──
draw_block(COL1, ROW1, BW, BH, "Initialization", "(Geometry, Materials, Actuators)", C_IO, C_IO_E)
draw_block(COL3, ROW1, BW, BH, "Post-Processing", "(3D Render, Stress, Energy Plots)", C_IO, C_IO_E)

# ── ROW 2: Physics ──
draw_block(COL1, ROW2, BW, BH, "ODE Solver", "(ode15s — stiff integrator)", C_PHYSICS, C_PHYSICS_E)
draw_block(COL2, ROW2, BW, BH, "System Dynamics", "(dynamics_fn.m)", C_PHYSICS, C_PHYSICS_E)
draw_block(COL3, ROW2, BW, BH, "Collision & Contact", "(seg_dist_fn.m)", C_PHYSICS, C_PHYSICS_E)

# ── ROW 3: Estimation & Control ──
draw_block(COL1, ROW3, BW, BH, "Sensor Model", "(Noisy Measurements)", C_ESTIM, C_ESTIM_E)
draw_block(COL2, ROW3, BW, BH, "State Estimator (UKF)", "(ukf_pred.m)", C_ESTIM, C_ESTIM_E)
draw_block(COL3, ROW3, BW, BH, "Nonlinear MPC", "(nmpc_solve.m)", C_CTRL, C_CTRL_E)

# ── ROW 4: Actuation ──
draw_block(COL2, ROW4, BW, BH, "Actuation Model", "(Cable Rest Length Control)", C_ACT, C_ACT_E)

# ══════════════════════════════════════════════════
#  ARROWS — all labels placed clearly outside paths
# ══════════════════════════════════════════════════

# Init → ODE (vertical down)
draw_arrow(COL1, ROW1 - BH/2, COL1, ROW2 + BH/2, "Initial State (x₀)", "right", 0.0)

# ODE ↔ Dynamics (horizontal, two arrows offset vertically)
draw_arrow(COL1 + BW/2, ROW2 + 0.15, COL2 - BW/2, ROW2 + 0.15, "State (x)", "above")
draw_arrow(COL2 - BW/2, ROW2 - 0.15, COL1 + BW/2, ROW2 - 0.15, "Derivatives (dx/dt)", "below")

# Dynamics → Collision (horizontal)
draw_arrow(COL2 + BW/2, ROW2, COL3 - BW/2, ROW2, "Distance Queries", "above")

# ODE → Sensor (vertical down)
draw_arrow(COL1, ROW2 - BH/2, COL1, ROW3 + BH/2, "True State + Noise", "right", 0.0)

# Sensor → UKF (horizontal)
draw_arrow(COL1 + BW/2, ROW3, COL2 - BW/2, ROW3, "z (measurements)", "above")

# UKF → NMPC (horizontal)
draw_arrow(COL2 + BW/2, ROW3, COL3 - BW/2, ROW3, "Estimated State (x̂)", "above")

# NMPC → Actuator (vertical down)
draw_arrow(COL3, ROW3 - BH/2, COL3, ROW4 + 0.55,
           "Optimal Control (u*)", "right", 0.0)
# Bend NMPC→Actuator to go down then left
draw_arrow(COL3, ROW4 + 0.55, COL2 + BW/2, ROW4, "", rad=0.0)

# Actuator → ODE (curved left + up)
ax.annotate("", xy=(COL1 - BW/2 - 0.1, ROW2 - 0.2),
            xytext=(COL2 - BW/2, ROW4),
            arrowprops=dict(arrowstyle='->', color=C_ARROW, lw=2.2,
                            connectionstyle='arc3,rad=0.5'),
            zorder=2)
ax.text(1.0, ROW4 + 0.9, "Actuation\nForces",
        ha='center', va='center', fontsize=9.5, fontweight='bold', color='#000000',
        bbox=dict(facecolor='white', edgecolor='#757575',
                  boxstyle='round,pad=0.25', alpha=0.95, linewidth=1.0),
        zorder=5)

# Dynamics → Post-Processing (up-right, dashed)
draw_arrow(COL2, ROW2 + BH/2, COL3, ROW1 - BH/2,
           "Trajectory & Force Logs", "left", 0.0, rad=-0.2, dashed=True)

# ── Legend ──
legend_items = [
    (C_IO,      C_IO_E,      "Input / Output"),
    (C_PHYSICS, C_PHYSICS_E, "Physics Engine"),
    (C_ESTIM,   C_ESTIM_E,   "State Estimation"),
    (C_CTRL,    C_CTRL_E,    "Control"),
    (C_ACT,     C_ACT_E,     "Actuation"),
]
lx, ly = 11.5, 1.8
ax.text(lx + 0.6, ly + 0.3, "Legend", fontsize=11, fontweight='bold',
        ha='center', color='#000000')
for i, (fc, ec, name) in enumerate(legend_items):
    yy = ly - i * 0.35
    box = patches.FancyBboxPatch((lx, yy - 0.1), 0.35, 0.25,
                                  boxstyle="round,pad=0.03",
                                  facecolor=fc, edgecolor=ec, linewidth=1.5, zorder=3)
    ax.add_patch(box)
    ax.text(lx + 0.55, yy + 0.02, name, fontsize=9.5, fontweight='bold', va='center', color='#000000', zorder=4)

# ── Title ──
ax.text(7, 10.6, "Tensegrity Tumbleweed Robot — System Architecture",
        ha='center', va='center', fontsize=18, fontweight='bold', color='#000000',
        path_effects=[pe.withStroke(linewidth=3, foreground='white')])

plt.tight_layout(pad=0.5)
plt.savefig('architecture_image.png', dpi=400, bbox_inches='tight', facecolor='white')
print("Saved architecture_image.png")
