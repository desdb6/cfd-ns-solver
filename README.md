# ns-cf-solver

A Python library for simulating 2D incompressible Newtonian fluid dynamics, built on finite-difference discretisation of the Navier-Stokes equations.
This project includes the toolbox, a demo script, and an in depth-report on various hydrodynamic problems investigated using the ns-cf-library.

---

## Features

- Finite-difference solver for the 2D incompressible Navier-Stokes equations
- Pressure-velocity coupling via a relaxation scheme
- Configurable boundary conditions: Dirichlet, Neumann, and Periodic
- Rectangular no-slip obstacle support with custom boundary conditions
- Vorticity calculation and symmetrical logarithmic scaling
- Built-in stability checker and adaptive timestep utility
- Rich visualisation suite: quiver, streamline, vorticity, and velocity profile plots

---

## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/desdb6/ns-cf-solver.git
cd ns-cf-solver
pip install numpy matplotlib
```

---

## Quick Start

```python
from ns-cf-solver import NewtonianFluid

# Create a 1x1 domain with dx=dy=0.05, dt=0.001
fluid = NewtonianFluid(
    Lx=1.0, Ly=1.0,
    dx=0.05, dy=0.05,
    dt=0.001,
    mu=0.01, rho=1.0, U=1.0
)

# Lid-driven cavity: moving top wall
fluid.set_boundary_conditions_u(BVC_up=1.0, BVC_down=0, BVC_left=0, BVC_right=0)

# Check stability before running
fluid.check_stability()

# Run the simulation
fluid.simulate_system()

# Visualise results
fluid.plot_subplots()
```

---

## Usage

### Creating a fluid

```python
fluid = NewtonianFluid(
    Lx=2.0, Ly=1.0,    # Domain dimensions
    dx=0.05, dy=0.05,   # Grid spacing
    dt=0.001,           # Timestep
    mu=0.01,            # Dynamic viscosity
    rho=1.0,            # Density
    U=1.0,              # Characteristic velocity (used for Re calculation)
    Fx=0.0, Fy=0.0      # Body forces
)
```

### Boundary conditions

Boundaries can be set independently for each edge. Accepted values are `"Neumann"`, `"Periodic"`, or a numeric Dirichlet value.

```python
fluid.set_boundary_conditions_u(BVC_up=1.0, BVC_down=0, BVC_left=0, BVC_right=0)
fluid.set_boundary_conditions_v(BVC_up=0, BVC_down=0, BVC_left=0, BVC_right=0)
fluid.set_boundary_conditions_p(BVC_up="Neumann", BVC_down="Neumann",
                                 BVC_left="Neumann", BVC_right="Neumann")
```

### Adding obstacles

```python
obstacle = fluid.add_obstacle()
obstacle.set_dimensions(x0=0.4, x1=0.6, y0=0.4, y1=0.6)

# Optional: override default no-slip conditions
obstacle.set_boundary_conditions_u(BVC_up=0, BVC_down=0, BVC_left=0, BVC_right=0)
```

### Running the simulation

```python
fluid.simulate_system(
    tol=5e-4,                    # Convergence tolerance
    max_iteratrions_initial=500, # Pre-solve velocity iterations
    max_iterations=10000         # Maximum total iterations
)
```

### Stability

```python
fluid.check_stability()   # Prints whether current dt satisfies stability criterion
fluid.update_dt(factor=10) # Set dt automatically based on grid spacing and viscosity
```

---

## Visualisation

| Method | Description |
|---|---|
| `plot_subplots()` | 2×2 overview: speed, pressure, streamlines, vorticity |
| `plot_quiver_speed()` | Velocity quiver plot with speed contours |
| `plot_quiver_pressure()` | Velocity quiver plot with pressure contours |
| `plot_streamline()` | Streamlines over vorticity field |
| `plot_vorticity()` | Vorticity contour plot |
| `plot_vorticity_symlog()` | Vorticity on a symmetric log scale |
| `plot_u_profile(x)` | Horizontal velocity profile at a given x |
| `plot_v_profile(y)` | Vertical velocity profile at a given y |
| `plot_profiles(x, y)` | Both profiles side by side |

---

## Example: Lid-Driven Cavity

The classic lid-driven cavity benchmark. The top wall moves at U=1, all other walls are stationary no-slip.

```python
fluid = NewtonianFluid(1.0, 1.0, 0.05, 0.05, 0.001, mu=0.01, rho=1.0, U=1.0)
fluid.set_boundary_conditions_u(BVC_up=1.0, BVC_down=0, BVC_left=0, BVC_right=0)
fluid.simulate_system()
fluid.plot_streamline()
```

## Example: Channel Flow with Obstacle

```python
fluid = NewtonianFluid(2.0, 1.0, 0.04, 0.04, 0.0005, mu=0.01, rho=1.0, U=0.0, Fx=0.1)
fluid.set_boundary_conditions_u(BVC_up=0, BVC_down=0, BVC_left="Neumann", BVC_right="Neumann")

obstacle = fluid.add_obstacle()
obstacle.set_dimensions(x0=0.7, x1=0.9, y0=0.3, y1=0.7)

fluid.simulate_system()
fluid.plot_subplots()
```

---

## Reynolds Number

The Reynolds number is computed on initialisation as:

```
Re = max(Lx, Ly) * U * rho / mu
```

It can be updated post-simulation based on the mean flow speed:

```python
fluid.update_reynolds_number()
print(fluid.reynolds_number)
```

---

## Requirements

- Python 3.8+
- NumPy
- Matplotlib

---

## License

MIT License. See `LICENSE` for details.
