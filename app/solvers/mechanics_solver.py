from app.models.spec import ProblemSpec
import os
import subprocess
import matplotlib.pyplot as plt

def solve_beam(spec: ProblemSpec):
    # Beam params (hardcoded for MVP)
    L = 2.0  # Length (m)
    F = 1000  # Load (N)
    E = 200e9  # Young's modulus (Pa)
    I = 1e-6   # Moment of inertia (m^4)

    # Subprocess to run FEniCSx in Conda env
    env_cmd = ['C:\\Users\\DELL\\miniconda3\\Scripts\\conda.exe', 'run', '-n', 'fenicsx-env', 'python', '-c']
    code = f"""
from dolfinx import fem, mesh
from dolfinx.fem import Function, FunctionSpace, dirichletbc
from dolfinx.fem.petsc import LinearProblem
from dolfinx.mesh import create_interval
from ufl import TrialFunction, TestFunction, dot, dx, grad, Constant
import numpy as np
import matplotlib.pyplot as plt
import os

# Create 1D mesh (serial, no MPI)
domain = create_interval(comm=None, nx=50, points=[0, {L}])
V = FunctionSpace(domain, ('CG', 1))

# Variational form (Euler-Bernoulli beam)
u = TrialFunction(V)
v = TestFunction(V)
E = Constant(domain, float({E}))
I = Constant(domain, float({I}))
f = Constant(domain, float(-{F}))
a = E * I * dot(grad(u), grad(v)) * dx
L = f * v * dx

# Boundary condition (fixed at x=0)
bc = dirichletbc(float(0), fem.locate_dofs_geometrical(V, lambda x: np.isclose(x[0], 0)), V)

# Solve
problem = LinearProblem(a, L, bcs=[bc])
uh = problem.solve()

# Plot deflection
x = np.linspace(0, {L}, 100)
y = [uh(xp) for xp in x]
plt.plot(x, y)
plt.title('Beam Deflection')
plt.xlabel('x (m)')
plt.ylabel('u (m)')
plt.grid(True)
os.makedirs('artifacts', exist_ok=True)
plt.savefig('artifacts/beam_deflection.png')
plt.close()
print('Max deflection:', max(y))
"""

    # Run in Conda env
    try:
        subprocess.run(env_cmd + [code], check=True, cwd=os.getcwd(), shell=True)
    except subprocess.CalledProcessError as e:
        return {"error": f"FEniCSx failed: {str(e)}", "message": "Beam sim error"}

    return {
        "solution": "Beam deflection solved with FEniCSx",
        "plot": "artifacts/beam_deflection.png",
        "explanation": f"Cantilever beam of {L}m, {F}N load. Max deflection at end."
    }