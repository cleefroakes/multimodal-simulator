from app.models.spec import ProblemSpec
import sympy as sp
import matplotlib.pyplot as plt
import numpy as np
import os

def solve_math(spec: ProblemSpec):
    try:
        x = sp.Symbol('x')
        equation = spec.equations[0]
        expr = sp.sympify(equation.split('=')[0]) - sp.sympify(equation.split('=')[1])
        solutions = sp.solve(expr, x)
        
        # Plot
        x_vals = np.linspace(-5, 5, 100)
        f = sp.lambdify(x, expr.lhs, 'numpy')
        y_vals = f(x_vals)

        os.makedirs("artifacts", exist_ok=True)
        plt.plot(x_vals, y_vals)
        plt.axhline(0, color='black', linewidth=0.5)
        plt.title(f"Solving {equation}")
        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid(True)
        plt.savefig("artifacts/plot.png")
        plt.close()
        
        return {
            "solution": [float(sol) for sol in solutions],
            "plot": "artifacts/plot.png",
            "explanation": f"Solved {equation} to get roots {solutions}"
        }
    except Exception as e:
        return {"error": str(e), "message": "Math solver failed"}