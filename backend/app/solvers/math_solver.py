from sympy import symbols, Eq, solve, sympify
import matplotlib.pyplot as plt
import numpy as np
from app.models.spec import ProblemSpec
import os

def solve_math(spec: ProblemSpec):
    try:
        x = symbols('x')
        eq_str = spec.equations[0]
        if '=' in eq_str:
            left, right = eq_str.split('=', 1)
            left_expr = sympify(left.strip())
            right_expr = sympify(right.strip())
            eq = Eq(left_expr, right_expr)
        else:
            eq = Eq(sympify(eq_str), 0)
        
        solution = solve(eq, x)
        
        plot_path = None
        if 'x**2' in eq_str or any(c in eq_str for c in ['+', '-', '*', '/']):
            xs = np.linspace(-5, 5, 100)
            expr = left_expr - right_expr if '=' in eq_str else sympify(eq_str)
            ys = [float(expr.subs(x, val)) for val in xs]
            os.makedirs("artifacts", exist_ok=True)
            plt.figure(figsize=(8, 6))
            plt.plot(xs, ys, label=f"f(x) = {eq_str}", color="#1e90ff", linewidth=2)
            for root in solution:
                plt.scatter([float(root)], [0], color="red", s=100, label=f"Root x={root}")
            plt.axhline(0, color="black", linestyle="--", alpha=0.5)
            plt.axvline(0, color="black", linestyle="--", alpha=0.5)
            plt.title(f"Equation Plot: {eq_str}", fontsize=14)
            plt.xlabel("x", fontsize=12)
            plt.ylabel("f(x)", fontsize=12)
            plt.grid(True, linestyle="--", alpha=0.7)
            plt.legend()
            plt.savefig("artifacts/plot.png", dpi=300, bbox_inches="tight")
            plt.close()
            plot_path = "artifacts/plot.png"
        
        return {
            "solution": [float(sol) for sol in solution],
            "plot": plot_path,
            "explanation": f"Solved {eq_str} to get roots {solution}"
        }
    except Exception as e:
        return {"error": str(e), "message": "Math solver failed"}