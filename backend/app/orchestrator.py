from app.models.spec import ProblemSpec
from app.solvers.math_solver import solve_math
from app.solvers.mechanics_solver import solve_beam

def orchestrate(spec: ProblemSpec):
    if spec.domain == "math" and spec.equations:
        return solve_math(spec)
    elif spec.domain == "mechanics":
        return solve_beam(spec)
    return {"error": "Invalid domain or no equations"}