from app.utils.parser import parse_text
from app.solvers.math_solver import solve_math
spec = parse_text('dummy')
result = solve_math(spec)
print(result)
