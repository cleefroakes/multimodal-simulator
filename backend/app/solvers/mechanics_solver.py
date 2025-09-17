from app.models.spec import ProblemSpec
import matplotlib.pyplot as plt
import numpy as np
import os

def solve_beam(spec: ProblemSpec):
    try:
        # Beam parameters (hardcoded for MVP, can be extended from spec later)
        L = 2.0  # Length (m)
        F = 1000  # Load (N)
        E = 200e9  # Young's modulus (Pa)
        I = 1e-6   # Moment of inertia (m^4)

        # Analytical solution for cantilever beam deflection
        x = np.linspace(0, L, 100)
        y = (F * x**2) / (6 * E * I) * (3 * L - x)
        max_deflection = (F * L**3) / (3 * E * I)

        # Plotting
        plt.figure(figsize=(8, 6))
        plt.plot(x, y, label="Deflection u(x)", color="#1e90ff", linewidth=2)
        plt.scatter([L], [max_deflection], color="red", s=100, label=f"Max Deflection: {max_deflection:.6f}m")
        plt.axhline(0, color="black", linestyle="--", alpha=0.5)
        plt.title("Beam Deflection (Analytical)", fontsize=14)
        plt.xlabel("x (m)", fontsize=12)
        plt.ylabel("u (m)", fontsize=12)
        plt.grid(True, linestyle="--", alpha=0.7)
        plt.legend()

        # Save plot
        os.makedirs("artifacts", exist_ok=True)
        plt.savefig("artifacts/beam_deflection.png", dpi=300, bbox_inches="tight")
        plt.close()

        return {
            "solution": "Beam deflection solved analytically",
            "plot": "artifacts/beam_deflection.png",
            "explanation": f"Cantilever beam of {L}m length, {F}N load. Max deflection at free end: {max_deflection:.6f}m."
        }
    except Exception as e:
        return {"error": f"Analytical solution failed: {str(e)}", "message": "Beam simulation error"}