from fastapi import FastAPI, Form, File, UploadFile
<<<<<<< HEAD
from typing import Optional
import os
import shutil
from app.utils.parser import parse_text
from app.solvers.math_solver import solve_math
from app.solvers.mechanics_solver import solve_beam
from app.models.spec import ProblemSpec
from app.utils.vlm import parse_image_and_text

app = FastAPI()

@app.get("/")
def root():
    return {"message": "Multimodal Simulation AI MVP - Ready for math and mechanics sims!"}

@app.post("/simulate")
async def simulate(text: str = Form(...), image: Optional[UploadFile] = File(None)):
    # Parse input (image first if provided, else text)
    if image is not None:
        try:
            os.makedirs("uploads", exist_ok=True)
            image_path = f"uploads/{image.filename}"
            with open(image_path, "wb") as buffer:
                shutil.copyfileobj(image.file, buffer)
            print(f"Image saved successfully to {image_path}")  # Debug
            spec = parse_image_and_text(image_path, text)  # VLM
            print(f"VLM spec: {spec.dict()}")  # Debug
        except Exception as e:
            print(f"Image handling error: {str(e)}")  # Log
            return {"error": str(e), "message": "Image upload failed"}
    else:
        spec = parse_text(text)  # Text-only

    # Select solver
    if spec.domain == "mechanics":
        result = solve_beam(spec)
    else:
        result = solve_math(spec)

    return {
        "spec": spec.dict(),
        "result": result,
        "message": f"Simulation complete! Check artifacts/{'beam_deflection.png' if spec.domain == 'mechanics' else 'plot.png'}"
    }

@app.post("/parse")
async def parse_input(text: str = Form(...)):
    spec = parse_text(text)
    return {"spec": spec.dict()}
=======
from app.utils.parser import parse_text
from app.utils.vlm import parse_image_and_text
from app.solvers.math_solver import solve_math
from app.solvers.mechanics_solver import solve_beam

app = FastAPI()

@app.post("/parse")
async def parse(text: str = Form(...)):
    spec = parse_text(text)
    return {"spec": spec}

@app.post("/simulate")
async def simulate(text: str = Form(...)):
    spec = parse_text(text)
    if spec.domain == "math":
        return solve_math(spec)
    elif spec.domain == "mechanics":
        return solve_beam(spec)
    return {"error": "Unknown domain"}

@app.post("/simulate_with_image")
async def simulate_with_image(text: str = Form(...), image: UploadFile = File(...)):
<<<<<<< HEAD
    spec = parse_image_and_text(image.filename, text) # type: ignore
    return solve_beam(spec)
>>>>>>> 2232274 (Initial commit of multimodal AI simulator)
=======
    spec = parse_image_and_text(image. filename, text) # type: ignore
    return solve_beam(spec)
>>>>>>> 434bbdf (Add type ignore to math_solver.py and update main.py for simulator functionality)
