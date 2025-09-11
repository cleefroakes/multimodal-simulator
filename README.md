# multimodal-simulator
Multimodal AI simulator for math and mechanics with matplotlib rendering

## About
This project is a FastAPI-based multimodal AI simulator that solves math and mechanics problems, rendering results with matplotlib. It supports text and image inputs via '/parse' and '/simulate' endpoints, generating plots stored in 'artifacts/'.


## Local Setup
'''bash
conda create -n sim_ai python=3.8
conda activate sim_ai
conda install -c conda-forge matplotlib sympy numpy fatsapi uvicorn pythin-multipart pillow
uvicorn main:app --reload