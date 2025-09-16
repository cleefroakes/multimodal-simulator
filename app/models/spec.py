from pydantic import BaseModel
from typing import List, Optional

class Variable(BaseModel):
    name: str
    value: float
    units: str

class ProblemSpec(BaseModel):
    domain: str
    objective: str
    inputs: Optional[List[str]] = []
    variables: List[Variable] = []
    equations: List[str] = []
    geometry: Optional[str] = None
    boundary_conditions: List[str] = []
    deliverables: List[str] = ["plot", "report"]