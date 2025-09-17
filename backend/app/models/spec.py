from pydantic import BaseModel
from typing import Dict, List, Optional

class Variable(BaseModel):
    name: str
    value: Optional[float] = None
    unit: Optional[str] = None

class ProblemSpec(BaseModel):
    input_text: Optional[str] = None
    domain: str
    objective: Optional[str] = None
    equations: List[str] = []
    entities: Optional[Dict] = {}
    math_expressions: Optional[List[str]] = []
    image_data: Optional[Dict] = None
    variables: Optional[List[Variable]] = []
    geometry: Optional[str] = None
    boundary_conditions: List[str] = []
    deliverables: List[str] = ["plot", "report"]