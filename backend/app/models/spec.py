from pydantic import BaseModel
<<<<<<< HEAD
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
=======
from typing import Dict, List, Optional

class Variable(BaseModel):
    name: str
    value: Optional[float] = None
    unit: Optional[str] = None

class ProblemSpec(BaseModel):
    input_text: str
    domain: str
    equations: List[str]
    entities: Optional[Dict] = {}
    math_expressions: Optional[List[str]] = []
    image_data: Optional[Dict] = None
    variables: Optional[List[Variable]] = []
>>>>>>> 2232274 (Initial commit of multimodal AI simulator)
