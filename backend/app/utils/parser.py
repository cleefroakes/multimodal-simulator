from app.models.spec import ProblemSpec, Variable

def parse_text(text: str) -> ProblemSpec:
    if "beam" in text.lower():
        return ProblemSpec(
            domain="mechanics",
            objective="Beam deflection",
            variables=[Variable(name="L", value=2.0, units="m"), Variable(name="F", value=1000, units="N")],
            equations=["-EI u'''' = f"]
        )
    elif any(c in text for c in ["x", "=", "+", "-"]):
        return ProblemSpec(
            domain="math",
            objective="Solve equation",
            variables=[Variable(name="x", value=1.0, units="")],
            equations=["x**2 + 2*x + 1 = 0"] if "=" not in text else [text]
        )
    else:
        equations = [text] if "=" in text else []
        variables = [Variable(name="x")] if "x" in text else []
        return ProblemSpec(
            input_text=text,
            domain="math" if any(c in text for c in ["x", "=", "+", "-"]) else "mechanics",
            equations=equations,
            entities={},
            math_expressions=[],
            image_data=None,
            variables=variables
        )