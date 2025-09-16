from app.models.spec import ProblemSpec, Variable
import os
from dotenv import load_dotenv
from huggingface_hub import InferenceClient

load_dotenv()
HF_TOKEN = os.getenv("HF_TOKEN")

def parse_image_and_text(image_path: str, text: str) -> ProblemSpec:
    if not HF_TOKEN:
        print("Warning: HF_TOKEN not set, using fallback parsing")
        return parse_text(text)
    try:
        client = InferenceClient(token=HF_TOKEN)
        response = client.image_to_text(image_path=image_path, prompt=text, model="llava-hf/llava-13b")
        if "beam" in response.lower():
            return ProblemSpec(
                domain="mechanics",
                objective="Beam deflection",
                variables=[Variable(name="L", value=2.0, units="m"), Variable(name="F", value=1000, units="N")],
                equations=["-EI u'''' = f"]
            )
        else:
            return ProblemSpec(
                domain="math",
                objective="Solve equation",
                variables=[Variable(name="x", value=1.0, units="")],
                equations=[response]
            )
    except Exception as e:
        print(f"VLM error: {str(e)}")
        return parse_text(text)