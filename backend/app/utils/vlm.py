from app.models.spec import ProblemSpec, Variable
import os
from dotenv import load_dotenv
from huggingface_hub import InferenceClient
import matplotlib.pyplot as plt
from PIL import Image

load_dotenv()
HF_TOKEN = os.getenv("HF_TOKEN")

def parse_image_and_text(image_path: str, text: str) -> ProblemSpec:
    if HF_TOKEN:
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
    
    # Fallback to basic image and text parsing
    variables = [Variable(name="x")] if "x" in text else []
    domain = "mechanics" if "beam" in text.lower() else "math"
    equations = [text] if "=" in text else []
    plot_path = None
    if "beam" in text.lower():
        img = Image.open(image_path)
        plt.figure(figsize=(8, 6))
        plt.imshow(img)
        plt.title(f"Beam Image: {os.path.basename(image_path)}", fontsize=14)
        plt.axis("off")
        os.makedirs("artifacts", exist_ok=True)
        plt.savefig("artifacts/beam_image.png", dpi=300, bbox_inches="tight")
        plt.close()
        plot_path = "artifacts/beam_image.png"
    
    return ProblemSpec(
        input_text=text,
        domain=domain,
        equations=equations,
        entities={},
        math_expressions=[],
        image_data={"path": image_path, "plot": plot_path},
        variables=variables
    )