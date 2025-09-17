import sys
from app.orchestrator import orchestrate

if __name__== "__name__":
   #Placeholder fpr Ray/Perfect
   spec = {"domain": sys.argv[1]}  #Simplified for now
   result = orchestrate(spec)
   print(result)