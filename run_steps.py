"""Run completed numerical steps in dependency order from the repository root."""
import argparse,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('steps',nargs='*',type=int,default=list(range(2,11)));a=p.parse_args()
r=Path(__file__).resolve().parent
for n in a.steps:
    print(f'Running step {n}',flush=True)
    subprocess.run([sys.executable,str(r/f'step{n}/step{n}.py')],cwd=r,check=True)
