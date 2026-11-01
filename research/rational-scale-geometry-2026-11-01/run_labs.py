#!/usr/bin/env python3
"""Run all four exact laboratories and record artifact hashes."""
from pathlib import Path
import subprocess,sys,json,hashlib

HERE=Path(__file__).resolve().parent
LABS=[('integer_geometry_checks.py','integer_geometry_results.json'),
      ('rational_lenses_lab.py','rational_lenses_results.json'),
      ('denominator_envelope_lab.py','denominator_envelope_results.json'),
      ('model_lab.py','model_results.json')]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    checks=[]
    for script,result in LABS:
        r=subprocess.run([sys.executable,str(HERE/script)],capture_output=True,text=True,timeout=120)
        if r.returncode:
            print(r.stdout);print(r.stderr,file=sys.stderr);raise SystemExit(r.returncode)
        json.loads((HERE/result).read_text())
        checks.append(dict(script=script,result=result,exit_code=0,
                           script_sha256=sha(HERE/script),result_sha256=sha(HERE/result)))
        print(json.dumps(dict(script=script,status='passed')),flush=True)
    r=subprocess.run([sys.executable,str(HERE/'plot_models.py')],capture_output=True,text=True,timeout=120)
    assert r.returncode==0,(r.stdout,r.stderr)
    figures=[dict(path=p.name,sha256=sha(p)) for p in sorted(HERE.glob('*.png'))]
    figures += [dict(path=p.name,sha256=sha(p)) for p in sorted(HERE.glob('*.svg'))]
    state=dict(status='passed',date='2026-11-01',labs=checks,figures=figures,
               scope='Finite exact models and symbolic proofs; no infinite-tail conclusion for pi.',
               dependencies=dict(python=sys.version.split()[0],matplotlib=__import__('matplotlib').__version__))
    (HERE/'verification.json').write_text(json.dumps(state,indent=2)+'\n')
    print(json.dumps(dict(status='passed',labs=len(checks),figures=len(figures))))

if __name__=='__main__':main()
