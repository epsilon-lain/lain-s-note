#!/usr/bin/env python3
"""Reproduce exact controls and record their hashes."""
from pathlib import Path
import subprocess,sys,json,hashlib

HERE=Path(__file__).resolve().parent
LABS=[('quantization_lab.py','results_quantization.json'),('mechanism_lab.py','mechanism_results.json')]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    rows=[]
    for script,result in LABS:
        r=subprocess.run([sys.executable,str(HERE/script)],capture_output=True,text=True,timeout=120)
        if r.returncode:
            print(r.stdout);print(r.stderr,file=sys.stderr);raise SystemExit(r.returncode)
        json.loads((HERE/result).read_text())
        rows.append(dict(script=script,result=result,exit_code=0,
                         script_sha256=sha(HERE/script),result_sha256=sha(HERE/result)))
        print(json.dumps(dict(script=script,status='passed')),flush=True)
    r=subprocess.run([sys.executable,str(HERE/'plot_models.py')],capture_output=True,text=True,timeout=120)
    assert r.returncode==0,(r.stdout,r.stderr)
    results=dict(status='passed',date='2026-11-01',labs=rows,
                 figure=dict(path='resolution_controls.png',sha256=sha(HERE/'resolution_controls.png')),
                 scope='Exact finite controls only; no physical field or pi endpoint theorem.',
                 python=sys.version.split()[0])
    (HERE/'verification.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(dict(status='passed',labs=2,figures=1)))

if __name__=='__main__':main()
