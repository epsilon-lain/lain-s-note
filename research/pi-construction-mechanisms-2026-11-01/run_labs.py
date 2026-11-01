"""Run all six independent finite labs without redirecting self-written JSON."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
JOBS=[('interpolation_lab.py','interpolation_lab_results.json',False),
      ('period_lab.py','results_period_lab.json',False),
      ('integral_lab.py','integral_lab_results.json',False),
      ('parameter_lab.py','results_parameters.json',True),
      ('check_classical_models.py','classical_model_checks.json',True),
      ('core_models.py','results_core_models.json',True)]


def run(job):
    script,output,stdout_json=job
    r=subprocess.run([sys.executable,str(HERE/script)],cwd=ROOT,
                     text=True,capture_output=True,timeout=300)
    if r.returncode:
        raise RuntimeError(f'{script} exited {r.returncode}: {r.stderr[-2000:]}')
    if stdout_json:
        data=json.loads(r.stdout)
        (HERE/output).write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    data=json.loads((HERE/output).read_text())
    b=(HERE/output).read_bytes()
    return dict(script=script,status='passed',exit_code=r.returncode,output=output,
                output_sha256=hashlib.sha256(b).hexdigest())


def main():
    results=[]
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures=[pool.submit(run,j) for j in JOBS]
        for f in as_completed(futures):
            result=f.result();results.append(result)
            print(result['script']+': passed',flush=True)
    from source_support import audit_source_files
    manifest=dict(status='passed',labs=sorted(results,key=lambda x:x['script']),
                  source_audit=audit_source_files(),
                  scope='Finite exact identities and stated numerical checks; no independent certification of mu(pi)=2 or pi BA/nonBA.')
    (HERE/'verification.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print('All six labs and JSON result files passed.')


if __name__=='__main__':
    main()
