#!/usr/bin/env python3
"""Reproduce finite checks, source-byte audit, and an explicit verification ledger."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PREVIOUS = HERE.parent/'pi-construction-mechanisms-2026-11-01'
sys.path.insert(0, str(PREVIOUS))
from source_support import audit_source_files

LABS = [
    ('arithmetic_rewrite_lab.py', 'arithmetic_rewrite_results.json'),
    ('selected_columns_lab.py', 'selected_columns_results.json'),
    ('tail_elimination_lab.py', 'tail_elimination_results.json'),
    ('exact_cauchy_lab.py', 'results_cauchy.json'),
    ('exact_budget_lab.py', 'results_adaptive_budget.json'),
    ('collision_transport_lab.py', 'collision_transport_results.json'),
    ('refined_transport_lab.py', 'refined_transport_results.json'),
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    records = []
    for script, output in LABS:
        process = subprocess.run([sys.executable, str(HERE/script)], cwd=ROOT,
                                 text=True, capture_output=True, timeout=120)
        if process.returncode:
            print(process.stdout)
            print(process.stderr, file=sys.stderr)
            raise SystemExit(process.returncode)
        data = json.loads((HERE/output).read_text())
        records.append(dict(script=script, exit_code=process.returncode,
                            output=output, script_sha256=sha(HERE/script),
                            output_sha256=sha(HERE/output), json_valid=True))
        print(json.dumps(dict(script=script, status='passed')), flush=True)
    result = dict(status='passed', date='2026-11-01', labs=records,
                  source_audit=audit_source_files(),
                  scope='Finite algebra, optimization and formula checks only; '
                        'no independent pi theorem, Lean kernel, or endpoint conclusion.',
                  dependencies=dict(python=sys.version.split()[0],
                                    sympy=__import__('sympy').__version__),
                  source_dependency='../pi-construction-mechanisms-2026-11-01/source_support.py')
    (HERE/'verification.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(dict(status='passed', labs=len(records),
                          source_audit=result['source_audit']['status'])))


if __name__ == '__main__':
    main()
