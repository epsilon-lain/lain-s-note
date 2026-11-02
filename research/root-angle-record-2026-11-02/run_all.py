"""Reproduce this directory's finite certificates using the Python standard library.

This does not reproduce the latest 653.52-billion-term author computation.
Run without -O: assertions are part of the verification procedure.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def main():
    if not __debug__:
        raise SystemExit('Run without -O: assertions are required.')
    commands = [
        'root-family/independent/build_cf_certificate.py',
        'root-family/independent/batch_floor_search.py',
        'root-family/verify.py',
        'root-family/verify_batch.py',
        'angle-growth/rotation/verify_angle_chain.py',
        'angle-growth/angle_kernel/trace_angle.py',
        'angle-growth/verify_independent.py',
        'record-review/verify_first_improvement.py',
        'record-review/verify_record_identities.py',
        'record-review/verify_static_formulas.py',
    ]
    checks = []
    for relative in commands:
        print('Checking ' + relative, flush=True)
        process = subprocess.run([sys.executable, str(ROOT / relative)],
                                 cwd=ROOT, text=True, capture_output=True)
        checks.append({'script': relative, 'returncode': process.returncode,
                       'stdout': process.stdout.replace(str(ROOT), '.'),
                       'stderr': process.stderr.replace(str(ROOT), '.')})
        if process.returncode:
            (ROOT / 'VERIFICATION.json').write_text(
                json.dumps({'passed': False, 'checks': checks}, indent=2), encoding='utf-8')
            print(process.stdout)
            print(process.stderr, file=sys.stderr)
            raise SystemExit(process.returncode)
    root = json.loads((ROOT / 'root-family/batch_verification.json').read_text())
    angle = json.loads((ROOT / 'angle-growth/independent_verification.json').read_text())
    record = json.loads((ROOT / 'record-review/first_improvement.json').read_text())
    source_hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in sorted(ROOT.rglob('*.py'))}
    output = {'passed': True, 'python': sys.version, 'checks': checks,
              'source_sha256': source_hashes,
              'root_batch_statistics': root['statistics'],
              'angle_independent': angle,
              'first_record_q_digits': record['first_better_convergent']['q_digits'],
              'scope': 'Finite exact certificates only; no Lean proof, no full latest-record rerun.'}
    (ROOT / 'VERIFICATION.json').write_text(json.dumps(output, ensure_ascii=False, indent=2)+'\n',
                                          encoding='utf-8')
    print('All 10 verification commands passed. See VERIFICATION.json.')


if __name__ == '__main__':
    main()
