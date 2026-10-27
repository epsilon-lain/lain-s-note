"""Run the research's exact checks with assertions enabled, from any cwd."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).parent
if not __debug__:
    raise RuntimeError('Assertions are required; do not use Python -O')

for name in ('verify.py', 'wave_compare.py', 'quadratic_barriers.py',
             'pure_root_barriers.py', 'envelope_audit.py', 'weighted_wave.py'):
    print('\nRunning', name, flush=True)
    subprocess.run([sys.executable, str(ROOT / name)], cwd=ROOT, check=True)
print('\nAll exact-check programs completed successfully.')
