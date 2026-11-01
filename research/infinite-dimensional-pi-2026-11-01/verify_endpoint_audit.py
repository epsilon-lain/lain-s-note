"""Run the exact/symbolic/numerical checks embedded in endpoint_audit.ipynb.

Requires Python 3, sympy and mpmath. The checks do not replace the proofs.
"""
from pathlib import Path
import contextlib
import io
import json

def main():
    folder = Path(__file__).resolve().parent
    notebook = json.loads((folder / 'endpoint_audit.ipynb').read_text(encoding='utf-8'))
    namespace = {'__name__': '__audit__'}
    results = []
    for cell in notebook['cells']:
        if cell['cell_type'] != 'code':
            continue
        source = cell['source']
        if isinstance(source, list):
            source = ''.join(source)
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            exec(compile(source, cell.get('id', '<cell>'), 'exec'), namespace)
        results.append({'cell_id': cell.get('id'), 'status': 'passed', 'stdout': output.getvalue()})
    result = {'status': 'passed', 'code_cells': len(results), 'results': results,
              'scope': 'Implementation checks only; general statements rely on proofs in endpoint_audit.md.'}
    (folder / 'endpoint_audit_verification.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(f"Endpoint audit: {len(results)} code cells passed.")

if __name__ == '__main__':
    main()
