"""Rerun batch-root search with the main verifier's independent pi enclosure."""
from pathlib import Path
from fractions import Fraction
import sys
import json

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE/'independent'))
from batch_floor_search import run_search
from verify import pi_bounds

lo, hi = pi_bounds()  # 220/75 terms, instead of the original 200/70.
result = run_search(10**100, lo, hi, Fraction(1, 250))
assert not result['unresolved']
assert [(x['p'], x['q'], x['gcd']) for x in result['retained']] == [('355','113','1')]
assert result['statistics']['node_count'] == 345
saved = json.loads((HERE/'independent/batch_floor_search_1e100.json').read_text())
assert [(n['A'],n['B'],n['count'],n['status']) for n in result['nodes']] == [
    (n['A'],n['B'],n['count'],n['status']) for n in saved['nodes']]
out = {'all_checks_passed': True, 'pi_interval_independent_of_original': True,
       'search_input_uses_continued_fractions': False,
       'same_complete_tree_and_counts': True, 'statistics': result['statistics'],
       'retained': result['retained'],
       'note': 'Same floor_sum implementation, independent pi enclosure; the separate CF certificate is an algorithmically independent cross-check.'}
(HERE/'batch_verification.json').write_text(json.dumps(out, indent=2))
print(json.dumps({'passed': True, 'statistics': result['statistics']}))
