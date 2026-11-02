"""Run only exact audits for the imported unique-root proof.

Uses Python's standard library; does not invoke the optional mpmath examples.
The imported source and its directory are never written. The JSON report is
saved beside this verifier as newton_exact_check_results.json.

Examples:
    python -B verify_newton_exact.py
    python -B verify_newton_exact.py --source /path/to/verify_unique_root.py
    python -B verify_newton_exact.py --source /path/to/imported-next-work
"""

import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys


def resolve_source(argument):
    if argument is not None:
        supplied = Path(argument).expanduser().resolve()
        candidates = [supplied]
        if supplied.is_dir():
            candidates = [
                supplied / "verify_unique_root.py",
                supplied / "root_geometry" / "verify_unique_root.py",
            ]
    else:
        package_root = Path(__file__).resolve().parent.parent
        candidates = [
            package_root / "imported-next-work" / "root_geometry"
            / "verify_unique_root.py",
            package_root / "joint-input" / "extracted" / "root_geometry"
            / "verify_unique_root.py",
        ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    attempted = ", ".join(str(candidate) for candidate in candidates)
    raise FileNotFoundError("No source verifier found; tried: " + attempted)


def import_source(source):
    # Importing a source file normally may create __pycache__ beside it.
    # Disable bytecode writing before executing the inspected input module.
    sys.dont_write_bytecode = True
    name = "_newton_exact_imported_unique_root"
    spec = importlib.util.spec_from_file_location(name, source)
    if spec is None or spec.loader is None:
        raise ImportError("Could not load source verifier: " + str(source))
    module = importlib.util.module_from_spec(spec)
    # dataclass decorators require the module to be registered during import.
    sys.modules[name] = module
    spec.loader.exec_module(module)
    for attribute in ("exact_audits", "B", "K"):
        if not callable(getattr(module, attribute, None)):
            raise TypeError("Source verifier lacks callable " + attribute)
    return module


def sharpening_audits(module):
    identity_checks = envelope_checks = 0
    for t in (Fraction(6), Fraction(44, 7), Fraction(7)):
        for n in range(4, 31):
            b2 = module.B(n, 2, t)
            assert b2 == module.K(n, t) * module.K(n - 1, t) / (2*n - 2)
            identity_checks += 1
            for k in range(2, n + 1):
                assert module.B(n, k, t) <= b2 * module.K(n - 1, t)**(k - 2)
                envelope_checks += 1
    assert identity_checks == 81
    assert envelope_checks == 1296
    return {
        "B_n2_exact_identity": identity_checks,
        "B_n2_based_geometric_envelope": envelope_checks,
        "rational_t_values": ["6", "44/7", "7"],
        "n_range_inclusive": [4, 30],
        "status": "PASS",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        help="Source verifier file, its root_geometry directory, or its parent directory.",
    )
    arguments = parser.parse_args()
    if not __debug__:
        parser.error("Do not use -O: the exact audits use assert statements.")
    try:
        source = resolve_source(arguments.source)
    except FileNotFoundError as error:
        parser.error(str(error))
    report = Path(__file__).resolve().with_name("newton_exact_check_results.json")
    if report.parent == source.parent:
        parser.error("Keep this verifier outside the input source directory.")

    module = import_source(source)
    result = {
        "status": "PASS",
        "source": str(source),
        "exact_checks": module.exact_audits(),
        "independent_sharpening_checks": sharpening_audits(module),
        "numeric_examples_run": False,
        "scope": "Finite exact algebra checks only; independent analytic review "
                 "supplies asymptotic conclusions.",
    }
    report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
