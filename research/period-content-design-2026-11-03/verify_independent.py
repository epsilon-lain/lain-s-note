"""Replay both independent verifiers using only the standard library."""
from pathlib import Path
from verify_weighted_independent import run as verify_candidates
from verify_global_independent import run as verify_global

if __name__=='__main__':
    folder=Path(__file__).resolve().parent
    verify_candidates(folder,folder/'independent_weighted_verification.json')
    verify_global()
