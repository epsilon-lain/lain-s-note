"""Algorithm checks and a q=1009 scale benchmark, with exact replay."""
from fractions import Fraction as Q
from pathlib import Path
import json,time,sys,hashlib
# Long exact orbit points are intentional research output, not decimal input.
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(100000)
from itinerary_certificate import propose,certify_itinerary,verify_certificate
HERE=Path(__file__).resolve().parent

def compact(c):
    return {k:v for k,v in c.items() if k not in ('word','period_point','phase_interval','multiplier')}

def run():
    cases=[('small_half',Q(51,100),Q(1,2),1,2),
           ('small_third',Q(17,50),Q(1,2),1,3),
           ('pi_bridge_known_fraction',Q(1,7),Q(493,5000),16,113)]
    report=[]
    for name,om,kap,p,q in cases:
        start=time.perf_counter();cert=propose(om,kap,1000,20000)
        assert cert is not None and (cert['p'],cert['q'])==(p,q)
        assert verify_certificate(json.loads(json.dumps(cert)))
        (HERE/(name+'_certificate.json')).write_text(json.dumps(cert,indent=2)+'\n')
        report.append({'name':name,**compact(cert),'seconds_display_only':time.perf_counter()-start})
    # Deliberately input a large q to benchmark validation only. This example is
    # NOT an approximation discovery, and does not use the target pi.
    q=1009;start=time.perf_counter()
    cert=certify_itinerary(Q(1,q),Q(1,q**3),1,q,1/(2*q))
    assert cert is not None and cert['q']==q
    assert verify_certificate(json.loads(json.dumps(cert)))
    # Retain only a fingerprint and bit costs of the large intermediate certificate.
    certificate_hash=hashlib.sha256(json.dumps(cert,sort_keys=True).encode()).hexdigest()
    report.append({'name':'q1009_validation_scale_only',**compact(cert),'certificate_sha256':certificate_hash,'seconds_display_only':time.perf_counter()-start})
    # Adversarial modifications are rejected by the float-free verifier.
    bad=json.loads(json.dumps(cert));bad['integer_lift']+=1
    try:verify_certificate(bad)
    except AssertionError:pass
    else:raise AssertionError('corrupted integer lift accepted')
    # A failed proposal is explicitly inconclusive, not a non-locking proof.
    unknown=propose(Q(1416,10000),Q(1,10),1000,20000)
    assert unknown is None
    result={'status':'PASS','certificates':report,'corrupted_certificate_rejected':True,
      'unresolved_proposal_example':{'omega':'177/1250','kappa':'1/10','max_period':1000,
         'scope':'No certificate found by this heuristic; not a proof that no period exists.'},
      'arithmetic':'Only itinerary proposals and runtimes use floats; exact certificates replay without floats or target pi.'}
    (HERE/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':run()
