"""Global Gram selection on the same integer chain, requiring q >= 33216."""
from fractions import Fraction as F
from pathlib import Path
from math import gcd
import json
import sympy as s
from integral_certificate import certificate

HERE=Path(__file__).resolve().parent
k=s.symbols('k')

def sf(v):
    v=s.Rational(v)
    return F(int(v.p),int(v.q))

def fl(v):return v.numerator//v.denominator
def ce(v):return -fl(-v)

def run():
    base=json.loads((HERE/'checks.json').read_text())
    data=json.loads((HERE/'inputs.json').read_text())
    qmin=33216;modulus=840
    P0,A0,P1,A1=base['period_pair']
    assert (P0,A0,P1,A1)==(313044,99645,-4948,-1575)
    raw=s.Poly(s.sympify(base['raw_score_polynomial'],locals={'k':k}),k)
    coeff=[sf(v) for v in raw.all_coeffs()]
    roots=s.Poly(s.diff(raw.as_expr(),k),k).intervals(eps=s.Rational(1,16))
    root_boxes=[(sf(a),sf(b)) for (a,b),mult in roots]
    def raw_value(kk):
        value=F(0)
        for c in coeff:value=value*kk+c
        return value
    def actual_score(kk):
        g=gcd(kk-3,840)
        assert gcd(P0+P1*kk,A0+A1*kk)==g
        q=abs(A0+A1*kk)//g
        return raw_value(kk)/g**4,q,g
    per_residue=[];all_ks=set()
    for residue in range(modulus):
        g=gcd(residue-3,modulus)
        js=set()
        # Every root of the derivative is common to all constant-gcd classes.
        for lo,hi in root_boxes:
            low=(lo-residue)/modulus;high=(hi-residue)/modulus
            js.update(range(fl(low)-1,ce(high)+2))
        # |A(k)| >= qmin*g cuts the lattice into two closed rays.
        boundaries=[F(sign*qmin*g-A0,A1) for sign in (-1,1)]
        for boundary in boundaries:
            point=(boundary-residue)/modulus
            js.update(range(fl(point)-1,ce(point)+2))
        valid=[]
        for j in js:
            kk=residue+modulus*j
            score,q,actual_g=actual_score(kk)
            if q>=qmin:
                valid.append((score,kk,q,actual_g));all_ks.add(kk)
        assert valid
        winner=min(valid)
        per_residue.append({'residue':residue,'gcd':g,'candidate_count':len(valid),
                            'best_k':winner[1],'best_q':winner[2],
                            'best_bound_squared':str(winner[0])})
    best=min(per_residue,key=lambda r:(F(r['best_bound_squared']),r['best_k']))
    kk=best['best_k'];score,q,g=actual_score(kk)
    a,b=data['rows'];R=[x+kk*y for x,y in zip(a['R'],b['R'])]
    alpha=int(a['alpha'])+kk*int(b['alpha'])
    rat=F(P0+P1*kk,A0+A1*kk)
    cert=certificate(data['c'],R,alpha,K=180)
    assert cert['q']==str(q) and q>=qmin
    assert F(cert['quality_upper'])**2<score
    initial=base['input_integral_certificates'][0]
    # Numeric displays are not part of any selection or proof comparison.
    is_real_improvement=F(cert['quality_upper'])<F(initial['quality_lower'])
    chosen_vs_initial_lower=F(cert['quality_lower'])/F(initial['quality_upper'])
    # Regression on a window larger than every sampled critical boundary.
    maximum=max(abs(x) for x in all_ks)
    checked=0
    for j in range(-max(20000,maximum+100),max(20000,maximum+100)+1):
        value,den,gg=actual_score(j)
        if den>=qmin:
            assert value>=score
            checked+=1
    output={'status':'PASS','q_min':qmin,'scope':'Exact global Gram minimum on the existing fixed chain under a denominator-growth constraint; no theorem that all other true qualities are worse.',
            'modulus':modulus,'residue_classes':len(per_residue),
            'derivative_root_boxes':[[str(a),str(b)] for a,b in root_boxes],
            'distinct_candidates':len(all_ks),'maximum_candidate_abs_k':maximum,
            'classes':per_residue,'winner':{'k':kk,'p':str(rat.numerator),'q':str(q),'gcd':g,
              'alpha':str(alpha),'R':R,'H':[x+kk*y for x,y in zip(a['H'],b['H'])],
              'bound_squared':str(score),'bound_display':float(score)**.5,
              'integral_certificate':cert},
            'strict_true_Q_improvement_over_initial':is_real_improvement,
            'true_Q_ratio_certified_lower':str(chosen_vs_initial_lower),
            'true_Q_ratio_lower_display':float(chosen_vs_initial_lower),
            'regression_cases':checked,'selection_uses_pi_or_CF':False}
    (HERE/'growth_checks.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'status':'PASS','q_min':qmin,'classes':len(per_residue),
      'candidates':len(all_ks),'k':kk,'rational':str(rat),'gcd':g,
      'Gram_Q_upper':float(score)**.5,'true_Q_upper':cert['quality_upper_display'],
      'true_Q_improved':is_real_improvement,'regression_cases':checked},indent=2))

if __name__=='__main__':run()
