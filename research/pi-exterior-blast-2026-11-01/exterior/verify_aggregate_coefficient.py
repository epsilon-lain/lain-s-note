"""Independent exact modular and interpolated boundary-coefficient checks."""
import json
from pathlib import Path
from verify_full_packet_review import P,ROWS,COLS,HEIGHT,entry,determinant_mod,lex_pivots

indices=json.loads(Path(__file__).with_name('full_packet_review_verification.json').read_text())['selected_column_indices']
n=len(ROWS)
b=[[entry(row,COLS[c],(1,1)) for c in indices] for row in ROWS]
e=[]
active=[]
for number,(j,s,beta) in enumerate(ROWS):
    target=(j,s-1,(beta[0]+1,beta[1]))
    source_cost=s+4*beta[0]+10*beta[1]
    if s and source_cost+3>=HEIGHT:
        # The matrices use s!-scaled rows; normalized boundary factor is
        # -(beta_1+1), so scaling the target to the source multiplies by s.
        e.append([-s*(beta[0]+1)*entry(target,COLS[c],(1,1))%P for c in indices])
        active.append(number)
    else:e.append([0]*n)
aug=[x+y for x,y in zip(b,e)]
for k in range(n):
    pivot=next(i for i in range(k,n) if aug[i][k])
    aug[k],aug[pivot]=aug[pivot],aug[k]
    inv=pow(aug[k][k],-1,P)
    aug[k]=[x*inv%P for x in aug[k]]
    for i in range(n):
        if i==k:continue
        f=aug[i][k]
        if f:aug[i]=[(x-f*y)%P for x,y in zip(aug[i],aug[k])]
x=[row[n:] for row in aug]
tr=sum(x[i][i] for i in range(n))%P
trsquare=sum(x[i][j]*x[j][i] for i in range(n) for j in range(n))%P
detb=determinant_mod(b)
coefficient=detb*(tr*tr-trsquare)*pow(2,-1,P)%P
result={'status':'passed','prime':P,'coefficient_eta2_zeta1_2':coefficient,'original_period_scaled_determinant':detb,'trace':tr,'trace_square':trsquare,'active_first_order_rows':len(active),'scope':'Aggregate coefficient for the corrected determinant before restoring b^D, and before its source row factorial normalization. Nonzero modular value certifies a nonzero formal-period polynomial; actual-period nonvanishing follows by transcendence.'}
assert coefficient

# A second computation extracts the degree-two coefficient by exact polynomial
# interpolation, without using the trace identity. E1 has nonzero rows only
# at `active`, so Sylvester reduces every determinant to 44 by 44.
r=len(active)
aug=[row+[int(i==j) for j in active] for i,row in enumerate(b)]
for k in range(n):
    pivot=next(i for i in range(k,n) if aug[i][k])
    aug[k],aug[pivot]=aug[pivot],aug[k]
    inv=pow(aug[k][k],-1,P)
    aug[k]=[v*inv%P for v in aug[k]]
    for i in range(n):
        if i==k:continue
        f=aug[i][k]
        if f:aug[i]=[(u-f*v)%P for u,v in zip(aug[i],aug[k])]
q=[row[n:] for row in aug]
kmat=[[sum(e[active[i]][ell]*q[ell][j] for ell in range(n))%P for j in range(r)] for i in range(r)]
values=[]
for eta in range(r+1):
    values.append(determinant_mod([[(int(i==j)+eta*kmat[i][j])%P for j in range(r)] for i in range(r)]))
# Check the Sylvester reduction against the unreduced matrix at two points.
for eta in (2,7):
    direct=determinant_mod([[(b[i][j]+eta*e[i][j])%P for j in range(n)] for i in range(n)])
    assert direct==detb*values[eta]%P
forward=[];current=values[:]
while current:
    forward.append(current[0]);current=[(current[j+1]-current[j])%P for j in range(len(current)-1)]
interpolated=0;harmonic=0
for degree in range(2,r+1):
    harmonic=(harmonic+pow(degree-1,-1,P))%P
    binomial_coefficient=pow(-1,degree,P)*harmonic*pow(degree,-1,P)%P
    interpolated=(interpolated+forward[degree]*binomial_coefficient)%P
assert detb*interpolated%P==coefficient
rank_e1=len(lex_pivots(e))
assert rank_e1==44
witness=next((i,j,(sum(x[i][k]*x[k][j] for k in range(n))+x[i][j])%P) for i in range(n) for j in range(n) if (sum(x[i][k]*x[k][j] for k in range(n))+x[i][j])%P)
result.update({'rank_first_order_boundary':rank_e1,'nonzero_entry_of_X_squared_plus_X':list(witness)})
result.update({'independent_interpolation_nodes':r+1,'interpolated_normalized_eta2_coefficient':interpolated,'full_determinant_cross_checks':[2,7],'determinant_values_after_sylvester_reduction':values})
Path(__file__).with_name('aggregate_coefficient_verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='determinant_values_after_sylvester_reduction'},indent=2))
