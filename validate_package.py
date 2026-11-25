"""Validate release checksums, current figure coverage and complete BH families."""
from pathlib import Path
import csv
import hashlib
import math
import sys

ROOT=Path(__file__).resolve().parent
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def bh(p):
    order=sorted(range(len(p)),key=p.__getitem__);q=[0.]*len(p);running=1.
    for rank in range(len(p),0,-1):
        i=order[rank-1];running=min(running,p[i]*len(p)/rank);q[i]=running
    return q
def check_bh(label,table,pcol,qcol):
    p=[float(r[pcol]) for r in table];q=[float(r[qcol]) for r in table]
    assert all(math.isfinite(x) and 0<=x<=1 for x in p+q),label
    err=max(abs(a-b) for a,b in zip(q,bh(p)))
    assert err<1e-10,(label,err)
    return len(table)
def main():
    expected={}
    for line in (ROOT/'checksums_sha256.txt').read_text(encoding='utf-8').splitlines():
        sha,relative=line.split('  ',1);path=(ROOT/relative).resolve()
        assert path.is_relative_to(ROOT),relative
        assert path.is_file(),relative+' missing'
        assert digest(path)==sha,relative+' checksum mismatch'
        expected[relative]=sha
    actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and
            '__pycache__' not in p.parts and '.git' not in p.relative_to(ROOT).parts and p.name!='checksums_sha256.txt'}
    assert actual==set(expected),('Unlisted or missing files',actual^set(expected))
    assert max((ROOT/p).stat().st_size for p in expected)<25*2**20
    figures=['Figure_'+str(i) for i in range(1,6)]+['Figure_S'+s for s in ['1','2','3','4','5a','5b','6','7']]
    assert all((ROOT/'reference/figures'/(f+'.pdf')).is_file() for f in figures)
    ref=ROOT/'reference/plot_data';tested=0;families=0
    for label,file,pcol,qcol in [
        ('Pre-inoculation four tests',ref/'pre_inoculation/Figure1d_Welch_effects.csv','p_value','fdr'),
        ('Bacterial two days',ref/'bacteria/Figure3c_daywise_Welch_effects.csv','p_value','fdr'),
        ('Growth four traits',ref/'growth/Growth_traits_average_inoculation_effects.csv','p_value','q_value'),
        ('All RNA genes',ROOT/'reference/statistics/rna/all_gene_statistics.csv','P.Value','adj.P.Val')]:
        table=rows(file);tested+=check_bh(label,table,pcol,qcol);families+=1
    for filename in ['CAMERA_GO_BP_primary_all.csv','CAMERA_GO_BP_estimated_correlation_sensitivity.csv','FRY_GO_BP_correlated_self_contained_all.csv']:
        tested+=check_bh(filename,rows(ROOT/'reference/statistics/rna'/filename),'PValue','FDR');families+=1
    table=rows(ref/'chemistry_global/feature_contrasts_primary_named_QC30.csv')
    for contrast in sorted({r['contrast'] for r in table}):
        part=[r for r in table if r['contrast']==contrast]
        assert len(part)==3126,contrast
        tested+=check_bh(contrast,part,'p','q');families+=1
    table=rows(ref/'chemistry_global/global_factorial_effects_family_sensitivity.csv')
    for group in sorted({r['family'] for r in table}):
        part=[r for r in table if r['family']==group];assert len(part)==7
        tested+=check_bh(group,part,'p','q');families+=1
    print(f'PASS: {len(expected)} files, 13 current figures, {families} complete BH families, {tested} tests checked')
    return 0
if __name__=='__main__':raise SystemExit(main())
