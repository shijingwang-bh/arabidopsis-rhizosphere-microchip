"""Compare recomputed complete test families with frozen manuscript references."""
from pathlib import Path
import csv
import json
import math
import sys

pkg=Path(sys.argv[1]).resolve();run=Path(sys.argv[2]).resolve();report=[]
def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def compare(label,new,old,keys,values,tolerance=1e-8):
    a=rows(new);b=rows(old)
    index={tuple(r[k] for k in keys):r for r in a}
    assert len(index)==len(a),label+' duplicate keys'
    assert len(a)==len(b),(label,len(a),len(b))
    maximum=0.
    for r in b:
        s=index[tuple(r[k] for k in keys)]
        for v in values:
            x,y=r[v],s[v]
            if x in ('','NA','NaN') and y in ('','NA','NaN'):continue
            d=abs(float(x)-float(y));assert math.isfinite(d),(label,v,x,y)
            maximum=max(maximum,d)
    assert maximum<tolerance,(label,maximum)
    report.append({'analysis':label,'rows_compared':len(a),'maximum_absolute_error':maximum,'tolerance':tolerance})
ref=pkg/'reference/plot_data';result=run/'results'
for stem,group,keys,values in [
 ('Figure1d_Welch_effects','pre_inoculation',['metric'],['mean_difference','ci95_low','ci95_high','p_value','fdr']),
 ('Figure3c_daywise_Welch_effects','bacteria',['day'],['mean_difference','ci95_low','ci95_high','p_value','fdr']),
 ('Growth_traits_TypeII','growth',['metric','term'],['f_value','p_value']),
 ('Growth_traits_HC3','growth',['metric','term'],['f_value','p_value']),
 ('Growth_traits_average_inoculation_effects','growth',['metric'],['estimate','ci95_low','ci95_high','p_value','q_value'])]:
    compare(stem,result/'phenotypes'/(stem+'.csv'),ref/group/(stem+'.csv'),keys,values)
compare('RNA genes',result/'rna/all_gene_statistics.csv',pkg/'reference/statistics/rna/all_gene_statistics.csv',
        ['gene_id'],['logFC','P.Value','adj.P.Val'])
compare('RNA GO sets',result/'rna/all_GO_statistics.csv',pkg/'reference/statistics/rna/all_GO_statistics.csv',
        ['GO_ID'],['FDR','FRY_q','estimated_rho_q','fgsea_NES'])
compare('Chemical global terms',result/'chemistry_global/tables/global_factorial_effects_family_sensitivity.csv',
        ref/'chemistry_global/global_factorial_effects_family_sensitivity.csv',['family','term'],['F','p','q','partial_R2'])
compare('All primary LCMS contrasts',result/'chemistry_global/tables/feature_contrasts_primary_named_QC30.csv',
        ref/'chemistry_global/feature_contrasts_primary_named_QC30.csv',['contrast','feature_id'],['log2_effect','se','df','p','q'])
compare('Primary interaction module',result/'chemistry_interactions/tables/all_primary_feature_contrasts.csv',
        ref/'chemistry_interactions/all_primary_feature_contrasts.csv',['contrast','feature_id'],['log2_effect','se','df','p','q'])
compare('Conditioned-medium curves',result/'conditioned_medium/tables/main_curve_mean_SD_SE.csv',
        ref/'conditioned_medium/main_curve_mean_SD_SE.csv',['collection_day','group','time_index'],['n_wells','mean','sd','se'])
(run/'reproduction_report.json').write_text(json.dumps({'passed':True,'checks':report},indent=2),encoding='utf-8')
print(json.dumps({'passed':True,'complete_families_and_summaries_checked':len(report),
                  'rows_compared':sum(r['rows_compared'] for r in report)},indent=2))
