from pathlib import Path
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter
import csv,json,hashlib,shutil,math,sys
PACKAGE=Path(sys.argv[1]).resolve();OUT=Path(sys.argv[2]).resolve()/'results/conditioned_medium'
for name in ['inputs','data','tables','qa']: (OUT/name).mkdir(parents=True,exist_ok=True)
records=[]; caches=[]; checks={}; hashes={}
groupmap={'1g':'stationary_1g','μg':'clinorotation','ms':'MS_control','lb':'LB_control'}
for day in [3,7]:
    src=PACKAGE/f'data/conditioned_medium/day{day}.xlsx'
    hashes[str(day)]={'original_name':src.name,'sha256':hashlib.sha256(src.read_bytes()).hexdigest()}
    shutil.copy2(src,OUT/'inputs'/f'growth_day{day}.xlsx')
    v=load_workbook(src,data_only=True).worksheets[0]
    f=load_workbook(src,data_only=False).worksheets[0]
    if day==3:
        blocks=[('assay',list(range(2,20)),3),('MS_blank',list(range(22,27)),3),('LB_blank',list(range(31,35)),3)]
        cache_rows={27:'MS_blank',35:'LB_blank'}
    else:
        blocks=[('assay',list(range(37,42))+list(range(45,50))+list(range(53,58))+list(range(61,66)),2),
                ('MS_blank',list(range(68,74)),2),('LB_blank',list(range(78,84)),2)]
        cache_rows={74:'MS_blank',84:'LB_blank'}
    for kind,rows,start in blocks:
        for r in rows:
            label=str(v.cell(r,1).value);well=label.rstrip('*')
            group=groupmap[str(v.cell(r,2).value)] if day==3 and kind=='assay' else ({'1':'stationary_1g','2':'clinorotation','3':'MS_control','4':'LB_control'}[well[-1]] if kind=='assay' else kind)
            if day==7:
                rawrow=next(rr for rr in range(2,34) if str(v.cell(rr,1).value)==well)
                assert all(v.cell(rawrow,c).value==v.cell(r,c).value for c in range(2,291))
            vals=[v.cell(r,c).value for c in range(start,start+289)]
            assert len(vals)==289 and all(isinstance(x,(float,int)) and math.isfinite(x) for x in vals)
            for i,x in enumerate(vals):
                records.append({'collection_day':day,'well':well,'source_label':label,'role':kind,'group':group,
                                'included_main':not(day==7 and well=='E2'),'time_index':i+1,'time_h':i*.25,
                                'od_raw':x,'source_sheet':'Sheet1','source_cell':f'{get_column_letter(start+i)}{r}'})
    for row,group in cache_rows.items():
        start=3 if day==3 else 2
        for i in range(289):
            cell=v.cell(row,start+i);formula=f.cell(row,start+i).value
            assert isinstance(formula,str) and formula.startswith('=AVERAGE(')
            caches.append({'collection_day':day,'kind':'blank_mean','group':group,'well':'','time_index':i+1,
                           'cached_value':cell.value,'source_cell':cell.coordinate,'formula':formula})
    if day==7:
        for row in range(86,105):
            well=v.cell(row,1).value
            for i in range(289):
                cell=v.cell(row,i+2)
                caches.append({'collection_day':day,'kind':'corrected','group':'','well':well,'time_index':i+1,
                               'cached_value':cell.value,'source_cell':cell.coordinate,'formula':f.cell(row,i+2).value})
    checks[f'day{day}_numeric_measurements']=sum(x['collection_day']==day for x in records)
    checks[f'day{day}_source_sha256']=hashes[str(day)]['sha256']
def save_csv(name,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8-sig') as h:
        w=csv.DictWriter(h,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
save_csv('data/growth_measurements_raw.csv',records)
save_csv('data/workbook_cached_calculations.csv',caches)
(OUT/'inputs/source_hashes.json').write_text(json.dumps(hashes,indent=2,ensure_ascii=False),encoding='utf-8')
checks['day7_duplicate_blocks_identical']=True
checks['day7_E2_exclusion']='Uninoculated E2 well; retain four inoculated clinorotation wells in the main summary.'
(OUT/'qa/extraction_checks.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(json.dumps(checks,ensure_ascii=False))
