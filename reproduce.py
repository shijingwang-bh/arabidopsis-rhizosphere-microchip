"""Recompute analyses and/or redraw the current manuscript figures in a new directory."""
from pathlib import Path
import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time

PACKAGE=Path(__file__).resolve().parent
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['all','analysis','figures'],default='all')
    parser.add_argument('--output',type=Path,required=True,help='New or empty directory outside this package')
    parser.add_argument('--rscript',default=os.environ.get('RSCRIPT',shutil.which('Rscript') or 'Rscript'))
    parser.add_argument('--r-library',type=Path,help='Optional existing R library directory')
    args=parser.parse_args()
    out=args.output.resolve()
    if out==PACKAGE or out.is_relative_to(PACKAGE) or PACKAGE.is_relative_to(out):
        parser.error('--output must be outside the package')
    if out.exists() and any(out.iterdir()):parser.error('--output must be new or empty')
    pydeps=set()
    rdeps={'jsonlite'}
    if args.mode in ['analysis','all']:
        pydeps.add('openpyxl')
        rdeps.update(['car','dplyr','vegan','edgeR','limma','statmod','fgsea','BiocParallel'])
    if args.mode in ['figures','all']:
        pydeps.update(['PIL','pypdf','reportlab'])
        rdeps.update(['ggplot2','patchwork','magick','ggrepel','circlize'])
    missing=[p for p in sorted(pydeps) if importlib.util.find_spec(p) is None]
    if missing:parser.error('Missing Python packages: '+', '.join(missing))
    env=os.environ.copy();env['PYTHONUTF8']='1'
    if args.r_library:env['R_LIBS_USER']=str(args.r_library.resolve(strict=True))
    expression='p<-commandArgs(TRUE);m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];if(length(m))stop(paste("Missing R packages:",paste(m,collapse=", ")))'
    check=subprocess.run([args.rscript,'--vanilla','-e',expression,*sorted(rdeps)],env=env)
    if check.returncode:return check.returncode
    checked=subprocess.run([sys.executable,str(PACKAGE/'validate_package.py')],env=env)
    if checked.returncode:return checked.returncode
    out.mkdir(parents=True,exist_ok=True)
    shutil.copytree(PACKAGE/'reference/plot_data',out/'plot_data')
    (out/'logs').mkdir()
    steps=[]
    if args.mode in ['all','analysis']:
        steps.extend(['analyse_phenotypes.R','extract_OD_workbooks.py','analyse_OD.R',
                      'analyse_RNA.R','analyse_LCMS_global.R','analyse_LCMS_interactions.R',
                      'verify_reproduction.py','sync_results.py'])
    if args.mode in ['all','figures']:
        steps.extend(['plot_platform_panels.R','plot_current_figures.R','compose_figures.py'])
    report=[]
    for i,name in enumerate(steps,1):
        script=PACKAGE/'code'/name
        command=([args.rscript,'--vanilla'] if script.suffix=='.R' else [sys.executable])+[str(script),str(PACKAGE),str(out)]
        log=out/'logs'/f'{i:02d}_{script.stem}.log'
        print('RUN',name,flush=True);start=time.time()
        with log.open('w',encoding='utf-8') as f:
            completed=subprocess.run(command,cwd=PACKAGE,env=env,stdout=f,stderr=subprocess.STDOUT)
        report.append({'step':name,'exit_code':completed.returncode,'seconds':round(time.time()-start,2),'log':str(log.relative_to(out))})
        (out/'run_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
        if completed.returncode:
            print('FAILED:',name,'; see',log,file=sys.stderr)
            return completed.returncode
    print('PASS:',len(report),'steps;',args.mode,'outputs:',out)
    return 0

if __name__=='__main__':raise SystemExit(main())
