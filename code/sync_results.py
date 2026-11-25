"""Place recomputed tables into an isolated run's plotting-data tree."""
from pathlib import Path
import shutil
import sys

pkg=Path(sys.argv[1]).resolve();run=Path(sys.argv[2]).resolve()
mapping={
 'Figure1d_Welch_effects':'pre_inoculation','Figure1d_group_summaries':'pre_inoculation',
 'Figure3c_daywise_Welch_effects':'bacteria','Figure3b_fluorescence_summaries':'bacteria',
 'Figure3c_cell_summaries':'bacteria','Growth_traits_TypeII':'growth','Growth_traits_HC3':'growth',
 'Growth_traits_average_inoculation_effects':'growth','Growth_traits_group_summaries':'growth'}
for stem,group in mapping.items():
    source=run/'results/phenotypes'/(stem+'.csv')
    shutil.copy2(source,run/'plot_data'/group/source.name)
for group,folder in [('chemistry_global','chemistry_global'),('chemistry_interactions','chemistry_interactions'),
                      ('conditioned_medium','conditioned_medium')]:
    for source in (run/'results'/folder/'tables').glob('*.csv'):
        target=run/'plot_data'/group/source.name
        if target.exists():shutil.copy2(source,target)
shutil.copy2(run/'plot_data/rna_degs/gene_statistics.csv',run/'plot_data/rna_degs/main_gene_statistics.csv')
print('Synchronized recomputed tables to isolated plotting inputs')
