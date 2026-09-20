from pathlib import Path
import argparse, json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, accuracy_score, balanced_accuracy_score, f1_score, recall_score

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description='Run marker-baseline, source-bias, and false-negative follow-up analyses.')
parser.add_argument('--data', type=Path, default=ROOT/'data/raw/phenotype_genotype_merged.csv.gz', help='Path to phenotype_genotype_merged.csv.gz')
args=parser.parse_args()
if not args.data.exists():
    raise SystemExit(f'Raw input not found: {args.data}\nSee data/README.md for download, placement, and checksum instructions.')
X=pd.read_csv(ROOT/'data/processed/feature_matrix.csv',index_col=0)
labels=pd.read_csv(ROOT/'data/processed/labels.csv').set_index('BioSample_ID')
y=labels.loc[X.index,'label'].astype(int)
raw=pd.read_csv(args.data,low_memory=False)
sub=raw[(raw['organism']=='Clostridioides difficile') & (raw['antibiotic_name']=='clindamycin') & raw['resistance_phenotype'].isin(['resistant','susceptible'])]
meta=sub[['BioSample_ID','database','country','ast_standard','measurement','measurement_sign','measurement_units','resistance_phenotype','assembly_ID','AMR_associated_publications']].drop_duplicates('BioSample_ID').set_index('BioSample_ID').reindex(X.index)
oof=pd.read_csv(ROOT/'results/oof_predictions.csv').set_index('BioSample_ID').reindex(X.index)

rules={
    'erm(B) only':(X['erm(B)']>0).astype(int),
    'cfr(C) only':(X['cfr(C)']>0).astype(int),
    'erm(B) OR cfr(C)':((X['erm(B)']>0)|(X['cfr(C)']>0)).astype(int),
    'Logistic regression OOF':oof['prediction'].astype(int),
}
rows=[]
for name,pred in rules.items():
    cm=confusion_matrix(y,pred,labels=[0,1])
    tn,fp,fn,tp=cm.ravel()
    rows.append([name,accuracy_score(y,pred),balanced_accuracy_score(y,pred),f1_score(y,pred,zero_division=0),recall_score(y,pred,zero_division=0),recall_score(y,pred,pos_label=0,zero_division=0),tn,fp,fn,tp])
pd.DataFrame(rows,columns=['method','accuracy','balanced_accuracy','f1','sensitivity','specificity','TN','FP','FN','TP']).to_csv(ROOT/'results/baseline_comparison.csv',index=False)

source_rows=[]
for db in sorted(meta['database'].dropna().unique()):
    ids=meta.index[meta['database']==db]
    yy=y.loc[ids]
    pred=rules['erm(B) only'].loc[ids]
    cm=confusion_matrix(yy,pred,labels=[0,1]); tn,fp,fn,tp=cm.ravel()
    source_rows.append([db,len(ids),int(yy.sum()),int((1-yy).sum()),accuracy_score(yy,pred),recall_score(yy,pred,zero_division=0),recall_score(yy,pred,pos_label=0,zero_division=0) if (yy==0).any() else np.nan,tn,fp,fn,tp])
pd.DataFrame(source_rows,columns=['database','n','resistant','susceptible','ermB_accuracy','ermB_sensitivity','ermB_specificity','TN','FP','FN','TP']).to_csv(ROOT/'results/source_stratified_analysis.csv',index=False)

fn_ids=oof.index[(oof['label']==1)&(oof['prediction']==0)]
fn=meta.loc[fn_ids].copy()
fn['cfr(C)']=X.loc[fn_ids,'cfr(C)']; fn['cplR']=X.loc[fn_ids,'cplR']; fn['erm(B)']=X.loc[fn_ids,'erm(B)']
fn.reset_index().to_csv(ROOT/'results/false_negative_audit.csv',index=False)

# Original-source linkage and feature nesting
ct=pd.crosstab(X['erm(B)'],X['cfr(C)'])
summary={
 'ermB_and_ml_predictions_identical': bool((rules['erm(B) only'].values==rules['Logistic regression OOF'].values).all()),
 'cfrC_positive_without_ermB': int(((X['cfr(C)']>0)&(X['erm(B)']==0)).sum()),
 'false_negative_ids': list(fn_ids),
 'source_counts': [{'database':a,'phenotype':b,'n':int(n)} for (a,b),n in meta.groupby(['database','resistance_phenotype']).size().items()],
 'interpretation':'At threshold 0.5, logistic-regression OOF classifications are identical to the erm(B)-only rule. cfr(C) is nested within erm(B) in this cohort and does not rescue any erm(B)-negative resistant isolate.'
}
(ROOT/'results/followup_summary.json').write_text(json.dumps(summary,indent=2,default=str))

# Figure
comp=pd.read_csv(ROOT/'results/baseline_comparison.csv')
fig,ax=plt.subplots(figsize=(7,4.5))
x=np.arange(len(comp)); w=0.25
ax.bar(x-w,comp['accuracy'],w,label='Accuracy')
ax.bar(x,comp['sensitivity'],w,label='Sensitivity')
ax.bar(x+w,comp['specificity'],w,label='Specificity')
ax.set_ylim(0,1.05); ax.set_ylabel('Score'); ax.set_xticks(x); ax.set_xticklabels(comp['method'],rotation=18,ha='right'); ax.legend(); ax.set_title('Marker rules versus machine-learning classification')
fig.tight_layout(); fig.savefig(ROOT/'figures/baseline_comparison.png',dpi=220); plt.close(fig)
print(json.dumps(summary,indent=2,default=str))
