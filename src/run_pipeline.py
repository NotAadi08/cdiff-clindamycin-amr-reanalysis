from pathlib import Path
import argparse, json, hashlib
import pandas as pd
import numpy as np
import joblib
import matplotlib.pyplot as plt
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.metrics import (roc_auc_score, average_precision_score, accuracy_score,
    balanced_accuracy_score, f1_score, recall_score, confusion_matrix,
    RocCurveDisplay, PrecisionRecallDisplay)

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description='Reproduce the C. difficile clindamycin AMR analysis.')
parser.add_argument('--data', type=Path, default=ROOT/'data/raw/phenotype_genotype_merged.csv.gz', help='Path to phenotype_genotype_merged.csv.gz')
args=parser.parse_args()
DATA=args.data
if not DATA.exists():
    raise SystemExit(f'Raw input not found: {DATA}\nSee data/README.md for download, placement, and checksum instructions.')
PROCESSED=ROOT/'data/processed'; RESULTS=ROOT/'results'; FIGURES=ROOT/'figures'; MODELS=ROOT/'models'
for d in [PROCESSED,RESULTS,FIGURES,MODELS]: d.mkdir(parents=True,exist_ok=True)

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

df=pd.read_csv(DATA,low_memory=False)
sub=df[(df['organism']=='Clostridioides difficile') &
       (df['antibiotic_name']=='clindamycin') &
       df['resistance_phenotype'].isin(['resistant','susceptible'])].copy()
lab=sub[['BioSample_ID','resistance_phenotype']].drop_duplicates().set_index('BioSample_ID')
sub['feature']=sub['amr_element_symbol'].fillna(sub['gene_symbol']).fillna(sub['reference_name'])
sub=sub[sub['feature'].notna()]
X=pd.crosstab(sub['BioSample_ID'],sub['feature']).clip(upper=1).reindex(lab.index,fill_value=0)
X=X.loc[:,X.sum(axis=0)>=2]
y=(lab.loc[X.index,'resistance_phenotype']=='resistant').astype(int)
X.to_csv(PROCESSED/'feature_matrix.csv')
pd.DataFrame({'BioSample_ID':X.index,'phenotype':lab.loc[X.index,'resistance_phenotype'].values,'label':y.values}).to_csv(PROCESSED/'labels.csv',index=False)
(PROCESSED/'feature_schema.json').write_text(json.dumps(list(X.columns),indent=2))

models={
 'Logistic regression':Pipeline([('scale',StandardScaler()),('model',LogisticRegression(max_iter=5000,class_weight='balanced',C=0.5))]),
 'Random forest':RandomForestClassifier(n_estimators=250,class_weight='balanced',min_samples_leaf=2,max_features='sqrt',random_state=42),
 'Extra Trees':ExtraTreesClassifier(n_estimators=250,class_weight='balanced',min_samples_leaf=2,max_features='sqrt',random_state=42)
}
cv=RepeatedStratifiedKFold(n_splits=5,n_repeats=5,random_state=42)
rows=[]
for name,m in models.items():
    vals={k:[] for k in ['roc_auc','pr_auc','accuracy','balanced_accuracy','f1','sensitivity','specificity']}
    for tr,te in cv.split(X,y):
        m.fit(X.iloc[tr],y.iloc[tr]); p=m.predict_proba(X.iloc[te])[:,1]; pred=(p>=.5).astype(int); yt=y.iloc[te]
        vals['roc_auc'].append(roc_auc_score(yt,p)); vals['pr_auc'].append(average_precision_score(yt,p))
        vals['accuracy'].append(accuracy_score(yt,pred)); vals['balanced_accuracy'].append(balanced_accuracy_score(yt,pred))
        vals['f1'].append(f1_score(yt,pred,zero_division=0)); vals['sensitivity'].append(recall_score(yt,pred,zero_division=0))
        vals['specificity'].append(recall_score(yt,pred,pos_label=0,zero_division=0))
    for metric,v in vals.items(): rows.append([name,metric,np.mean(v),np.std(v),np.percentile(v,2.5),np.percentile(v,97.5)])
res=pd.DataFrame(rows,columns=['model','metric','mean','sd','p2.5','p97.5'])
res.to_csv(RESULTS/'cv_metrics.csv',index=False)
best=res[res.metric=='roc_auc'].sort_values('mean',ascending=False).iloc[0]['model']; bm=models[best]
skf=StratifiedKFold(n_splits=5,shuffle=True,random_state=2026)
p=cross_val_predict(bm,X,y,cv=skf,method='predict_proba')[:,1]; pred=(p>=.5).astype(int)
pd.DataFrame({'BioSample_ID':X.index,'label':y.values,'probability':p,'prediction':pred}).to_csv(RESULTS/'oof_predictions.csv',index=False)
cm=confusion_matrix(y,pred); np.savetxt(RESULTS/'confusion_matrix.csv',cm,delimiter=',',fmt='%d')
bm.fit(X,y)
joblib.dump({'model':bm,'features':list(X.columns),'target':'clindamycin resistance','organism':'Clostridioides difficile','threshold':0.5},MODELS/'cdiff_clindamycin_model.joblib')
if best=='Logistic regression': imp=np.abs(bm.named_steps['model'].coef_[0])
else: imp=bm.feature_importances_
fi=pd.DataFrame({'feature':X.columns,'importance':imp}).sort_values('importance',ascending=False); fi.to_csv(RESULTS/'feature_importance.csv',index=False)
prev=[]
for f in fi.head(15).feature:
    for label,val in [('Resistant',1),('Susceptible',0)]:
        z=X.loc[y==val,f]; prev.append([f,label,int(z.sum()),len(z),float(z.mean())])
pd.DataFrame(prev,columns=['feature','phenotype','present','n','prevalence']).to_csv(RESULTS/'top_feature_prevalence.csv',index=False)
fig,ax=plt.subplots(figsize=(6,5)); RocCurveDisplay.from_predictions(y,p,ax=ax,name=best); ax.plot([0,1],[0,1],'--'); ax.set_title('Out-of-fold ROC curve'); fig.tight_layout(); fig.savefig(FIGURES/'roc_curve.png',dpi=220); plt.close(fig)
fig,ax=plt.subplots(figsize=(6,5)); PrecisionRecallDisplay.from_predictions(y,p,ax=ax,name=best); ax.set_title('Out-of-fold precision-recall curve'); fig.tight_layout(); fig.savefig(FIGURES/'pr_curve.png',dpi=220); plt.close(fig)
top=fi.head(12).sort_values('importance'); fig,ax=plt.subplots(figsize=(7,5)); ax.barh(top.feature,top.importance); ax.set_xlabel('Feature importance'); ax.set_title('Top genomic predictors'); fig.tight_layout(); fig.savefig(FIGURES/'feature_importance.png',dpi=220); plt.close(fig)
summary={'data_sha256':sha256(DATA),'n_isolates':int(len(y)),'resistant':int(y.sum()),'susceptible':int((1-y).sum()),'n_features':int(X.shape[1]),'features':list(X.columns),'best_model':best,'oof_roc_auc':float(roc_auc_score(y,p)),'oof_pr_auc':float(average_precision_score(y,p)),'oof_accuracy':float(accuracy_score(y,pred)),'oof_balanced_accuracy':float(balanced_accuracy_score(y,pred)),'oof_f1':float(f1_score(y,pred)),'oof_sensitivity':float(recall_score(y,pred)),'oof_specificity':float(recall_score(y,pred,pos_label=0)),'confusion_matrix':cm.tolist()}
(RESULTS/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
