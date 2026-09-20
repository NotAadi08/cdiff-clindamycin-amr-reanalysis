from pathlib import Path
import argparse, json, joblib, pandas as pd
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description='Predict clindamycin resistance from presence/absence of packaged AMR features.')
parser.add_argument('--features',required=True,help='Comma-separated present feature names, e.g. erm(B),cfr(C)')
a=parser.parse_args()
bundle=joblib.load(ROOT/'models/cdiff_clindamycin_model.joblib')
present={x.strip() for x in a.features.split(',') if x.strip()}
X=pd.DataFrame([[int(f in present) for f in bundle['features']]],columns=bundle['features'])
p=float(bundle['model'].predict_proba(X)[:,1][0]); pred='resistant' if p>=bundle['threshold'] else 'susceptible'
print(json.dumps({'present_features':sorted(present),'model_features':bundle['features'],'resistance_probability':p,'prediction':pred},indent=2))
