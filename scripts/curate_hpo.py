"""Historical first-table mapping with reviewer-controlled family-history term.
This is a transcription rule, not automated clinical interpretation.
"""
import csv,json,re
from project_paths import project_root,required_env
R=project_root();family=required_env('MVA_FAMILY_HISTORY_HPO')
x=json.loads((R/'results/phenotype/clinical_document_extracted.json').read_text());rows=[]
for row in x['tables'][0][1:]:
    if len(row)!=4:raise ValueError('Expected four columns in clinical table 1')
    feature,term,hpo,note=row
    if not re.fullmatch(r'HP:\d{7}',hpo):raise ValueError('Invalid HPO identifier in source table')
    rows.append({'hpo_id':hpo,'source_term':term,'feature':feature,'status':'present','subject':'family_history' if hpo==family else 'proband','source':'clinical DOCX table 1','note':note})
if not rows or family not in {r['hpo_id'] for r in rows}:raise ValueError('Reviewer-selected family-history term not found')
(R/'results/phenotype/curated_hpo.json').write_text(json.dumps({'all_challenge_hpo':[r['hpo_id'] for r in rows],'proband_only_hpo':[r['hpo_id'] for r in rows if r['subject']=='proband'],'explicit_absent_hpo':[],'rows':rows,'interpretation':'All supplied table terms retained; sensitivity view excludes the reviewer-selected family-history term. No unreported phenotype is treated as absent.'},indent=2)+'\n')
with (R/'results/phenotype/curated_hpo.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
print('Curated HPO transcription saved locally; review subject assignments before ranking.')
