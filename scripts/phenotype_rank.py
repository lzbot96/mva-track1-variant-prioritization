"""Transparent phenotype-only baseline using gene-frequency information content.
Not Exomiser; query-to-gene best common ancestor semantic similarity.
"""
from pathlib import Path
import csv,json,math
from collections import defaultdict,Counter
from functools import lru_cache
from project_paths import project_root, required_env
R = project_root()
parents={};names={};term=None
for line in (R/'data/references/hpo/hp.obo').read_text().splitlines():
 if line=='[Term]':term=None
 elif line.startswith('id: HP:'):term=line[4:];parents.setdefault(term,set())
 elif term and line.startswith('name: '):names[term]=line[6:]
 elif term and line.startswith('is_a: HP:'):parents[term].add(line.split()[1])
@lru_cache(None)
def ancestors(term):
 out={term}
 for p in parents.get(term,[]):out.update(ancestors(p))
 return frozenset(out)
genes=defaultdict(set)
with (R/'data/references/hpo/genes_to_phenotype.txt').open() as f:
 for row in csv.DictReader(f,delimiter='\t'):
  if row['gene_symbol'] not in ('-', '') and row['ncbi_gene_id'].isdigit():
   genes[row['gene_symbol']].add(row['hpo_id'])
prop={g:set().union(*(ancestors(t) for t in ts)) for g,ts in genes.items()}
freq=Counter(t for ts in prop.values() for t in ts)
ic={t:-math.log(n/len(genes)) for t,n in freq.items()}
phen=json.loads((R/'results/phenotype/curated_hpo.json').read_text())
queries={'all_eight':phen['all_challenge_hpo'],'proband_seven':phen['proband_only_hpo'],'without_cancer':[t for t in phen['proband_only_hpo'] if t!=required_env('MVA_EXCLUDED_HPO')]}
allresults={}
for label,terms in queries.items():
 rows=[]
 for gene,ts in prop.items():
  scores=[max((ic.get(a,0) for a in ancestors(t)&ts),default=0) for t in terms]
  row={'gene':gene,'mean_resnik':sum(scores)/len(scores),'exact_query_terms':sum(t in genes[gene] for t in terms),'query_terms':len(terms),'per_term':dict(zip(terms,scores))};rows.append(row)
 rows.sort(key=lambda r:(-r['mean_resnik'],-r['exact_query_terms'],r['gene']))
 for rank,row in enumerate(rows,1):row['rank']=rank
 allresults[label]=rows
 with (R/f'results/phenotype/ranking_{label}.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=['rank','gene','mean_resnik','exact_query_terms','query_terms'],delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
(R/'results/phenotype/semantic_ranking.json').write_text(json.dumps({'method':'Query-to-gene mean Resnik similarity; IC from ancestor-propagated unique gene counts; no genotype input','genes':len(genes),'query_names':{t:names.get(t) for t in queries['all_eight']},'rankings':allresults},indent=2)+'\n')
print('Phenotype ranking complete; controlled outputs saved locally.')
