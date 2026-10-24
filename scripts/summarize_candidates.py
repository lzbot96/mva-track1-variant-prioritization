"""Summarize local consequences and public frequencies; missing AF stays unknown."""
from pathlib import Path
import csv,json,re,collections
import pysam
from project_paths import project_root, required_env
def called_alt_indices(genotype, alternate_count):
 """Carried ALT indices in a complete called genotype; phase is irrelevant.

 Missing or partial genotypes are excluded rather than treated as confident
 carriers. Check indices against the actual ALT list, including multiallelic GT.
 """
 if not genotype or any(a is None for a in genotype):return ()
 if any(not isinstance(a,int) or a<0 or a>alternate_count for a in genotype):
  raise ValueError('Genotype allele index is outside the record alleles')
 return tuple(i for i in range(1,alternate_count+1) if i in genotype)

R = project_root()
review=set(json.loads((R/'results/candidates/review_gene_manifest.json').read_text())['genes'])
ranks={k:{r['gene']:r['rank'] for r in v} for k,v in json.loads((R/'results/phenotype/semantic_ranking.json').read_text())['rankings'].items()}
freq={};loaded=set()
paths=[R/'data/references/bub1b_gnomad.json']+list((R/'data/references/gnomad_review_genes').glob('*.json'))
for p in paths:
 try:
  d=json.loads(p.read_text());gene=d.get('data',{}).get('gene')
  if not gene:continue
  loaded.add(p.stem)
  for v in gene.get('variants',[]):
   vals=[v[x]['af'] for x in ('exome','genome') if v.get(x) and v[x].get('af') is not None]
   freq[v['variant_id']]={'max_overall_af':max(vals) if vals else None,'exome_ac':(v.get('exome') or {}).get('ac'),'exome_an':(v.get('exome') or {}).get('an'),'genome_ac':(v.get('genome') or {}).get('ac'),'genome_an':(v.get('genome') or {}).get('an')}
 except (ValueError,AttributeError):continue
fields=['chrom','pos','ref','alt','gene','consequence','transcript','amino_acid_change','gt','dp','gq','ad','filter','baseline_qc_pass','max_overall_af','frequency_status','clinvar','clinvar_review','all_eight_rank','proband_seven_rank','without_cancer_rank','original_record']
rows=[];clinical=[];count=collections.Counter()
vcf=pysam.VariantFile(R/'work/primary.csq.clinvar.vcf.gz')
for v in vcf:
 count['records']+=1
 if len(v.samples)!=1:raise ValueError('Candidate summary requires one sample')
 sample=next(iter(v.samples.values()));gt=sample.get('GT');dp=sample.get('DP');gq=sample.get('GQ')
 carried=called_alt_indices(gt,len(v.alts or ()))
 if not carried:
  count['noncarrier_or_incomplete_gt_records_excluded']+=1
  continue
 if len(v.alts or ())!=1:
  raise ValueError('Candidate summary requires split biallelic records; normalize with bcftools norm -m -any before annotation')
 count['called_alt_records']+=1
 base={'chrom':v.contig,'pos':v.pos,'ref':v.ref,'alt':','.join(v.alts),'gt':('|' if sample.phased else '/').join(str(a) for a in gt),'dp':dp,'gq':gq,'ad':','.join(str(a) for a in (sample.get('AD') or [])),'filter':';'.join(v.filter),'baseline_qc_pass':1 in carried and list(v.filter)==['PASS'] and dp is not None and gq is not None and dp>=10 and gq>=20,'original_record':v.info.get('ORIG','')}
 cl=v.info.get('CLNSIG',());cl=','.join(cl) if isinstance(cl,tuple) else str(cl)
 cr=v.info.get('CLNREVSTAT',());cr=','.join(cr) if isinstance(cr,tuple) else str(cr)
 base.update({'clinvar':cl,'clinvar_review':cr})
 pathogenic=not 'Conflicting' in cl and any(s in ('Pathogenic','Likely_pathogenic') for s in re.split(r'[/|]',cl))
 if pathogenic:
  extra={'geneinfo':v.info.get('GENEINFO'),'disease':v.info.get('CLNDN'),'consequences':v.info.get('BCSQ')}
  clinical.append({**base,**extra})
 seen=set()
 for csq in v.info.get('BCSQ',()):
  parts=csq.split('|')
  if len(parts)<4:continue
  consequence,gene,transcript=parts[:3]
  if gene not in review:continue
  if not any(effect in consequence for effect in ['missense','stop_gained','frameshift','splice_acceptor','splice_donor','start_lost','stop_lost','inframe','in_frame']):continue
  key=(gene,consequence,transcript)
  if key in seen:continue
  seen.add(key)
  fid=f'{v.contig}-{v.pos}-{v.ref}-{v.alts[0]}';f=freq.get(fid)
  row={**base,'gene':gene,'consequence':consequence,'transcript':transcript,'amino_acid_change':parts[5] if len(parts)>5 else '', 'max_overall_af':f['max_overall_af'] if f else None, 'frequency_status':'observed_public_frequency' if f and f['max_overall_af'] is not None else 'unknown_or_not_in_downloaded_gene_summary'}
  row.update({label+'_rank':rank.get(gene,'') for label,rank in ranks.items()})
  rows.append(row)
count['clinvar_pathogenic_or_likely_pathogenic']=len(clinical)
count['review_gene_consequence_rows']=len(rows)
with (R/'results/candidates/review_genes_coding.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
rare=[r for r in rows if r['baseline_qc_pass'] and (r['max_overall_af'] is None or r['max_overall_af']<=0.01)]
with (R/'results/candidates/review_genes_rare_or_unknown.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rare)
(R/'results/candidates/clinvar_primary_chromosome_findings.json').write_text(json.dumps(clinical,indent=2)+'\n')
(R/'results/candidates/annotation_summary.json').write_text(json.dumps({'counts':dict(count),'frequency_gene_files_loaded':len(loaded),'review_genes':len(review),'rare_or_unknown_rows':len(rare),'limitations':'Only complete called genotypes carrying ALT are summarized; missing/partial GT and reference-only records are excluded. Input must be split biallelic before annotation. Frequency unknown is not rare; overall AF not ancestry-specific; multiple rows can describe one variant across transcripts. Primary nuclear chromosomes only. Site-local consequences do not establish phase.'},indent=2)+'\n')
print(json.dumps(dict(count)));print('Rare-or-unknown candidate rows:',len(rare))
