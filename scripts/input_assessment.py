"""Local input QC summary, with no patient variants emitted to stdout."""
from pathlib import Path
import collections
import json
import re
import subprocess

from project_paths import project_root, required_env
ROOT = project_root()
raw=ROOT/required_env('MVA_INPUT_VCF')
header=subprocess.check_output(['bcftools','view','-h',str(raw)],text=True)
summary={'reference_headers':[], 'contigs':{},'filter_definitions':[], 'info_descriptions':{}, 'sample_count':0}
for line in header.splitlines():
    if line.startswith('##reference='):summary['reference_headers'].append(line.split('=',1)[1])
    if line.startswith('##contig='):
        m=re.search(r'ID=([^,>]+).*?length=(\d+)',line)
        if m:summary['contigs'][m[1]]=int(m[2])
    if line.startswith('##FILTER='):summary['filter_definitions'].append(line)
    if line.startswith('##INFO='):
        m=re.search(r'ID=([^,>]+).*?Description="(.*?)"',line)
        if m:summary['info_descriptions'][m[1]]=m[2]
    if line.startswith('#CHROM'):summary['sample_count']=len(line.split('\t'))-9
format_string='%FILTER\t[%GT\t%DP\t%GQ\t%AD\t%PGT\t%PID]\n'
proc=subprocess.Popen(['bcftools','query','-f',format_string,str(raw)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
counts=collections.Counter();filters=collections.Counter();depth=collections.Counter();gqhist=collections.Counter();ab=collections.Counter()
for line in proc.stdout:
    fields=line.rstrip('\n').split('\t')
    filt,gt,dp,gq,ad,pgt,pid=fields
    counts['records']+=1;filters[filt]+=1
    called='.' not in gt
    alleles=re.split(r'[/|]',gt)
    nonref=called and any(int(a)>0 for a in alleles)
    counts['called_genotypes']+=called
    counts['nonreference_genotypes']+=nonref
    counts['phased_separator']+='|' in gt
    counts['physical_phase_fields_present']+=pgt!='.' and pid!='.'
    if dp!='.':depth[min(int(dp),500)]+=1
    else:counts['dp_missing']+=1
    if gq!='.':gqhist[min(int(float(gq)),100)]+=1
    else:counts['gq_missing']+=1
    baseline=filt=='PASS' and dp!='.' and gq!='.' and int(dp)>=10 and float(gq)>=20 and nonref
    counts['baseline_pass_dp10_gq20_nonref']+=baseline
    if called and len(alleles)==2 and set(alleles)=={'0','1'} and ad!='.':
        vals=ad.split(',')
        if len(vals)==2 and '.' not in vals and sum(map(int,vals))>0:
            r,a=map(int,vals);ab[min(int(100*a/(r+a)),100)]+=1
            counts['biallelic_het_with_allele_depth']+=1
            counts['biallelic_het_ad_sum_differs_dp']+=dp!='.' and r+a!=int(dp)
err=proc.stderr.read();code=proc.wait()
if code:raise RuntimeError('bcftools query failed; inspect local log')
summary.update({'counts':dict(counts),'filters':dict(filters),'dp_histogram_capped_500':dict(depth),'gq_histogram_capped_100':dict(gqhist),'biallelic_het_alt_fraction_percent_histogram':dict(ab),'warnings':err})
out=ROOT/'results/qc/input_assessment.json';out.write_text(json.dumps(summary,indent=2)+'\n')
print('Local VCF assessment completed; results saved without emitting patient content.')
