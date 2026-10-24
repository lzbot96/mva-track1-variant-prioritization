#!/usr/bin/env python3
"""Recheck a two-SNV candidate pair using existing private calls and selected reads.

Requires pysam and BWA. All outputs, including selected sequences, local reference,
SAM alignments and case-specific diagnostics, are PRIVATE controlled derivatives.
This is not a full-read rescan, independent call set, genome-wide mapping check,
causal calibration or automated cis/trans inference. The original data are read-only.
"""
import argparse,collections,csv,datetime,hashlib,json,re,subprocess
from pathlib import Path


def normalize_pair_id(header):
    """Normalize CASAVA headers and legacy /1,/2 suffixes to a pair identifier."""
    fields=header.strip().split()
    if not fields or not fields[0].startswith('@'):
        raise ValueError('FASTQ header must start with @')
    token=fields[0][1:]
    token=re.sub(r'/[12]$','',token)
    if not token:
        raise ValueError('FASTQ pair identifier is empty')
    return token


def validate_variants(rows):
    """Require one candidate row containing two distinct, same-chromosome SNVs."""
    if len(rows)!=1:
        raise ValueError('Exactly one two-SNV candidate row is required')
    variants=[]
    for slot in ('1','2'):
        row=rows[0]
        chrom=row.get('chrom_'+slot,'').removeprefix('chr')
        try: pos=int(row.get('pos_'+slot,''))
        except ValueError: raise ValueError('Candidate position must be a positive integer') from None
        ref=row.get('ref_'+slot,'');alt=row.get('alt_'+slot,'')
        if not chrom or pos<16 or len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT' or ref==alt:
            raise ValueError('Both candidates must be distinct-base SNVs with complete 31-base reference windows')
        variants.append(dict(label='candidate_'+slot,chrom=chrom,pos=pos,ref=ref,alt=alt))
    if variants[0]['chrom']!=variants[1]['chrom'] or variants[0]['pos']==variants[1]['pos']:
        raise ValueError('Require two distinct positions on the same chromosome')
    return variants


def read_patterns(path):
    patterns={}
    for line in path.read_text().splitlines():
        fields=line.split()
        if len(fields)!=2:
            raise ValueError('Each anchor line must have a name and sequence')
        name,seq=fields
        if name in patterns or len(seq)!=31 or not set(seq)<=set('ACGT'):
            raise ValueError('Anchor names must be unique and sequences must be 31 uppercase A/C/G/T bases')
        patterns[name]=seq
    if not patterns:raise ValueError('No anchors supplied')
    return patterns


def validate_candidate_patterns(patterns,variants,reference):
    if len(patterns)!=8:
        raise ValueError('Expect REF/ALT forward/reverse anchors for exactly two SNVs')
    expected_names=set()
    rc=lambda seq:seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
    for v in variants:
        ref=reference.fetch(v['chrom'],v['pos']-16,v['pos']+15)
        if len(ref)!=31 or ref[15]!=v['ref']:
            raise ValueError('Reference window or central candidate base disagrees')
        names=[k for k,s in patterns.items() if k.endswith('_REF_F') and s==ref]
        if len(names)!=1:raise ValueError('A candidate has no unique matching forward reference anchor')
        prefix=names[0][:-6]
        alt=ref[:15]+v['alt']+ref[16:]
        expected={prefix+'_REF_F':ref,prefix+'_ALT_F':alt,prefix+'_REF_R':rc(ref),prefix+'_ALT_R':rc(alt)}
        if any(patterns.get(k)!=s for k,s in expected.items()):
            raise ValueError('An anchor does not match the candidate SNV or reverse complement')
        expected_names.update(expected)
    if expected_names!=set(patterns):raise ValueError('Unexpected or ambiguous candidate anchor groups')


def bounded_assets(project,reference_path):
    """Inventory named raw/work folders; do not infer family or sequencing platform."""
    folders=[project/'data/raw',project/'work']
    alignment_paths=[]
    for folder in folders:
        if folder.exists():
            alignment_paths.extend(str(x.relative_to(project)) for x in folder.rglob('*')
                                   if x.is_file() and x.suffix.lower() in {'.bam','.cram','.sam'})
    parts=[Path(str(reference_path)+suffix) for suffix in ['.amb','.ann','.bwt','.pac','.sa']]
    return {'inventory_scope':['data/raw','work'],
            'alignment_files_found':sorted(alignment_paths),
            'alignment_file_count':len(alignment_paths),
            'classic_bwa_reference_index_complete':all(x.is_file() for x in parts),
            'classic_bwa_reference_index_parts_found':sum(x.is_file() for x in parts),
            'parental_or_long_read_data':'Not inferred by this bounded file inventory'}


def main():
    import pysam
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project',type=Path,required=True)
    p.add_argument('--csv',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--source-vcf',type=Path)
    p.add_argument('--source-index',type=Path)
    p.add_argument('--bwa',default='bwa',help='BWA executable, resolved through PATH by default')
    a=p.parse_args()
    if a.out.exists() and any(a.out.iterdir()):
     raise ValueError('Output directory must be new or empty; prior receipts are never overwritten')
    a.out.mkdir(parents=True,exist_ok=True)
    def sha(x):return hashlib.sha256(x.read_bytes()).hexdigest()
    def write(name,d):(a.out/name).write_text(json.dumps(d,indent=2)+'\n')
    rows=list(csv.DictReader(a.csv.open()))
    variants=validate_variants(rows)
    fa_path=a.project/'data/references/ucsc/GRCh38.primary.no_chr.fa';fa=pysam.FastaFile(str(fa_path))
    asset_inventory=bounded_assets(a.project,fa_path)
    source_vcf=a.source_vcf
    if source_vcf is None:
     candidates=list((a.project/'data/raw').glob('*.vcf.gz'))
     if len(candidates)!=1:raise ValueError('Pass --source-vcf when the raw directory does not contain exactly one VCF')
     source_vcf=candidates[0]
    vcfs={'supplied':source_vcf,'normalized':a.project/'work/primary.normalized.vcf.gz','annotated':a.project/'work/primary.csq.clinvar.vcf.gz'}
    vcf_checks={};source_hashes={str(a.csv):sha(a.csv)}
    for label,path in vcfs.items():
     vf=pysam.VariantFile(str(path), index_filename=str(a.source_index) if label=='supplied' and a.source_index else None);
     if len(vf.header.samples)!=1:raise ValueError('This check requires a single-sample VCF')
     checks=[]
     for v in variants:
      chrom=v['chrom'] if v['chrom'] in vf.header.contigs else 'chr'+v['chrom']
      hits=[r for r in vf.fetch(chrom,v['pos']-1,v['pos']) if r.pos==v['pos'] and r.ref==v['ref'] and v['alt'] in (r.alts or ())]
      
      if len(hits)!=1:raise ValueError('A candidate allele is not uniquely found in the VCF')
      r=hits[0];s=next(iter(r.samples.values()));ai=r.alts.index(v['alt'])+1;gt=s.get('GT')
      checks.append({'candidate':v['label'],'gt':gt,'called_alt_index':ai,'alt_is_present_in_gt':ai in (gt or ()),'heterozygous_ref_alt':sorted(gt)==[0,ai] if gt and None not in gt else False,'phased':s.phased,'phase_fields':{k:s.get(k) for k in ('PS','PGT','PID','HP') if k in vf.header.formats},'filter':list(r.filter),'dp':s.get('DP'),'gq':s.get('GQ'),'ad':s.get('AD'),'reference_agrees':fa.fetch(v['chrom'],v['pos']-1,v['pos']-1+len(v['ref']))==v['ref']})
     
     if label=='supplied':
      source_hashes['supplied_vcf_current_bytes']=sha(path)
      if a.source_index:source_hashes['private_source_index']=sha(a.source_index)
     vcf_checks[label]={'sample_count':len(vf.header.samples),'phase_header_fields':[k for k in ('PS','PGT','PID','HP') if k in vf.header.formats],'checks':checks};vf.close()
    patterns=read_patterns(a.project/'work/allele_kmers.tsv')
    validate_candidate_patterns(patterns,variants,fa)
    source_hashes['work/allele_kmers.tsv']=sha(a.project/'work/allele_kmers.tsv')
    old=json.loads((a.project/'results/read_validation/summary.json').read_text());source_hashes['results/read_validation/summary.json']=sha(a.project/'results/read_validation/summary.json')
    counts={k:collections.Counter(reads=0,q20_reads=0,q30_reads=0) for k in patterns};ids=collections.defaultdict(set);read_positions=collections.defaultdict(list);read_lengths=collections.Counter();selected=[];per_file={};malformed=[]
    for file in sorted((a.project/'results/read_validation').glob('*.matched.fastq')):
     source_hashes[str(file.relative_to(a.project))]=sha(file);n=0;lines=file.read_text().splitlines();
     file_counts={k:collections.Counter(reads=0,q20_reads=0,q30_reads=0) for k in patterns}
     if len(lines)%4:raise ValueError('Truncated saved FASTQ record')
     for ix in range(0,len(lines),4):
      h,seq,plus,q=lines[ix:ix+4]
      if not h.startswith('@') or not plus.startswith('+') or not seq or len(seq)!=len(q):raise ValueError('Malformed saved FASTQ record')
      pair_id=normalize_pair_id(h)
      n+=1;read_lengths[len(seq)]+=1;found=[]
      for name,anchor in patterns.items():
       pos=seq.find(anchor)
       if pos<0:continue
       found.append(name);qual=ord(q[pos+15])-33
       counts[name]['reads']+=1;counts[name]['q20_reads']+=qual>=20;counts[name]['q30_reads']+=qual>=30
       file_counts[name]['reads']+=1;file_counts[name]['q20_reads']+=qual>=20;file_counts[name]['q30_reads']+=qual>=30
       ids[name.rsplit('_',1)[0]].add(pair_id);read_positions[name].append(pos+16)
      
      if not found:raise ValueError('Saved selected read has no current 31-base anchor')
      selected.append({'anon_name':f'read_{len(selected)+1:04d}','sequence':seq,'qualities':q,'patterns':found,'pair_identifier_hash':hashlib.sha256(pair_id.encode()).hexdigest()})
     per_file[file.name]={'read_records':n,'counts_by_pattern':{k:dict(v) for k,v in file_counts.items()}}
    all_count_match={k:dict(counts[k])==old['counts_by_pattern'][k] for k in counts}
    allele_groups=sorted(ids);between=[]
    for i,l in enumerate(allele_groups):
     for rr in allele_groups[i+1:]:
      if l.rsplit('_',1)[0]!=rr.rsplit('_',1)[0]:between.append({'allele_1':l,'allele_2':rr,'shared_saved_read_pair_ids':len(ids[l]&ids[rr])})
    read_check={'selected_read_records':len(selected),'counts_by_file':per_file,'read_lengths':dict(read_lengths),'counts_by_pattern':{k:dict(v) for k,v in counts.items()},'all_saved_pattern_counts_match_prior_summary':all(all_count_match.values()),'pattern_count_comparisons':all_count_match,'unique_read_pair_identifiers_by_allele':{k:len(v) for k,v in ids.items()},'cross_site_pair_linkage':between,'same_pair_ids_across_sites':sum(v['shared_saved_read_pair_ids'] for v in between),'scan_provenance':'Fresh re-read of existing matched FASTQ subset and recomputation; original full-FASTQ scan not rerun.','ascertainment_bias':'Reads were selected by exact candidate anchors; counts do not estimate unbiased depth or allele fraction.','completion_receipt_validation':'Only fields and legacy shape inspected; no fullscan receipt is regenerated'}
    complete=json.loads((a.project/'results/read_validation/complete.json').read_text());read_check['original_completion_receipt_fields']=sorted(complete);read_check['original_completion_receipt_is_legacy']=sorted(complete)==['completed']
    # A deliberately regional alignment of saved selected reads; asset availability is inventoried below.
    chroms={v['chrom'] for v in variants}
    start=max(0,min(v['pos'] for v in variants)-10001);end=max(v['pos'] for v in variants)+10000
    ref=a.out/'private_regional_reference.fa';region=fa.fetch(variants[0]['chrom'],start,end);ref.write_text('>candidate_interval\n'+region+'\n')
    fq=a.out/'private_selected_reads.fastq';fq.write_text(''.join('@'+x['anon_name']+'\n'+x['sequence']+'\n+\n'+x['qualities']+'\n' for x in selected))
    cmds=[]
    with (a.out/'alignment_commands.log').open('w') as log:
     cmd=[a.bwa,'index',str(ref)];cmds.append(cmd);subprocess.run(cmd,check=True,stdout=log,stderr=log)
     cmd=[a.bwa,'mem','-t','2',str(ref),str(fq)];cmds.append(cmd)
     with (a.out/'private_regional_alignments.sam').open('w') as sam:subprocess.run(cmd,check=True,stdout=sam,stderr=log)
    alignments=pysam.AlignmentFile(str(a.out/'private_regional_alignments.sam'),'r');primary=[];candidate_pileups={v['label']:collections.Counter() for v in variants};cigars=collections.Counter();nms=[];localmapq=collections.Counter();softclips=0
    for rec in alignments:
     if rec.is_secondary or rec.is_supplementary:continue
     if rec.is_unmapped:continue
     primary.append(rec.query_name);cigars[rec.cigarstring]+=1;nms.append(rec.get_tag('NM') if rec.has_tag('NM') else None);localmapq[rec.mapping_quality]+=1
     softclips+=any(op==4 for op,length in (rec.cigartuples or []))
     pairmap={rp:qp for qp,rp in rec.get_aligned_pairs() if rp is not None}
     for v in variants:
      rel=v['pos']-1-start
      if rel not in pairmap:continue
      qp=pairmap[rel]
      if qp is None:candidate_pileups[v['label']]['deletion']+=1
      else:
       base=rec.query_sequence[qp]
       category='ref' if base==v['ref'] else ('alt' if base==v['alt'] else 'other')
       candidate_pileups[v['label']][category]+=1
       if rec.query_qualities[qp]>=30:candidate_pileups[v['label']][category+'_q30']+=1
    alignment={'method':'Fresh BWA-MEM single-end regional alignment of saved exact-anchor-selected reads to a local GRCh38 interval.','selected_reads':len(selected),'primary_mapped_reads':len(primary),'reads_with_soft_clipping':softclips,'edit_distance_distribution':dict(collections.Counter(nms)),'cigar_distribution':dict(cigars),'regional_mapq_distribution':dict(localmapq),'candidate_aligned_base_counts':{k:dict(v) for k,v in candidate_pileups.items()},'limits':['Not a genome-wide alignment; regional MAPQ is not evidence of genome-wide uniqueness.','Input is candidate-anchor-selected; not independent variant calling or unbiased depth/allele balance.','Single-end alignment of this subset does not establish read-pair fragment orientation, molecular deduplication or phase.'],'commands':cmds}
    assets=asset_inventory
    assets.update(source_sample_count=vcf_checks['supplied']['sample_count'],saved_subset_read_lengths=sorted(read_lengths),candidate_distance_bases=abs(variants[1]['pos']-variants[0]['pos']))
    result={'performed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'new_computations':['Indexed recheck of supplied, normalized and annotated candidate genotypes/phase fields and GRCh38 REF agreement','Fresh reparse/count/base-quality check of saved matched reads and shared-pair-identifier test','Fresh regional full-read alignment of the saved selected subset'],'vcf_checks':vcf_checks,'saved_read_recheck':read_check,'regional_alignment':alignment,'assets':assets,'conclusions':{'submitted_alts_are_called_heterozygous':all(c['heterozygous_ref_alt'] and c['alt_is_present_in_gt'] for x in vcf_checks.values() for c in x['checks']),'cis_trans_inference_performed':False,'phase_assessment_scope':'Phase fields and shared pair IDs are reported, not converted into a cis/trans inference','independent_variant_calling_completed':False,'genome_wide_coverage_completed':False,'genome_wide_alignment_completed':False,'legacy_fullscan_receipt_regenerated':False},'source_hashes':source_hashes,'variants_private':variants}
    write('PRIVATE_PHASE_READ_RECHECK.json',result)
    summary={'performed_at_utc':result['performed_at_utc'],'genotypes_called_heterozygous':result['conclusions']['submitted_alts_are_called_heterozygous'],'saved_reads':len(selected),'saved_counts_match_previous':read_check['all_saved_pattern_counts_match_prior_summary'],'cross_site_shared_pair_ids':read_check['same_pair_ids_across_sites'],'regional_primary_mapped_reads':len(primary),'regional_soft_clipped_reads':softclips,'regional_candidate_counts':alignment['candidate_aligned_base_counts'],'phase_established_by_this_workflow':False,'full_genome_alignment_or_independent_calling':False,'old_fullscan_receipt_regenerated':False}
    write('SHAREABLE_SUMMARY.json',summary);print(json.dumps(summary))

if __name__=="__main__":
    main()
