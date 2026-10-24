"""Synthetic regression tests for genotype-aware candidate and CDS summaries.

Uses temporary generated records only; no patient data or network. Requires
pysam and Biopython from the documented analysis environment.
"""
import ast,json,os,subprocess,sys,tempfile,unittest
from pathlib import Path
import pysam
R=Path(__file__).resolve().parents[1]

def helper(script):
    tree=ast.parse((R/'scripts'/script).read_text())
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='called_alt_indices')
    scope={};exec(compile(ast.Module(body=[node],type_ignores=[]),script,'exec'),scope)
    return scope['called_alt_indices']

def fixture(root, multiallelic=True):
    for p in ['data/references/ucsc','data/references/gnomad_review_genes','work','results/candidates','results/phenotype']:(root/p).mkdir(parents=True)
    seq='ATG'+'GCT'*1049+'TAA'
    fa=root/'data/references/ucsc/GRCh38.primary.no_chr.fa';fa.write_text('>15\n'+seq+'\n');pysam.faidx(str(fa))
    gene={'start':1,'end':3153,'canonical_transcript':'TX1.1','Transcript':[{'id':'TX1','strand':1,'Translation':{'start':1,'end':3153},'Exon':[{'start':1,'end':3153}]}]}
    (root/'data/references/bub1b_ensembl.json').write_text(json.dumps(gene))
    (root/'data/references/bub1b_gnomad.json').write_text(json.dumps({'data':{'gene':{'variants':[]}}}))
    (root/'results/candidates/review_gene_manifest.json').write_text(json.dumps({'genes':['SYNTHETIC']}))
    (root/'results/phenotype/semantic_ranking.json').write_text(json.dumps({'rankings':{label:[{'gene':'SYNTHETIC','rank':1}] for label in ['all_eight','proband_seven','without_cancer']}}))
    h=pysam.VariantHeader();h.contigs.add('15',length=len(seq));h.add_sample('synthetic')
    for key,number,typ in [('GT',1,'String'),('DP',1,'Integer'),('GQ',1,'Integer')]:h.formats.add(key,number,typ,key)
    for key in ['BCSQ','CLNSIG','CLNREVSTAT','ORIG','GENEINFO','CLNDN']:h.info.add(key,'.','String',key)
    cases=[(4,('A',),(0,0),False),(7,('A',),(None,None),False),(10,('A',),(None,1),False),(13,('A',),(0,1),False),(16,('A',),(1,1),False),(25,('A',),(1,0),True)]
    if multiallelic:cases += [(19,('A','C'),(1,2),False),(22,('A','C'),(0,2),False),(28,('*',),(0,1),False)]
    plain=root/'work/synthetic.vcf'
    with pysam.VariantFile(str(plain),'w',header=h) as out:
        for pos,alts,gt,phased in sorted(cases):
            v=out.new_record(contig='15',start=pos-1,alleles=('G',)+alts);v.filter.add('PASS')
            v.samples['synthetic']['GT']=gt;v.samples['synthetic']['DP']=35;v.samples['synthetic']['GQ']=50;v.samples['synthetic'].phased=phased
            v.info['BCSQ']=('missense|SYNTHETIC|TX1|protein_coding|+|synthetic',);v.info['CLNSIG']=('Likely_pathogenic',)
            out.write(v)
    for name in ['primary.normalized.vcf.gz','primary.csq.clinvar.vcf.gz']:
        target=root/'work'/name;pysam.tabix_compress(str(plain),str(target),force=True);pysam.tabix_index(str(target),preset='vcf',force=True)

def run(script, root):
    return subprocess.run([sys.executable,'-B',str(R/'scripts'/script)],env={**os.environ,'MVA_PROJECT_ROOT':str(root),'PYTHONDONTWRITEBYTECODE':'1'},capture_output=True,text=True)

class GenotypeCarrierTests(unittest.TestCase):
    def test_complete_called_allele_membership(self):
        for script in ['summarize_candidates.py','reconstruct_bub1b.py']:
            f=helper(script)
            for gt,n,expected in [(None,1,()),((None,None),1,()),((None,1),1,()),((0,0),1,()),((0,1),1,(1,)),((1,0),1,(1,)),((1,1),1,(1,)),((1,2),2,(1,2)),((0,2),2,(2,))]:
                with self.subTest(script=script,gt=gt):self.assertEqual(f(gt,n),expected)
            with self.assertRaises(ValueError):f((0,2),1)
    def test_reconstruction_only_carried_snv_alleles(self):
        with tempfile.TemporaryDirectory(prefix='mva-gt-test-') as td:
            p=Path(td);fixture(p);r=run('reconstruct_bub1b.py',p);self.assertEqual(r.returncode,0,r.stderr)
            rows=json.loads((p/'results/candidates/bub1b_cds_reconstruction.json').read_text())['variants']
            self.assertEqual([(x['position'],x['alt']) for x in rows],[(13,'A'),(16,'A'),(19,'A'),(19,'C'),(22,'C'),(25,'A')])
    def test_summary_excludes_reference_and_incomplete_genotypes(self):
        with tempfile.TemporaryDirectory(prefix='mva-gt-test-') as td:
            p=Path(td);fixture(p,False);r=run('summarize_candidates.py',p);self.assertEqual(r.returncode,0,r.stderr)
            d=json.loads((p/'results/candidates/annotation_summary.json').read_text())
            self.assertEqual(d['rare_or_unknown_rows'],3);self.assertEqual(d['counts']['noncarrier_or_incomplete_gt_records_excluded'],3)
            c=json.loads((p/'results/candidates/clinvar_primary_chromosome_findings.json').read_text());self.assertEqual(len(c),3)
    def test_summary_refuses_unsplit_multiallelic_consequence_assignment(self):
        with tempfile.TemporaryDirectory(prefix='mva-gt-test-') as td:
            p=Path(td);fixture(p);r=run('summarize_candidates.py',p)
            self.assertNotEqual(r.returncode,0);self.assertIn('split biallelic',r.stderr)
            self.assertFalse((p/'results/candidates/review_genes_rare_or_unknown.tsv').exists())

if __name__=='__main__':unittest.main()
