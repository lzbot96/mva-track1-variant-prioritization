"""Synthetic checks for mate identity, two-SNV assumptions and bounded inventory."""
import copy,hashlib,importlib.util,json,subprocess,sys,tempfile,unittest
from pathlib import Path
SCRIPTS=Path(__file__).resolve().parents[1]/'scripts'
spec=importlib.util.spec_from_file_location('phase_recheck',SCRIPTS/'recheck_phase_reads.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class PhaseRecheckTests(unittest.TestCase):
    def row(self):
        return dict(chrom_1='chr1',pos_1='101',ref_1='A',alt_1='G',chrom_2='chr1',pos_2='301',ref_2='C',alt_2='T')

    def test_casava_mates_share_pair_identity(self):
        self.assertEqual(m.normalize_pair_id('@synthetic:lane:tile:x:y 1:N:0:index'),m.normalize_pair_id('@synthetic:lane:tile:x:y 2:N:0:index'))

    def test_legacy_suffix_and_unrelated_suffix(self):
        self.assertEqual(m.normalize_pair_id('@synthetic/1'),m.normalize_pair_id('@synthetic/2'))
        self.assertNotEqual(m.normalize_pair_id('@synthetic/10'),m.normalize_pair_id('@synthetic/1'))
        with self.assertRaises(ValueError):m.normalize_pair_id('not-fastq')
        with self.assertRaises(ValueError):m.normalize_pair_id('@/1')

    def test_reject_indel_missing_pair_and_different_chromosome(self):
        self.assertEqual(len(m.validate_variants([self.row()])),2)
        for edit in [dict(alt_1='AG'),dict(ref_2=''),dict(chrom_2='chr2'),dict(pos_2='101'),dict(alt_1='A'),dict(pos_1='1')]:
            r=self.row();r.update(edit)
            with self.subTest(edit=edit),self.assertRaises(ValueError):m.validate_variants([r])
        with self.assertRaises(ValueError):m.validate_variants([])

    def test_anchor_format_rejects_invalid_length_alphabet_duplicate(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'anchors.tsv'
            for text in ['x '+('A'*30), 'x '+('N'*31), 'x '+('A'*31)+'\nx '+('C'*31)]:
                p.write_text(text)
                with self.subTest(text=text),self.assertRaises(ValueError):m.read_patterns(p)
            p.write_text('x '+('A'*31)+'\n');self.assertEqual(m.read_patterns(p)['x'],'A'*31)

    def test_candidate_anchor_center_and_reverse_complement(self):
        variants=m.validate_variants([self.row()]);refs={'candidate_1':'T'*15+'A'+'C'*15,'candidate_2':'A'*15+'C'+'G'*15}
        class Ref:
            def fetch(self,chrom,start,end):return refs['candidate_1' if start==85 else 'candidate_2']
        patterns={};rc=lambda x:x.translate(str.maketrans('ACGT','TGCA'))[::-1]
        for v in variants:
            ref=refs[v['label']];alt=ref[:15]+v['alt']+ref[16:]
            for suffix,seq in [('REF_F',ref),('ALT_F',alt),('REF_R',rc(ref)),('ALT_R',rc(alt))]:patterns[v['label']+'_'+suffix]=seq
        m.validate_candidate_patterns(patterns,variants,Ref())
        wrong=copy.copy(patterns);wrong['candidate_1_ALT_F']=wrong['candidate_1_REF_F']
        with self.assertRaises(ValueError):m.validate_candidate_patterns(wrong,variants,Ref())

    def test_asset_flags_are_measured_not_assumed(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'work').mkdir();ref=root/'ref.fa'
            empty=m.bounded_assets(root,ref);self.assertEqual(empty['alignment_file_count'],0);self.assertFalse(empty['classic_bwa_reference_index_complete'])
            (root/'work'/'synthetic.bam').write_bytes(b'')
            for ext in ['.amb','.ann','.bwt','.pac','.sa']:Path(str(ref)+ext).write_bytes(b'')
            observed=m.bounded_assets(root,ref);self.assertEqual(observed['alignment_file_count'],1);self.assertTrue(observed['classic_bwa_reference_index_complete'])
            self.assertEqual(observed['parental_or_long_read_data'],'Not inferred by this bounded file inventory')

    def test_source_hash_detects_changed_bytes_and_preserves_receipts(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);source=root/'synthetic.vcf.gz';source.write_bytes(b'synthetic');h=hashlib.sha256(source.read_bytes()).hexdigest()
            manifest=root/'manifest.json';receipt=root/'receipt.json';out=root/'result.json'
            manifest.write_text(json.dumps({'files':[dict(path=source.name,size=9,sha256=h)]}))
            receipt.write_text(json.dumps({'files':[dict(path=source.name,bytes=9,sha256=h,upstream_sha256_verified=False)]}))
            cmd=[sys.executable,str(SCRIPTS/'check_source_checksum.py'),'--source-vcf',str(source),'--manifest',str(manifest),'--receipt',str(receipt),'--output',str(out)]
            self.assertEqual(subprocess.run(cmd,capture_output=True).returncode,0)
            saved=out.read_bytes();self.assertNotEqual(subprocess.run(cmd,capture_output=True).returncode,0);self.assertEqual(out.read_bytes(),saved)
            source.write_bytes(b'changed');cmd[-1]=str(root/'changed.json')
            self.assertNotEqual(subprocess.run(cmd,capture_output=True).returncode,0)
            self.assertFalse(json.loads((root/'changed.json').read_text())['matches_authenticated_manifest'])

if __name__=='__main__':unittest.main()
