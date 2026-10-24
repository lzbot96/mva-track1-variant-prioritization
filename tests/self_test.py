"""No patient input, network, package install, or analysis execution."""
import ast, csv, gzip, json, os, subprocess, sys, tempfile
from pathlib import Path
R=Path(__file__).resolve().parents[1]
for p in (R/'scripts').glob('*.py'): ast.parse(p.read_text(), filename=p.name)
for p in (R/'scripts').glob('*.sh'): subprocess.run(['bash','-n',str(p)],check=True)
with tempfile.TemporaryDirectory(prefix='mva-code-test-') as td:
    t=Path(td);exe=t/'scanner'
    subprocess.run(['c++','-std=c++17','-O2',str(R/'scripts/scan_allele_kmers.cpp'),'-lz','-o',str(exe)],check=True)
    pattern='ACGTTGCACTGATCGATGCTAGCTAACGTGC'
    assert len(pattern)==31
    reverse=pattern.translate(str.maketrans('ACGT','TGCA'))[::-1]
    (t/'patterns.tsv').write_text('synthetic_F\t'+pattern+'\nsynthetic_R\t'+reverse+'\n')
    seqs=[('NN'+pattern+'NN','I'),(pattern,'5'),(pattern,'!'),(pattern+'N'+pattern,'I'),(reverse,'I'),('N'*40,'I')]
    with gzip.open(t/'reads.gz','wt') as f:
        for i,(seq,q) in enumerate(seqs):f.write(f'@synthetic_{i}\n{seq}\n+\n{q*len(seq)}\n')
    args=[str(exe),str(t/'patterns.tsv'),str(t/'reads.gz'),str(t/'counts.tsv'),str(t/'matched.fastq')]
    subprocess.run(args,check=True,capture_output=True)
    with (t/'counts.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
    assert [int(x['reads']) for x in rows]==[4,1]
    assert [int(x['q20_reads']) for x in rows]==[3,1]
    assert [int(x['q30_reads']) for x in rows]==[2,1]
    assert all(int(x['file_records'])==6 for x in rows)
    assert len((t/'matched.fastq').read_text().splitlines())==20
    (t/'patterns.tsv').write_text('bad\tACGT\n')
    assert subprocess.run(args,capture_output=True).returncode==3
    (t/'patterns.tsv').write_text('synthetic\t'+pattern+'\n')
    with gzip.open(t/'reads.gz','wt') as f:f.write('@bad\n'+pattern+'\n+\n')
    assert subprocess.run(args,capture_output=True).returncode==5
    env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','MVA_PROJECT_ROOT':td}
    pre=subprocess.run([sys.executable,str(R/'scripts/preflight.py'),'--root',td,'--stage','scan'],env=env,capture_output=True,text=True)
    assert pre.returncode==2 and 'FASTQ files' in json.loads(pre.stdout)['stages']['scan']['missing']
    # Ensure missing root configuration fails before writing any output.
    env.pop('MVA_PROJECT_ROOT')
    missing=subprocess.run([sys.executable,str(R/'scripts/run_kmer_validation.py')],env=env,capture_output=True,text=True)
    assert missing.returncode!=0 and 'Set MVA_PROJECT_ROOT' in missing.stderr
print('PASS: Python AST, shell syntax, scanner counts/orientations/ambiguous bases/quality/repeats, invalid pattern, truncated FASTQ, preflight and missing-root guard.')
