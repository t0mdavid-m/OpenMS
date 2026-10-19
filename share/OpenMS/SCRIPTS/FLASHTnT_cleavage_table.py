"""Learn a FLASHTnT cleavage table (-ex:cleavage_table) from one run.

For every bond of an annotated, unmodified target proteoform: is it covered by a terminal fragment of the series
(b/y-like: a b x y; c/z-like: c z zp1 zp2)? The table holds log2(rate / mean) per residue before (P1) and after (P1') the bond.

  FLASHTnT -in run_deconv.mzML -fasta db.fasta -out_prsm prsm.tsv ...           (first pass, built-in tables)
  FLASHAnnotator -in run.mzML -in_prsm prsm.tsv -out annot.tsv -out_annotation env.tsv -ion_type <the same ion types>
  python3 FLASHTnT_cleavage_table.py run_cleavage.tsv annot.tsv env.tsv [annot2.tsv env2.tsv ...]
  FLASHTnT ... -ex:cleavage_table run_cleavage.tsv                              (second pass)

A series with fewer than --min_prsms PrSMs is not written and keeps its built-in values.
--pseudo bonds at the mean rate are added to every residue, which shrinks rare residues towards no preference.
"""
import argparse, collections, csv, math

SERIES = {'a': 'by', 'b': 'by', 'x': 'by', 'y': 'by', 'c': 'cz', 'z': 'cz', 'zp1': 'cz', 'zp2': 'cz'}
AA = 'ACDEFGHIKLMNPQRSTVWY'


def count(pairs):
    total = {g: [collections.Counter(), collections.Counter()] for g in ('by', 'cz')}
    cleaved = {g: [collections.Counter(), collections.Counter()] for g in ('by', 'cz')}
    prsms = collections.Counter()
    for annot, env in pairs:
        seqs = {}
        for r in csv.DictReader(open(annot), delimiter='\t'):
            if 'DECOY' in r['ProteinAccession'] or r['ModCount'] not in ('0', ''): continue
            seqs[r['PrSMIndex']] = r['ProteinSequence']
        bonds = collections.defaultdict(lambda: collections.defaultdict(set))
        for e in csv.DictReader(open(env), delimiter='\t'):
            if e['PrSMIndex'] not in seqs or float(e['DecoyShift']) != 0.0 or e['IonType'].startswith('int'): continue
            g = SERIES.get(e['IonType'])
            if g: bonds[e['PrSMIndex']][g].add(int(e['Bond']))
        for k, s in seqs.items():
            for g, hit in bonds[k].items():
                if not hit: continue
                prsms[g] += 1
                for b in range(1, len(s)):
                    total[g][0][s[b - 1]] += 1
                    total[g][1][s[b]] += 1
                    if b in hit:
                        cleaved[g][0][s[b - 1]] += 1
                        cleaved[g][1][s[b]] += 1
    return total, cleaved, prsms


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('out', help='cleavage table to write (TSV)')
    ap.add_argument('files', nargs='+', help='pairs of FLASHAnnotator -out and -out_annotation files')
    ap.add_argument('--min_prsms', type=int, default=30)
    ap.add_argument('--pseudo', type=float, default=200)
    a = ap.parse_args()
    if len(a.files) % 2: ap.error('give FLASHAnnotator -out and -out_annotation files in pairs')
    total, cleaved, prsms = count(list(zip(a.files[::2], a.files[1::2])))
    with open(a.out, 'w') as f:
        f.write('# series\tside\tresidue\tlog2(rate/mean)\n')
        for g in ('by', 'cz'):
            n, c = sum(total[g][0].values()), sum(cleaved[g][0].values())
            if n == 0 or c == 0: continue
            keep = prsms[g] >= a.min_prsms
            f.write(f'# {g}: {prsms[g]} PrSMs, {n} bonds, {c} cleaved' + ('' if keep else ' - too few PrSMs, built-in values kept') + '\n')
            if not keep: continue
            mean = c / n
            for side, i in (('before', 0), ('after', 1)):
                for r in AA:
                    rate = (cleaved[g][i][r] + a.pseudo * mean) / (total[g][i][r] + a.pseudo)
                    f.write(f'{g}\t{side}\t{r}\t{math.log2(rate / mean):+.3f}\n')


if __name__ == '__main__':
    main()
