import os
import sys

import motif_find

CMP_DIR = 'cmp'
OUT_DIR = 'out'

MOTIF_FILE = 'my_motif.tsv'
SEQ_FILE = 'ARG81.fa'
PROBS = ('0.01', '0.1', '0.99')


def _test_all(capsys, alg, p, q):
    for seq in open(SEQ_FILE, 'r'):
        seq = seq.strip()
        sys.argv = [sys.argv[0], '--alg', alg, seq, MOTIF_FILE, p, q]
        motif_find.main()
        captured = capsys.readouterr()
        cmp_name = f'{seq}_p{p}_q{q}.{alg}'
        # open(os.path.join(CMP_DIR, cmp_name), 'w').write(captured.out)
        cmp = open(os.path.join(CMP_DIR, cmp_name)).read()
        if alg == "forward" or alg == "backward":
            assert round(float(captured.out.strip()), 2) == round(float(cmp.strip()), 2)
        else:
            assert captured.out.strip() == cmp.strip()



def _test_alg(capsys, alg):
    for p in PROBS:
        for q in PROBS:
            _test_all(capsys, alg, p, q)


def test_forward(capsys):
    _test_alg(capsys, 'forward')


def test_backward(capsys):
    _test_alg(capsys, 'backward')


def test_posterior(capsys):
    _test_alg(capsys, 'posterior')


def test_viterbi(capsys):
    _test_alg(capsys, 'viterbi')
