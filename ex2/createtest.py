import numpy as np
from random import choice
def create_test( ):
    def _create_test( ):
        tau1, tau2 = { 'A' : 0.4 , 'G' : 0, 'C' : 0.6 , 'T' : 0   }, { 'A' : 1 , 'G' : 0, 'C' : 0 , 'T' : 0   }

        seq = ""
        m = 300
        for _ in range(m):
            seq +=choice(['A', 'G', 'C', 'T'])

        p,q = 0.5, 0.5

        _sum = 0
        for i,w in enumerate(seq[1:-2]): 
            xi = tau1[w] * tau2[seq[i+1]]
            _sum += xi * 0.25**(m-2) * (1-p)**(m-2)*(p)**2 *q
        _sum += (1-q)*(0.25**m) * (1-p)**(m-1)*p

        print("python3 ./motif_find.py --alg forward {0} ./othertest.tsv {1} {2} ".format(seq, p, q))
        return str(np.log(_sum)) 
    _l = [ ]
    for _ in range(3):
        _l.append(_create_test())
    
    print("\n".join(_l))

if __name__ == "__main__":
	create_test()