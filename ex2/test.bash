python3 ./createtest.py > ./simpletest
bash ./simpletest


# python3 ./motif_find.py --alg forward GAAG ./mytest.tsv .9 .5
# python3 ./motif_find.py --alg backward GAAG ./mytest.tsv .9 .5
# python3 ./motif_find.py --alg posterior GAAG ./mytest.tsv .9 .5
# python3 ./motif_find.py --alg posterior GAAGAAG ./mytest.tsv .9 .5
# python3 ./motif_find.py --alg viterbi GAAGAAG ./mytest.tsv .9 .5 
# python3 ./motif_find.py --alg backward ACTGGACTACGTCATGCA ./initial_emision.tsv .1 .99