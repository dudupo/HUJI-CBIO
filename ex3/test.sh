# python3 ZOOPS_EM.py ./seqs_CCGG.fasta CCGG 0.05 0.9 0.1 0.1
rm generated_tests.sh
python3  -W ignore create_zoops_test.py ./seqs_CCGG.fasta AAAAAAA 0.3 0.99 0.0 1e-10
# python3  -W ignore create_zoops_test.py ./seqs_CCGG.fasta CCGG 0.05 0.1 0.1 0.1
# python3  -W ignore create_zoops_test.py ./seqs_CCGG.fasta CCGG 0.05 0.2 0.1 0.1
# python3  -W ignore create_zoops_test.py ./seqs_CCGG.fasta CCGG 0.05 0.3 0.1 0.1

# python3  -W ignore create_zoops_test.py ./seqs_CCGG.fasta CCGG 0.01 0 0.1 0.1
# python3  -W ignore create_zoops_test.py ./seqs_CCGG.fasta CCGG 0.02 0.1 0.1 0.1
# python3  -W ignore create_zoops_test.py ./seqs_CCGG.fasta CCGG 0.03 0.2 0.1 0.1
# python3  -W ignore create_zoops_test.py ./seqs_CCGG.fasta CCGG 0.07 0.3 0.1 0.1

# python3  -W ignore create_zoops_test.py ./seqs_CCGG.fasta CCGG 0.01 0.3 0.1 0.1

cat ./generated_tests.sh 
chmod 777 ./generated_tests.sh
./generated_tests.sh