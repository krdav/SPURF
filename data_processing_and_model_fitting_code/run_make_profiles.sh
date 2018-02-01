
# Run the public data:
python make_profile.py --glob_pat "public_data/*partitions-cluster-annotations.csv" --nproc 16 --dataset_name public --totN_sub 1 --nsubs 500 --MAX_REP 30 --allele_finding true

# Run the juno data:
python make_profile.py --glob_pat "juno_data/*partitions-cluster-annotations.csv" --nproc 16 --dataset_name juno_data --minN_mixing 28 --totN_sub 1 --nsubs 100 --MAX_REP 30 --allele_finding true

