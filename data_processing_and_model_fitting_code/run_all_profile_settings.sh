

# Run with all settings on Symphogen data:
#python make_profile.py --glob_pat "sdata/*partitions-cluster-annotations.csv" --nproc 5 --dataset_name symphogen --totN_sub 50 --nsubs 200
#python make_profile.py --glob_pat "sdata/*partitions-cluster-annotations.csv" --nproc 5 --dataset_name symphogen --totN_sub 50 --nsubs 1000
#python make_profile.py --glob_pat "sdata/*partitions-cluster-annotations.csv" --nproc 5 --dataset_name symphogen --totN_sub 10 --nsubs 200
#python make_profile.py --glob_pat "sdata/*partitions-cluster-annotations.csv" --nproc 5 --dataset_name symphogen --totN_sub 10 --nsubs 1000
#python make_profile.py --glob_pat "sdata/*partitions-cluster-annotations.csv" --nproc 16 --dataset_name symphogen --totN_sub 1 --nsubs 200 --MAX_REP 30
#    python make_profile.py --glob_pat "sdata/*partitions-cluster-annotations.csv" --nproc 10 --dataset_name symphogen --totN_sub 1 --nsubs 1000

# Run with all settings on public data:
#python make_profile.py --glob_pat "public_data/*partitions-cluster-annotations.csv" --nproc 5 --dataset_name public --totN_sub 50 --nsubs 200
#python make_profile.py --glob_pat "public_data/*partitions-cluster-annotations.csv" --nproc 5 --dataset_name public --totN_sub 50 --nsubs 1000
#python make_profile.py --glob_pat "public_data/*partitions-cluster-annotations.csv" --nproc 5 --dataset_name public --totN_sub 10 --nsubs 200
#python make_profile.py --glob_pat "public_data/*partitions-cluster-annotations.csv" --nproc 5 --dataset_name public --totN_sub 10 --nsubs 1000
#python make_profile.py --glob_pat "public_data/*partitions-cluster-annotations.csv" --nproc 16 --dataset_name public --totN_sub 1 --nsubs 200 --MAX_REP 30 --allele_finding true
#python make_profile.py --glob_pat "public_data/*partitions-cluster-annotations.csv" --nproc 16 --dataset_name public --totN_sub 1 --nsubs 500 --MAX_REP 30 --allele_finding true

# Run the juno data:
python make_profile.py --glob_pat "juno_data/*partitions-cluster-annotations.csv" --nproc 16 --dataset_name juno_data --minN_mixing 28 --totN_sub 1 --nsubs 100 --MAX_REP 30 --allele_finding true


# Move everything to an ordered folder:
mv *aammp_profiles*.txt profiles/

