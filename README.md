Pull this GitHub repo recursively to get the necessary submodules:
```shell
git clone --recursive https://github.com/krdav/SPURF.git
git pull --recurse-submodules https://github.com/krdav/SPURF.git
```

Use the INSTALL executable to install the required python environment and partis (via `./INSTALL`).

After installation, the conda environment needs to be loaded every time before use, like this:
```shell
cd SPURF
source activate SPURF
```

Then, use the Rscript to infer a substitution profile for your input sequence:
```
Rscript --vanilla run_SPURF.R <input_sequence> <output_base>
E.g.:
Rscript --vanilla run_SPURF.R CGCAGGACTGTTGANGCCTTCGGAGACCCTGTCCCTCACCTGCGTTGTCTCTGGCGGGTCCTTCAGTGATTACTACTGGAGCTGGATCCATCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCAATCATAGTGGGAGCACCAACTACAACCCGTCCCTCGAAAGTCGAGCCACCATATCAGTAGACACGTCCCAGAACAACCTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACTCGGCTGTGTATTACTGTGCGAGAGGCCCGACTACAATGGCTCACGACTTTGACTACTGGGGCCAGGGAACCCTGGTCACC seqXYZ_SPURF_output
```

The output is a PSSM and a logo plot of this PSSM.


