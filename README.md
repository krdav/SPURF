Pull this GitHub repo recursively to get the necessary submodules:
```shell
git clone --recursive https://github.com/krdav/AMP.git
git pull --recurse-submodules https://github.com/krdav/AMP.git
```

Use the INSTALL executable to install the required python environment and partis (via `./INSTALL`).

After installation, the conda environment needs to be loaded every time before use, like this:
```shell
source activate AMP
```

Then, enter into an interactive R session and input your test sequence:
```R
# This is the input sequence.
input_seq = "CGCAGGACTGTTGANGCCTTCGGAGACCCTGTCCCTCACCTGCGTTGTCTCTGGCGGGTCCTTCAGTGATTACTACTGGAGCTGGATCCATCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCAATCATAGTGGGAGCACCAACTACAACCCGTCCCTCGAAAGTCGAGCCACCATATCAGTAGACACGTCCCAGAACAACCTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACTCGGCTGTGTATTACTGTGCGAGAGGCCCGACTACAATGGCTCACGACTTTGACTACTGGGGCCAGGGAACCCTGGTCACC"
```

Predict the clonal family substitution profile corresponding to the above input sequence as follows:
```R
source("AMP.R")
pred.prof = predict.prof(input_seq)
```
