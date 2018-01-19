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
The output is a single row matrix consisting of the predicted amino acid frequnencies at each AHo position.
Each entry in the matrix is named using the `p_X_a_Y` format, where `X` refers to the AHo position and `Y` represents the amino acid index at that position.
The amino acid ordering used here is `('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-')`, where `'-'` indicates a gap character.

