### Substitution Profiles Using Related Families (SPURF)

This code repository contains a command line implementation of SPURF, that takes a single antibody heavy chain DNA sequence and returns its inferred Position-Specific Scoring Matrix (PSSM) and a logo plot of this PSSM.
SPURF uses cached data from a large-scale Rep-Seq dataset as input to a statistical model made to determine a detailed clonal family specific substitution profile for a single input sequence.
Source code to fit the SPURF model from scratch using another dataset is also provided.


### Cloning this repo

Clone this GitHub repo recursively to get the necessary submodules:
```shell
git clone --recursive https://github.com/krdav/SPURF.git
cd SPURF
git pull --recurse-submodules https://github.com/krdav/SPURF.git
```

### Installation

There are two supported ways of installing the command line implementation of SPURF 1) using conda on Linux and 2) using Docker and the provided Dockerfile.
The conda installation has been tested on our own servers and a fresh Ubuntu installation on a VirtualBox.
Using VirtualBox SPURF can be installed on both MAC and Windows.
Alternatively, Docker can also be used on multiple platforms.


#### Using conda

Installation conda using the guide here: https://conda.io/docs/user-guide/install/linux.html

Use the INSTALL executable to install the required python environment and partis (via `./INSTALL`). Notice that HMMER3 is required in the PATH since this is a dependency of ANARCI (for AHo annotation). Installing partis may require extra attention. We have tested the installation on a fresh Ubuntu installation and running the following will satisfy the extra partis requirements:
```
sudo apt-get update
sudo apt-get upgrade -y
sudo apt-get install -y libz-dev cmake scons libgsl0-dev libncurses5-dev libxml2-dev libxslt1-dev mafft hmmer
```

With Docker:
```
cd SPURF
docker build -t spurf .
sudo docker run -it spurf bash
```

Detach using `ctrl-p ctrl-q`


After installation, the conda environment needs to be loaded every time before use, like this:
```shell
source activate SPURF
```

Then, use the Rscript to infer a substitution profile for your input sequence:
```
Rscript --vanilla run_SPURF.R <input_sequence> <output_base>
E.g.:
Rscript --vanilla run_SPURF.R CGCAGGACTGTTGANGCCTTCGGAGACCCTGTCCCTCACCTGCGTTGTCTCTGGCGGGTCCTTCAGTGATTACTACTGGAGCTGGATCCATCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCAATCATAGTGGGAGCACCAACTACAACCCGTCCCTCGAAAGTCGAGCCACCATATCAGTAGACACGTCCCAGAACAACCTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACTCGGCTGTGTATTACTGTGCGAGAGGCCCGACTACAATGGCTCACGACTTTGACTACTGGGGCCAGGGAACCCTGGTCACC seqXYZ_SPURF_output
```

The output is a PSSM and a logo plot of this PSSM.


