### Substitution Profiles Using Related Families (SPURF)

This code repository contains a command line implementation of SPURF, that takes a single antibody heavy chain DNA sequence and returns its inferred Position-Specific Scoring Matrix (PSSM) and a logo plot of this PSSM.
SPURF uses cached data from a large-scale Rep-Seq dataset as input to a statistical model made to determine a detailed clonal family specific substitution profile for a single input sequence.
Source code to fit the SPURF model from scratch using another dataset is also provided.


### Cloning this repo

Clone this GitHub repo recursively to get the necessary submodules:
```
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

Install conda for Python2 using the guide here: https://conda.io/docs/user-guide/install/linux.html
Miniconda is sufficient and much faster at installing.
Remember to `source ~/.bashrc` if continuing installing in the same terminal window.


Install dependencies with `apt-get`:
```
sudo apt-get update
sudo apt-get upgrade -y
sudo apt-get install -y libz-dev cmake scons libgsl0-dev libncurses5-dev libxml2-dev libxslt1-dev mafft hmmer
```

Use the INSTALL executable to install the required python environment and partis (via `./INSTALL`).
After installation, the conda environment needs to be loaded every time before use, like this:
```
source activate SPURF
```


#### Using Docker

Install Docker following the guide here: https://docs.docker.com/engine/installation/

Inside the repository folder build the container and run it:
```
sudo docker build -t spurf .
sudo docker run -it spurf bash
```

Detach using `ctrl-p ctrl-q`.

We also have a docker image on Docker Hub that can be pulled and used directly:
```
sudo docker pull krdav/spurf
sudo docker run -it -v /:/host krdav/spurf /bin/bash
```



### Running SPURF

SPURF is wrapped into an Rscript name `run_SPURF.R` that takes three inputs 1) an antibody heavy chain DNA sequence, and optionally, 2) the basename for the two output files which are a PSSM and a logo plot, and 3) the model type (i.e. `l2` or `jaccard`).
Example of a run:
```
Rscript --vanilla run_SPURF.R <input_sequence> <output_base> <model_type>
E.g.:
Rscript --vanilla run_SPURF.R CGCAGGACTGTTGANGCCTTCGGAGACCCTGTCCCTCACCTGCGTTGTCTCTGGCGGGTCCTTCAGTGATTACTACTGGAGCTGGATCCATCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCAATCATAGTGGGAGCACCAACTACAACCCGTCCCTCGAAAGTCGAGCCACCATATCAGTAGACACGTCCCAGAACAACCTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACTCGGCTGTGTATTACTGTGCGAGAGGCCCGACTACAATGGCTCACGACTTTGACTACTGGGGCCAGGGAACCCTGGTCACC seqXYZ_SPURF_output l2
```
By default, the model type `l2` is used.


