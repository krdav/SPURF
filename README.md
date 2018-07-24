<p align="left">
  <img src="https://img.shields.io/docker/automated/krdav/spurf.svg" />
  <img src="https://img.shields.io/docker/build/krdav/spurf.svg" />
</p>


### Substitution Profiles Using Related Families (SPURF)

This code repository contains a command line implementation of SPURF, that takes a single antibody heavy chain DNA sequence and returns its inferred substitution profile and a logo plot of this.
SPURF uses cached data from a large-scale Rep-Seq dataset as input to a statistical model made to determine a detailed clonal family specific substitution profile for a single input sequence.
Source code to fit the SPURF model from scratch using another dataset is also provided.
Results and methods are described in our [preprint](https://arxiv.org/abs/1802.06406).
The dataset used in the paper is available in our [Zenodo bucket](https://doi.org/10.5281/zenodo.1289984).


### Cloning this repo

Clone this GitHub repo recursively to get the necessary submodules:
```
git clone --recursive https://github.com/krdav/SPURF.git
cd SPURF
git pull --recurse-submodules https://github.com/krdav/SPURF.git
```


### Installation

There are two supported ways of installing the command line implementation of SPURF: 1) using Conda on Linux and 2) using Docker and the provided Dockerfile.
The Conda installation has been tested on our own servers and a fresh Ubuntu installation on a VirtualBox.
Using VirtualBox, SPURF can be installed on both Mac and Windows.
Alternatively, Docker can also be used on any platform that supports it.


#### Using Conda

First, [install Conda](https://conda.io/docs/user-guide/install/linux.html) for Python 2.
Miniconda is sufficient and much faster at installing.
Remember to `source ~/.bashrc` if continuing installing in the same terminal window.

Install dependencies with `apt-get`:
```
sudo apt-get update
sudo apt-get upgrade -y
sudo apt-get install -y libz-dev cmake scons libgsl0-dev libncurses5-dev libxml2-dev libxslt1-dev mafft hmmer
```

Use the INSTALL executable to install the required python environment and partis (via `./INSTALL`).
After installation, the Conda environment needs to be loaded every time before use, like this:
```
source activate SPURF
```


#### Using Docker

First [install Docker](https://docs.docker.com/engine/installation/).

We have a Docker [image on Docker Hub](https://hub.docker.com/r/krdav/spurf/) that is automatically kept up to date with the master branch of this repository.
It can be pulled and used directly:
```
sudo docker pull krdav/spurf
```

Alternatively you can build the container yourself from inside the main repository directory:
```
sudo docker build -t spurf .
```

To run this container, use a command such as (see modifications below)

```
sudo docker run -it -v host-dir:/host krdav/spurf /bin/bash
```

* replace `host-dir` with the local directory to which you would like access inside your container 
* replace `/host` with the place you would like this directory to be mounted
* if you built your own container, use `spurf` in place of `krdav/spurf`

Detach using `ctrl-p ctrl-q`.


### Running SPURF

SPURF is wrapped into an Rscript named `run_SPURF.R` that takes three inputs: 

1. an antibody heavy chain DNA sequence
2. (optional) the basename for the two output files which are a substitution profile and a logo plot
3. the model type (i.e. `l2` or `jaccard`).

Example run:
```
Rscript --vanilla run_SPURF.R <input_sequence> <output_base> <model_type>
```
E.g.:
```
Rscript --vanilla run_SPURF.R CGCAGGACTGTTGANGCCTTCGGAGACCCTGTCCCTCACCTGCGTTGTCTCTGGCGGGTCCTTCAGTGATTACTACTGGAGCTGGATCCATCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCAATCATAGTGGGAGCACCAACTACAACCCGTCCCTCGAAAGTCGAGCCACCATATCAGTAGACACGTCCCAGAACAACCTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACTCGGCTGTGTATTACTGTGCGAGAGGCCCGACTACAATGGCTCACGACTTTGACTACTGGGGCCAGGGAACCCTGGTCACC seqXYZ_SPURF_output l2
```
By default, the model type `l2` is used.


### Output example

![Output logo plot](/output_examples/logo_output_example.png)

Zooming in on CDR2 and its flanking frameworks:
![Output logo plot cut](/output_examples/logo_output_example_cut.png)

