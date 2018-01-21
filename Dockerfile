#FROM python:2.7-slim
FROM continuumio/miniconda


RUN conda create -y -n SPURF
RUN source activate SPURF

RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y libz-dev cmake scons libgsl0-dev libncurses5-dev libxml2-dev libxslt1-dev mafft hmmer

WORKDIR /spurf
ADD . /spurf
RUN ./INSTALL_docker


