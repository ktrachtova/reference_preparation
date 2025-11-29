# Reference preparation repository

[![Bash](https://img.shields.io/badge/Bash-4EAA25?logo=gnubash&logoColor=fff)](#)
[![R](https://img.shields.io/badge/R-%23276DC3.svg?logo=r&logoColor=white)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=fff)](https://www.python.org/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

This repository contains code and data for reference preparation of various projects. Current projects:
* [PerSeqPIPE](https://github.com/ktrachtova/perseqpipe): bioinformatic pipeline for analysis of small non-coding RNAs using sequence-centric approach

## PerSeqPIPE

* `docker/` folder contains Dockerfile and Makefile that builds docker image used to prepare various databases and subsequently final custom sncRNA GTF file for PerSeqPIPE pipeline
* `scripts/` folder contains simple scripts that prepare various databases for PerSeqPIPE pipeline as described [here](https://github.com/ktrachtova/perseqpipe/blob/main/docs/reference_databases.md)