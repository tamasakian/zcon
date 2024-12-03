# ZCON

ZCON: A Zsh-driven Multifunctional Library for Genomics Analysis.

## Contents

1. [Installation instructions](#installation)
2. [Dependencies](#dependencies)
3. [Directory structure](#directory)
4. [Usage](#usage)

## Installation

You can clone ZCON with `git clone`.

```zsh
git clone https://github.com/tamasakian/zcon.git
```

## Dependencies

### CommandLineTools

- BLAST+
- [datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
- [DIAMOND](https://github.com/bbuchfink/diamond/wiki)
- [gs2](https://github.com/MotomuMatsui/gs)
- [HMMER3](https://github.com/EddyRivasLab/hmmer)
- LAST
- [MAFFT](https://github.com/GSLBiotech/mafft)
- MCScanX
- MEME
- MUMmer4
- [OrthoFinder3](https://github.com/davidemms/OrthoFinder)
- [rapidNJ](https://github.com/somme89/rapidNJ)
- [RAxML Next Generation](https://github.com/amkozlov/raxml-ng/wiki)
- [SonicParanoid2](https://gitlab.com/salvo981/sonicparanoid2)
- trimAl

### PythonPackages

- [BIOTP](https://github.com/tamasakian/biotp)
- [FASP](https://github.com/tamasakian/fasp)
- [JCVI](https://github.com/tanghaibao/jcvi/wiki)

If you use `venv` to manage packages; 

```zsh
cd zcon
python3 -m venv venv
source venv/bin/activate
pip3 install git+https://github.com/tamasakian/biotp.git
pip3 install git+https://github.com/tamasakian/fasp.git
pip3 install jcvi
deactivate
```

## Directory

### Overview

```
zcon/
├── data/
│   └── input_files
├── doc/
│   └── archive_files
├── download/
│   └── raw_data
├── job/
│   ├── profile.zsh
│   ├── example.zsh
│   └── job_files
├── lib/
│   └── function_files
├── script/
│   ├── data/
│   │   └── input_files
│   ├── doc/
│   │   └── output_files
│   └── job_files
├── taskfile/
│   └── output_files
└── venv/
```

#### DATA

Manage the input data.

#### DOC

Archive the most important files of the output data.

#### DOWNLOAD

Manage the raw data.

#### JOB

Manage the actual zsh scripts to run.

#### LIB

Manage the zsh scripts with functions defined for Genomics Analysis.

#### SCRIPT

Perform individual analysis or visualization.

#### TASKFILE

Manage the output data.

#### VENV

Manage the python packages.


## Usage

```
zsh job/example.zsh
```