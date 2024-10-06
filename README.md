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
- datasets
- DIAMOND
- LAST
- MAFFT
- MCScanX
- MEME
- MUMmer4
- rapidNJ
- RAxML Next Generation
- SonicParanoid2
- trimAl

### PythonPackages

- BIOTP
- FASP
- JCVI

If you use `venv` to manage packages; 

```zsh
## prepare venv
cd zcon
python3 -m venv venv
source venv/bin/activate

## install BIOTP
pip3 install git+https://github.com/tamasakian/biotp.git

## install FASP
pip3 install git+https://github.com/tamasakian/fasp.git

## install JCVI
pip3 install jcvi
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