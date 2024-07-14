# ZCON

ZCON: A Zsh-driven Multifunctional Library for Genomics Analysis.

## Installation

You can clone ZCON with `git clone`.

```zsh
git clone https://github.com/tamasakian/zcon.git
```

## Dependencies

### CommandLineTools

- BLAST+
- datasets
- LAST
- MAFFT
- MCScanX
- MEME
- MUMmer4
- rapidNJ
- RAxML Next Generation
- trimAl

### PythonPackages

- BIOTP
- JCVI

If you use `venv` to manage packages; 

```zsh
## prepare venv
cd zcon
python3 -m venv venv
. venv/bin/activate

## install BIOTP
cd ..
git clone https://github.com/tamasakian/biotp.git
cd zcon
pip3 install ../biotp

## install JCVI
pip3 install jcvi
```

## Design

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

## Usage

```
zsh job/example.zsh
```