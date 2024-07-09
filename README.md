# zcon

ZCON: A Zsh-driven Multifunctional Library for Genomics Analysis.

## Dependencies

### CommandLineTools

- blast
- datasets
- mafft
- meme
- mummer4
- raxml-ng
- trimal

### PythonPackages

- jcvi
- biotp

If you use `venv` to manage packages; 

```zsh
## prepare venv
python3 -m venv venv
. venv/bin/activate

## install jcvi
pip3 install jcvi

## install biotp
cd ..
git clone https://github.com/tamasakian/biotp.git
cd zcon
pip3 install ../biotp
```

## Design

```
zcon/
├── data/
│   └── input_files
│
├── doc/
│   └── archive_files
│
├── download/
│   └── raw_data
│
├── job/
│   ├── profile.zsh
│   │
│   ├── example.zsh
│   │
│   └── job_files
│
├── lib/
│   └── function_files
│
├── script/
│   ├── data/
│   │   └── input_files
│   │
│   ├── doc/
│   │   └── output_files
│   │
│   └── job_files
│
├── taskfile/
│   └── output_files
│
└── venv
```

## Usage
```
zsh job/example.zsh
```