# zcon

ZCON: A Zsh-driven Multifunctional Library for Genomics Analysis.

## Dependencies

CommandLineTools

- blast
- datasets
- mafft
- meme
- mummer4
- raxml-ng
- trimal

PythonPackages

- jcvi

## DirectoryStructureDiagram

```
.
└── zcon/
    ├── data/
    │   └── input_files
    ├── doc/
    │   └── archive_files
    ├── download/
    │   └── raw_data
    ├── job/
    │   ├── example.zsh
    │   └── profile.zsh 
    ├── lib/
    │   └── hoge.zsh 
    ├── script/
        ├── data/
        │   └── input_files
        ├── doc/
        │   └── output_files
        └── hoge.R
    ├── taskfile/
    │   └── output_files
    ├── venv/
    ├── .gitignore
    └── README.md
```