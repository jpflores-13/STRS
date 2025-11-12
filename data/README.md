# data/

This directory organizes **raw inputs** (FASTQs and metadata) and **processed outputs** used in analyses and figures.

---
data/
├── raw/ # Not tracked in Git; downloaded from SRA/GEO (see raw/README)
└── processed/ # Derived files used by analyses & figures
├── hic/ # .hic, loops
├── rna/ # count matrices, DE results
└── cutntag/ # bigWigs, peaks
---

## Access

- Raw data acquired from **SRA/GEO** (accessions in `raw/README.md`).
- Processed files generated with scripts

See the subdirectory READMEs for formats and how to regenerate.
