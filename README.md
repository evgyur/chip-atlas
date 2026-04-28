# chip-atlas

Public OpenClaw skill + standalone CLI for generating a Russian HTML health plan from Atlas.ru / 23andMe-style raw genotype files.

It has two modes:

1. **Fast human report** — works immediately from a curated marker panel.
2. **Full-reference report** — offers to download large public GWAS/ClinVar references and adds local reference matches to the report.

## What it does

`chip-atlas` reads a raw DNA TSV file like:

```text
# rsid	chromosome	position	genotype
rs9939609	16	53820527	AA
rs7903146	10	114758349	CT
```

and creates a practical browser report:

- medication cautions;
- metabolism / weight / glucose notes;
- folate / B-vitamin checks;
- thrombosis common-marker notes;
- iron / HFE notes;
- lactose, alcohol, caffeine;
- stress, sport, skin/sun, ancestry markers;
- optional GWAS Catalog / ClinVar reference matches;
- concrete checklist.

No API keys, no external upload of your DNA data.

## Install

```bash
git clone https://github.com/evgyur/chip-atlas.git
cd chip-atlas
python3 bin/chip_atlas.py --help
```

No dependencies beyond Python 3.10+ standard library.

## Quick use: human report

Generate HTML:

```bash
python3 bin/chip_atlas.py analyze ~/Downloads/atlas_raw.txt --out output/my_report.html
```

Generate HTML + Markdown:

```bash
python3 bin/chip_atlas.py analyze ~/Downloads/atlas_raw.txt \
  --out output/my_report.html \
  --md output/my_report.md
```

Show stats only:

```bash
python3 bin/chip_atlas.py stats ~/Downloads/atlas_raw.txt
```

## Full-reference mode: download the big public databases

The curated report is useful, but it is intentionally small. If you want a richer report like a serious local research pass, download public reference databases.

First see what will be downloaded and why:

```bash
python3 bin/chip_atlas.py refs-info
```

Current references:

- **GWAS Catalog associations TSV**  
  URL: `https://www.ebi.ac.uk/gwas/api/search/downloads/alternative`  
  Why: broad research associations between rsID markers and published traits.

- **ClinVar GRCh38 VCF**  
  URL: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz`  
  Why: clinical-reference submissions that should be confirmed with a doctor/lab.

These can be hundreds of MB and may become >1GB as sources grow or if you later add more references. They are saved under `data/` and ignored by git.

Download all references:

```bash
python3 bin/chip_atlas.py download-refs --all --dir data
```

Non-interactive version:

```bash
python3 bin/chip_atlas.py download-refs --all --dir data --yes
```

Then generate a full report using them:

```bash
python3 bin/chip_atlas.py analyze ~/Downloads/atlas_raw.txt \
  --refs-dir data \
  --out output/full_report.html \
  --md output/full_report.md
```

What happens after download:

- your raw DNA file is still processed locally;
- `chip-atlas` scans downloaded GWAS/ClinVar files for rsIDs present in your raw file;
- the HTML gets an extra **“Расширенные публичные базы”** section;
- the section is deliberately cautious: it gives research/clinical-reference matches, not diagnoses.

## Suggested “good report” workflow

```bash
# 1) Check the raw file
python3 bin/chip_atlas.py stats ~/Downloads/atlas_raw.txt

# 2) Generate the fast readable report
python3 bin/chip_atlas.py analyze ~/Downloads/atlas_raw.txt \
  --out output/human_report.html \
  --md output/human_report.md

# 3) Download public references for deeper local matching
python3 bin/chip_atlas.py refs-info
python3 bin/chip_atlas.py download-refs --all --dir data

# 4) Generate the expanded report
python3 bin/chip_atlas.py analyze ~/Downloads/atlas_raw.txt \
  --refs-dir data \
  --out output/full_reference_report.html \
  --md output/full_reference_report.md
```

Open the HTML in a browser and use the Markdown as a compact text version for follow-up analysis.

## Privacy

The tool runs locally. It does not upload your DNA data anywhere. Downloaded references are public files; your genotype file stays on your machine.

## Disclaimer

This is educational software, not a medical device. Do not change medication, supplements, pregnancy plans, or medical treatment based only on this report. Confirm important findings with a doctor and a certified clinical lab.
