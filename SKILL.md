---
name: chip-atlas
aliases:
  - chip atlas
  - atlas-dna
  - atlas-genetics
  - /chipatlas
  - /atlas_report
description: Standalone Atlas/23andMe raw DNA analyzer that generates a Russian HTML health plan and can optionally download large GWAS/ClinVar references for a deeper local report. No secrets, no external upload of DNA data.
metadata:
  clawdbot:
    emoji: "🧬"
    command: "/chipatlas"
---

# chip-atlas

Generate a **human-readable HTML genetic health plan** from an Atlas.ru / 23andMe-style raw genotype file.

The skill has two levels:

1. **Fast curated report** — portable, no downloads, embedded marker panel.
2. **Full-reference report** — explicitly offers the user to download large public GWAS/ClinVar references, explains why, gives concrete links, and then uses downloaded files to enrich the report.

## Operator behavior

When a user asks for a serious/complete report, do **not** silently stay in lightweight mode. Explain:

- default mode is quick and useful;
- full-reference mode downloads large public files;
- these files are needed to scan the raw DNA against broad GWAS/ClinVar references;
- the user's DNA stays local;
- this is educational, not medical advice.

Then offer the full workflow:

```bash
python3 bin/chip_atlas.py refs-info
python3 bin/chip_atlas.py download-refs --all --dir data
python3 bin/chip_atlas.py analyze path/to/raw_dna.txt \
  --refs-dir data \
  --out output/full_reference_report.html \
  --md output/full_reference_report.md
```

Use `--yes` only for non-interactive runs where the user already approved downloading large files:

```bash
python3 bin/chip_atlas.py download-refs --all --dir data --yes
```

## Input format

Expected TSV-like raw genotype file:

```text
# rsid	chromosome	position	genotype
rs9939609	16	53820527	AA
rs4988235	2	136608646	AG
```

The parser accepts tab or whitespace separated columns and ignores comment lines starting with `#`.

## CLI

Fast report:

```bash
python3 bin/chip_atlas.py analyze path/to/atlas_raw.txt --out output/report.html
```

Markdown too:

```bash
python3 bin/chip_atlas.py analyze path/to/atlas_raw.txt --out output/report.html --md output/report.md
```

Stats:

```bash
python3 bin/chip_atlas.py stats path/to/atlas_raw.txt
```

Reference explanation and concrete links:

```bash
python3 bin/chip_atlas.py refs-info
```

Download references:

```bash
python3 bin/chip_atlas.py download-refs --all --dir data
```

Expanded report after download:

```bash
python3 bin/chip_atlas.py analyze path/to/atlas_raw.txt \
  --refs-dir data \
  --out output/full_reference_report.html \
  --md output/full_reference_report.md
```

## Reference sources

- GWAS Catalog associations TSV  
  `https://www.ebi.ac.uk/gwas/api/search/downloads/alternative`  
  Use: broad research associations between rsID markers and published traits.

- ClinVar GRCh38 VCF  
  `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz`  
  Use: clinical-reference submissions that require doctor/lab confirmation.

These downloads can be hundreds of MB and may become >1GB as references grow. They live in `data/` and are ignored by git.

## Output style

The generated HTML follows the practical “human report” style:

- no rsID-heavy section titles;
- plain Russian explanations;
- colored labels: `опасно`, `осторожно`, `можно`, `уточнить`;
- every section says what it means and what to do;
- optional “Расширенные публичные базы” section when refs are provided;
- action checklist at the end;
- strong medical disclaimer.

## Done criteria

A good complete run should:

1. parse the raw file;
2. show SNP count and curated-marker count;
3. generate fast standalone HTML/Markdown;
4. offer full-reference download with concrete URLs and reason;
5. after download, generate expanded HTML/Markdown with GWAS/ClinVar matches;
6. avoid uploading DNA data;
7. keep medical claims cautious.

## Medical disclaimer

This is **not medical advice**, not diagnosis, and not a clinical genetic test. Consumer raw DNA files can have strand/build/provider quirks. Any medication, pregnancy, oncology, thrombosis, or treatment decision must be confirmed with a physician and a certified lab.
