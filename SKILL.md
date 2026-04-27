---
name: chip-atlas
aliases:
  - chip atlas
  - atlas-dna
  - atlas-genetics
  - /chipatlas
  - /atlas_report
description: Standalone Atlas/23andMe raw DNA analyzer that generates a simple human-readable HTML health plan from a TSV genotype file. No secrets, no external APIs.
metadata:
  clawdbot:
    emoji: "🧬"
    command: "/chipatlas"
---

# chip-atlas

Generate a **human-readable HTML genetic health plan** from an Atlas.ru / 23andMe-style raw genotype file.

This skill is intentionally portable and public:

- no API keys;
- no private databases;
- no network calls;
- no personal data bundled;
- Python standard library only.

## Input format

Expected TSV-like raw genotype file:

```text
# rsid	chromosome	position	genotype
rs9939609	16	53820527	AA
rs4988235	2	136608646	AG
```

The parser accepts tab or whitespace separated columns and ignores comment lines starting with `#`.

## CLI

```bash
python3 bin/chip_atlas.py analyze path/to/atlas_raw.txt --out output/report.html
```

Optional Markdown summary:

```bash
python3 bin/chip_atlas.py analyze path/to/atlas_raw.txt --out output/report.html --md output/report.md
```

Print file statistics only:

```bash
python3 bin/chip_atlas.py stats path/to/atlas_raw.txt
```

## Output style

The generated HTML follows the practical “human report” style:

- no rsID-heavy section titles;
- plain Russian explanations;
- colored labels: `опасно`, `осторожно`, `можно`, `уточнить`;
- every section says what it means and what to do;
- action checklist at the end;
- strong medical disclaimer.

## Medical disclaimer

This is **not medical advice**, not diagnosis, and not a clinical genetic test. Consumer raw DNA files can have strand/build/provider quirks. Any medication, pregnancy, oncology, thrombosis, or treatment decision must be confirmed with a physician and a certified lab.

## Done criteria

A good run should:

1. parse the raw file;
2. show SNP count and chromosome coverage;
3. detect key curated markers when present;
4. generate a standalone `.html` report that opens in a browser;
5. avoid exposing or storing secrets.
