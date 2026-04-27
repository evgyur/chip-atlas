# chip-atlas

Public OpenClaw skill + standalone CLI for generating a simple Russian HTML health plan from Atlas.ru / 23andMe-style raw genotype files.

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
- concrete checklist.

No API keys, no external services, no private data.

## Install

```bash
git clone https://github.com/evgyur/chip-atlas.git
cd chip-atlas
python3 bin/chip_atlas.py --help
```

No dependencies beyond Python 3.10+ standard library.

## Use

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

## Privacy

The tool runs locally. It does not upload your DNA data anywhere.

## Disclaimer

This is educational software, not a medical device. Do not change medication, supplements, pregnancy plans, or medical treatment based only on this report. Confirm important findings with a doctor and a certified clinical lab.
