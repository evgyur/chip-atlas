#!/usr/bin/env python3
"""chip-atlas: standalone Atlas/23andMe raw DNA HTML report generator.

Educational only. Not medical advice.
"""
from __future__ import annotations

import argparse
import gzip
import html
import sys
from collections import Counter
from urllib.request import urlopen
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class Marker:
    rsid: str
    gene: str
    title: str
    genotype: str | None
    status: str
    tone: str
    meaning: str
    actions: tuple[str, ...]
    category: str


CURATED: dict[str, dict[str, str]] = {
    # Medication
    "rs1142345": {"gene": "TPMT", "risk": "T", "label": "thiopurines"},
    "rs1799853": {"gene": "CYP2C9", "risk": "T", "label": "warfarin_cyp2c9_2"},
    "rs1057910": {"gene": "CYP2C9", "risk": "A", "label": "warfarin_cyp2c9_3"},
    "rs4244285": {"gene": "CYP2C19", "risk": "A", "label": "clopidogrel"},
    "rs3892097": {"gene": "CYP2D6", "risk": "T", "label": "cyp2d6"},
    "rs1045642": {"gene": "ABCB1", "risk": "T", "label": "abcb1"},
    # Metabolic
    "rs9939609": {"gene": "FTO", "risk": "A", "label": "weight"},
    "rs7903146": {"gene": "TCF7L2", "risk": "T", "label": "glucose"},
    "rs12255372": {"gene": "TCF7L2", "risk": "T", "label": "glucose"},
    "rs1801282": {"gene": "PPARG", "risk": "G", "label": "insulin"},
    # Folate / methylation
    "rs1801133": {"gene": "MTHFR", "risk": "T", "label": "mthfr_c677t"},
    "rs1801131": {"gene": "MTHFR", "risk": "C", "label": "mthfr_a1298c"},
    "rs1801394": {"gene": "MTRR", "risk": "G", "label": "b12"},
    "rs1805087": {"gene": "MTR", "risk": "G", "label": "b12"},
    # Thrombosis / iron / cardiovascular
    "rs6025": {"gene": "F5", "risk": "T", "label": "factor_v_leiden"},
    "rs1799963": {"gene": "F2", "risk": "A", "label": "prothrombin"},
    "rs1800562": {"gene": "HFE", "risk": "A", "label": "hfe_c282y"},
    "rs1799945": {"gene": "HFE", "risk": "G", "label": "hfe_h63d"},
    "rs7412": {"gene": "APOE", "risk": "T", "label": "apoe"},
    "rs429358": {"gene": "APOE", "risk": "C", "label": "apoe"},
    "rs4343": {"gene": "ACE", "risk": "G", "label": "blood_pressure"},
    "rs699": {"gene": "AGT", "risk": "G", "label": "blood_pressure"},
    "rs5186": {"gene": "AGTR1", "risk": "C", "label": "blood_pressure"},
    # Immune / neuro / lifestyle
    "rs7574865": {"gene": "STAT4", "risk": "T", "label": "autoimmune"},
    "rs4988235": {"gene": "MCM6/LCT", "risk": "G", "label": "lactose"},
    "rs671": {"gene": "ALDH2", "risk": "A", "label": "alcohol"},
    "rs762551": {"gene": "CYP1A2", "risk": "C", "label": "caffeine"},
    "rs4680": {"gene": "COMT", "risk": "A", "label": "stress"},
    "rs6265": {"gene": "BDNF", "risk": "T", "label": "stress"},
    "rs53576": {"gene": "OXTR", "risk": "A", "label": "social_low_confidence"},
    "rs1815739": {"gene": "ACTN3", "risk": "T", "label": "sport"},
    # Skin / ancestry
    "rs12913832": {"gene": "HERC2", "risk": "A", "label": "eye_pigment"},
    "rs1426654": {"gene": "SLC24A5", "risk": "A", "label": "skin_pigment"},
    "rs16891982": {"gene": "SLC45A2", "risk": "G", "label": "skin_pigment"},
    "rs1805007": {"gene": "MC1R", "risk": "T", "label": "sun"},
    "rs3827760": {"gene": "EDAR", "risk": "A", "label": "east_asian_marker"},
    "rs17822931": {"gene": "ABCC11", "risk": "A", "label": "earwax"},
    "rs2814778": {"gene": "DARC", "risk": "C", "label": "duffy"},
}



REF_SOURCES = {
    "gwas": {
        "name": "GWAS Catalog associations TSV",
        "url": "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
        "file": "gwas_catalog.tsv",
        "size": "~150–500 MB (changes over time)",
        "why": "broad research associations between rsID markers and published traits",
    },
    "clinvar": {
        "name": "ClinVar GRCh38 VCF",
        "url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
        "file": "clinvar.vcf.gz",
        "size": "~100–300 MB compressed, larger when processed",
        "why": "clinical-reference submissions that should be confirmed with a doctor/lab",
    },
}
REF_URLS = {k: v["url"] for k, v in REF_SOURCES.items()}


def download_url(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url) as response, dest.open("wb") as f:
        total = int(response.headers.get("Content-Length", 0) or 0)
        done = 0
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)
            done += len(chunk)
            if total:
                print(f"{dest.name}: {done/1024/1024:.1f} / {total/1024/1024:.1f} MB", end="\r")
        if total:
            print()


def parse_raw(path: Path) -> dict[str, str]:
    snps: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 4:
                continue
            rsid, _chrom, _pos, genotype = parts[:4]
            if not rsid.startswith("rs"):
                continue
            genotype = genotype.strip().upper().replace("--", "")
            if genotype:
                snps[rsid] = genotype
    return snps


def genotype_has_allele(genotype: str | None, allele: str) -> bool:
    return bool(genotype and allele and allele.upper() in genotype.upper())


def build_markers(snps: dict[str, str]) -> list[Marker]:
    g = lambda rs: snps.get(rs)
    markers: list[Marker] = []

    # Medication
    tpmt = g("rs1142345")
    if tpmt:
        has = genotype_has_allele(tpmt, "T")
        markers.append(Marker("rs1142345", "TPMT", "Азатиоприн и 6-меркаптопурин", tpmt, "опасно" if has else "можно", "red" if has else "green", "Есть сигнал, что организм может плохо перерабатывать тиопуриновые препараты." if has else "По этому маркеру сильного сигнала плохой переработки тиопуринов не видно.", ("Не принимать такие препараты без отдельного фармакогенетического теста", "Перед назначением иммуносупрессантов сказать врачу", "Если препарат нужен — дозу подбирать только с врачом"), "meds"))

    cyp2c9_hits = [rs for rs in ("rs1799853", "rs1057910") if genotype_has_allele(g(rs), CURATED[rs]["risk"])]
    if g("rs1799853") or g("rs1057910"):
        markers.append(Marker("rs1799853/rs1057910", "CYP2C9", "Варфарин и сильные препараты для свёртываемости", f"rs1799853={g('rs1799853') or 'нет'}, rs1057910={g('rs1057910') or 'нет'}", "осторожно" if cyp2c9_hits else "уточнить", "yellow" if cyp2c9_hits else "blue", "Есть сигнал, который может влиять на переработку варфарина. По raw-файлу нельзя назначать дозу." if cyp2c9_hits else "По основным маркерам явного сигнала не видно, но это не полный клинический тест.", ("Не подбирать дозу по этому отчёту", "Перед варфарином — клинический фармакогенетический тест + INR", "Показать раздел врачу, если антикоагулянты реально назначают"), "meds"))

    clop = g("rs4244285")
    if clop:
        has = genotype_has_allele(clop, "A")
        markers.append(Marker("rs4244285", "CYP2C19", "Клопидогрел", clop, "осторожно" if has else "можно", "yellow" if has else "green", "Есть сигнал возможной слабой активации клопидогрела." if has else "По одному из главных маркеров плохой переработки клопидогрела опасного сигнала не видно.", ("Если препарат назначают после стента/инфаркта — обсудить CYP2C19 с врачом",), "meds"))

    # Metabolic
    fto = g("rs9939609")
    if fto:
        has = genotype_has_allele(fto, "A")
        markers.append(Marker("rs9939609", "FTO", "Вес, аппетит и тяга к калорийной еде", fto, "внимание" if has else "можно", "yellow" if has else "green", "Есть генетический сигнал склонности к набору веса; эффект хорошо компенсируется режимом." if has else "По этому маркеру выраженного сигнала склонности к набору веса не видно.", ("Белок в каждый приём пищи", "Овощи и клетчатка каждый день", "Силовые 2–3 раза в неделю + ходьба/кардио", "Ограничить сладкие напитки и ультрапереработанные снеки"), "food"))

    glucose_hits = [rs for rs in ("rs7903146", "rs12255372") if genotype_has_allele(g(rs), "T")]
    if g("rs7903146") or g("rs12255372"):
        markers.append(Marker("rs7903146/rs12255372", "TCF7L2", "Сахар и инсулин", f"rs7903146={g('rs7903146') or 'нет'}, rs12255372={g('rs12255372') or 'нет'}", "контроль" if glucose_hits else "можно", "yellow" if glucose_hits else "green", "Есть повод внимательнее относиться к углеводному обмену и риску инсулинорезистентности." if glucose_hits else "По этим маркерам выраженного сигнала по углеводному обмену не видно.", ("Сдать глюкозу натощак и HbA1c", "При лишнем весе/тяге к сладкому добавить инсулин и HOMA-IR", "Держать стабильную активность и сон"), "food"))

    lact = g("rs4988235")
    if lact:
        has = genotype_has_allele(lact, "G")
        markers.append(Marker("rs4988235", "MCM6/LCT", "Молочные продукты", lact, "по самочувствию" if has else "можно", "blue" if has else "green", "Есть промежуточный/возможный сигнал чувствительности к лактозе." if has else "Генетически лактоза во взрослом возрасте должна переноситься лучше среднего.", ("Если после молока вздутие — попробовать lactose-free", "Кефир, йогурт и сыр часто переносятся лучше", "Если симптомов нет — специально исключать молочку не нужно"), "food"))

    aldh = g("rs671")
    if aldh:
        has = genotype_has_allele(aldh, "A")
        markers.append(Marker("rs671", "ALDH2", "Алкоголь", aldh, "осторожно" if has else "можно умеренно", "yellow" if has else "green", "Есть сигнал генетической непереносимости алкоголя / накопления ацетальдегида." if has else "Типичной генетической непереносимости алкоголя по главному маркеру не видно.", ("Умеренность важнее генетики", "Не пить вместе с потенциально токсичными лекарствами", "При flush-реакции или плохом самочувствии — лучше избегать"), "food"))

    caf = g("rs762551")
    if caf:
        has = genotype_has_allele(caf, "C")
        markers.append(Marker("rs762551", "CYP1A2", "Кофеин", caf, "после обеда осторожно" if has else "можно", "yellow" if has else "green", "Профиль похож на более медленную/промежуточную переработку кофеина." if has else "Профиль похож на более быструю переработку кофеина.", ("Если сон или тревожность страдают — кофе только до обеда", "Убрать кофеин после 14:00 на неделю и сравнить сон"), "food"))

    # Folate
    m677 = g("rs1801133")
    m1298 = g("rs1801131")
    if m677 or m1298:
        has = genotype_has_allele(m677, "T") or genotype_has_allele(m1298, "C")
        markers.append(Marker("rs1801133/rs1801131", "MTHFR", "Фолаты и витамин B9", f"C677T={m677 or 'нет'}, A1298C={m1298 or 'нет'}", "уточнить" if has else "можно", "blue", "Есть мягкий/неполный сигнал по фолатному циклу. По raw-файлу нельзя делать жёсткие выводы." if has else "По найденным MTHFR-маркерам выраженного сигнала не видно или ключевой маркер отсутствует.", ("Сдать B12, фолат и гомоцистеин", "При планировании беременности дозировки фолата обсуждать с врачом", "Не пить большие дозы метилфолата вслепую"), "food"))
    else:
        markers.append(Marker("rs1801133", "MTHFR", "Фолаты и витамин B9", None, "уточнить", "blue", "Главный маркер MTHFR C677T не найден в файле, поэтому вывод неполный.", ("Сдать B12, фолат и гомоцистеин", "При необходимости проверить MTHFR клинически"), "food"))

    # Cardio / blood / iron
    f5 = g("rs6025")
    f2 = g("rs1799963")
    if f5 or f2:
        has = genotype_has_allele(f5, "T") or genotype_has_allele(f2, "A")
        markers.append(Marker("rs6025/rs1799963", "F5/F2", "Тромбозы", f"F5={f5 or 'нет'}, F2={f2 or 'нет'}", "осторожно" if has else "частые маркеры спокойные", "yellow" if has else "green", "Есть сигнал по одному из частых маркеров тромбофилии." if has else "По двум самым известным частым маркерам тромбофилии явного сигнала не видно.", ("Это не исключает все причины тромбозов", "При семейной истории/беременности/гормональной терапии обсуждать с врачом", "Не начинать антикоагулянты по raw DNA"), "heart"))

    hfe1, hfe2 = g("rs1800562"), g("rs1799945")
    if hfe1 or hfe2:
        has = genotype_has_allele(hfe1, "A") or genotype_has_allele(hfe2, "G")
        markers.append(Marker("rs1800562/rs1799945", "HFE", "Железо и ферритин", f"C282Y={hfe1 or 'нет'}, H63D={hfe2 or 'нет'}", "контролировать" if has else "можно", "blue" if has else "green", "Есть мягкий маркер, который может быть связан с обменом железа. Сам по себе это не диагноз." if has else "Сильного сигнала наследственного перегруза железом не видно.", ("Сдать ферритин, железо, трансферрин/ОЖСС", "Не принимать железо без подтверждённого дефицита"), "heart"))

    apoe1, apoe2 = g("rs7412"), g("rs429358")
    if apoe1 or apoe2:
        complete = bool(apoe1 and apoe2)
        markers.append(Marker("rs7412/rs429358", "APOE", "Холестерин и сосуды", f"rs7412={apoe1 or 'нет'}, rs429358={apoe2 or 'нет'}", "профилактика" if complete else "уточнить", "yellow" if complete else "blue", "APOE-статус частично/полностью виден, но практический вывод всё равно лучше делать по липидам." if complete else "Полный APOE-статус определить нельзя: один из двух нужных маркеров отсутствует.", ("Сдать липидограмму, ApoB и Lp(a)", "Измерять давление дома", "Держать сон, активность, клетчатку и омега‑3 по показаниям"), "heart"))

    # Sport/stress/skin/ancestry
    actn = g("rs1815739")
    if actn:
        markers.append(Marker("rs1815739", "ACTN3", "Тип нагрузки", actn, "универсально" if actn == "CT" else "подобрать", "green", "Профиль не определяет судьбу в спорте. Чаще всего лучше работает смесь силы и выносливости.", ("Силовые 2–3 раза в неделю", "Кардио/ходьба 150 минут в неделю", "Главная цель — стабильность, мышцы и чувствительность к инсулину"), "sport"))

    comt = g("rs4680")
    if comt:
        has = genotype_has_allele(comt, "A")
        markers.append(Marker("rs4680", "COMT", "Стресс-чувствительность", comt, "особенность" if has else "обычно", "yellow" if has else "green", "Есть профиль, который часто связывают с более острой реакцией на стресс и стимуляторы." if has else "По этому маркеру выраженного slow-COMT сигнала не видно.", ("Сон 7–8 часов", "Осторожнее с кофеином, предтрениками и сильными стимуляторами", "Дыхание/прогулка/растяжка при остром стрессе"), "stress"))

    skin_hits = [rs for rs in ("rs1426654", "rs16891982", "rs12913832", "rs1805007") if g(rs)]
    if skin_hits:
        markers.append(Marker("pigmentation", "SLC24A5/SLC45A2/HERC2/MC1R", "Кожа и солнце", ", ".join(f"{rs}={g(rs)}" for rs in skin_hits), "SPF", "yellow", "Маркеры пигментации помогают понять вероятную чувствительность к солнцу, но реальный фенотип важнее SNP.", ("SPF при активном солнце", "Не загорать до ожога", "Если много родинок или семейная история — дерматоскопия по графику врача"), "skin"))

    ancestry_hits = [rs for rs in ("rs3827760", "rs17822931", "rs2814778") if g(rs)]
    if ancestry_hits:
        markers.append(Marker("ancestry", "EDAR/ABCC11/DARC", "Происхождение", ", ".join(f"{rs}={g(rs)}" for rs in ancestry_hits), "инфо", "blue", "Это отдельные ancestry-маркеры, а не полноценный ancestry-тест.", ("EDAR A часто связан с восточноазиатским компонентом", "ABCC11 и DARC дают отдельные популяционные подсказки", "Не использовать как точную карту происхождения"), "ancestry"))

    auto = g("rs7574865")
    if auto and genotype_has_allele(auto, "T"):
        markers.append(Marker("rs7574865", "STAT4", "Аутоиммунные ассоциации", auto, "только при симптомах", "blue", "Есть ассоциация с некоторыми аутоиммунными заболеваниями, но это не диагноз и не повод скрининговать всё подряд.", ("При сухости глаз/рта, утренней скованности суставов, сыпях — к ревматологу", "Без симптомов не делать панель анализов из страха"), "heart"))

    return markers


def stats_for(snps: dict[str, str]) -> dict[str, object]:
    chroms = Counter()
    # second pass not available; stats based on curated only + total. CLI parse is rsid->gt only by design.
    return {"total": len(snps), "curated_found": sum(1 for rs in CURATED if rs in snps)}


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _split_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in info.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v.replace("_", " ").replace("%2C", ",")
    return out


def load_reference_hits(snps: dict[str, str], refs_dir: Path, max_hits: int = 80) -> dict[str, list[dict[str, str]]]:
    """Best-effort local matching of raw rsIDs against downloaded public references."""
    rsids = set(snps)
    hits: dict[str, list[dict[str, str]]] = {"gwas": [], "clinvar": []}

    gwas = refs_dir / REF_SOURCES["gwas"]["file"]
    if gwas.exists():
        try:
            with _open_text(gwas) as f:
                header = f.readline().rstrip("\n").split("\t")
                idx = {name: i for i, name in enumerate(header)}
                snp_i = idx.get("SNPS")
                trait_i = idx.get("DISEASE/TRAIT")
                gene_i = idx.get("MAPPED_GENE")
                risk_i = idx.get("STRONGEST SNP-RISK ALLELE")
                p_i = idx.get("P-VALUE")
                if snp_i is not None:
                    for line in f:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) <= snp_i:
                            continue
                        found = [rs.strip() for rs in parts[snp_i].replace(";", ",").split(",") if rs.strip() in rsids]
                        if not found:
                            continue
                        hits["gwas"].append({
                            "rsid": ", ".join(found[:3]),
                            "genotype": ", ".join(f"{rs}={snps.get(rs, '')}" for rs in found[:3]),
                            "trait": parts[trait_i] if trait_i is not None and len(parts) > trait_i else "",
                            "gene": parts[gene_i] if gene_i is not None and len(parts) > gene_i else "",
                            "risk": parts[risk_i] if risk_i is not None and len(parts) > risk_i else "",
                            "p": parts[p_i] if p_i is not None and len(parts) > p_i else "",
                        })
                        if len(hits["gwas"]) >= max_hits:
                            break
        except Exception as e:
            hits["gwas"].append({"rsid": "error", "trait": f"Could not parse GWAS reference: {e}", "genotype": "", "gene": "", "risk": "", "p": ""})

    clinvar = refs_dir / REF_SOURCES["clinvar"]["file"]
    if clinvar.exists():
        try:
            with _open_text(clinvar) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 8:
                        continue
                    ids = [x for x in parts[2].split(";") if x in rsids]
                    if not ids:
                        continue
                    info = _split_info(parts[7])
                    hits["clinvar"].append({
                        "rsid": ", ".join(ids[:3]),
                        "genotype": ", ".join(f"{rs}={snps.get(rs, '')}" for rs in ids[:3]),
                        "trait": info.get("CLNDN", ""),
                        "gene": info.get("GENEINFO", ""),
                        "significance": info.get("CLNSIG", ""),
                    })
                    if len(hits["clinvar"]) >= max_hits:
                        break
        except Exception as e:
            hits["clinvar"].append({"rsid": "error", "trait": f"Could not parse ClinVar reference: {e}", "genotype": "", "gene": "", "significance": ""})

    return hits


def render_reference_section(ref_hits: dict[str, list[dict[str, str]]] | None) -> str:
    if not ref_hits:
        return ""
    blocks = []
    gwas = ref_hits.get("gwas") or []
    if gwas:
        rows = []
        for h in gwas[:40]:
            rows.append(f"<li><b>{html.escape(h.get('trait','')[:120])}</b> — {html.escape(h.get('rsid',''))} ({html.escape(h.get('genotype',''))}); gene: {html.escape(h.get('gene','')[:80])}; risk: {html.escape(h.get('risk','')[:80])}</li>")
        blocks.append(f"<h3>GWAS Catalog: исследовательские ассоциации</h3><ul>{''.join(rows)}</ul>")
    clinvar = ref_hits.get("clinvar") or []
    if clinvar:
        rows = []
        for h in clinvar[:40]:
            rows.append(f"<li><b>{html.escape(h.get('trait','')[:120])}</b> — {html.escape(h.get('rsid',''))} ({html.escape(h.get('genotype',''))}); gene: {html.escape(h.get('gene','')[:80])}; ClinVar: {html.escape(h.get('significance','')[:80])}</li>")
        blocks.append(f"<h3>ClinVar: совпадения, которые требуют проверки</h3><ul>{''.join(rows)}</ul>")
    if not blocks:
        return ""
    return f'''
<div class="card">
  <h2>🔬 Расширенные публичные базы</h2>
  <div class="tip">Это не диагноз. Это локальные совпадения raw DNA с публичными справочниками. Их задача — дать материал для глубокого разбора и список тем, которые можно подтвердить с врачом/лабораторией.</div>
  {''.join(blocks)}
</div>'''


def pill_class(tone: str) -> str:
    return {"red": "pill-red", "yellow": "pill-yellow", "green": "pill-green", "blue": "pill-blue"}.get(tone, "pill-blue")


def render_item(m: Marker) -> str:
    genotype = f'<div class="tip">Маркер: {html.escape(m.gene)} · {html.escape(m.rsid)} · генотип: {html.escape(m.genotype or "не найден")}</div>'
    lis = "\n".join(f"<li>{html.escape(a)}</li>" for a in m.actions)
    return f'''
  <div class="item">
    <h3>{html.escape(m.title)} <span class="pill {pill_class(m.tone)}">{html.escape(m.status)}</span></h3>
    <p>{html.escape(m.meaning)}</p>
    <ul>
      {lis}
    </ul>
    {genotype}
  </div>'''


def render_section(title: str, markers: Iterable[Marker]) -> str:
    items = "\n".join(render_item(m) for m in markers)
    if not items:
        return ""
    return f'''
<div class="card">
  <h2>{html.escape(title)}</h2>
{items}
</div>'''


def build_checklist() -> list[str]:
    return [
        "Сохранить раздел «Лекарства» и показать врачу при назначении новых препаратов",
        "Сдать: глюкоза натощак, HbA1c, липидограмма",
        "Добавить: ApoB и Lp(a) — для сосудистого риска",
        "Сдать: B12, фолат, гомоцистеин, ферритин",
        "Железо + трансферрин/ОЖСС — если ферритин высокий или низкий",
        "Проверить домашнее давление тонометром",
        "На 7 дней убрать кофеин после 14:00 и посмотреть сон",
        "Собрать простой режим: белок + овощи + ежедневная ходьба",
        "При планировании беременности обсудить фолат, железо и тромбозный анамнез с врачом",
    ]


def render_html(snps: dict[str, str], markers: list[Marker], ref_hits: dict[str, list[dict[str, str]]] | None = None) -> str:
    now = datetime.now().strftime("%Y-%m-%d")
    by_cat: dict[str, list[Marker]] = {k: [] for k in ["meds", "food", "heart", "sport", "stress", "skin", "ancestry"]}
    for m in markers:
        by_cat.setdefault(m.category, []).append(m)
    checklist = "\n".join(f"<li>☐ {html.escape(x)}</li>" for x in build_checklist())
    return f'''<!DOCTYPE html>
<html lang="ru">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>🧬 Генетический план здоровья</title>
<style>
* {{ box-sizing: border-box; margin: 0; padding: 0; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background: #0a0e1a; color: #e2e8f0; line-height: 1.7; padding: 16px; }}
.wrap {{ max-width: 700px; margin: 0 auto; }}
.hero {{ text-align: center; padding: 24px 0 32px; border-bottom: 2px solid #3b82f6; margin-bottom: 28px; }}
.hero h1 {{ font-size: 1.8em; color: #60a5fa; margin-bottom: 6px; }}
.hero p {{ color: #94a3b8; font-size: 1em; }}
.card {{ background: #111827; border-radius: 14px; padding: 22px; margin-bottom: 20px; }}
.card h2 {{ font-size: 1.25em; color: #60a5fa; margin-bottom: 14px; display: flex; align-items: center; gap: 10px; }}
.pill {{ display: inline-block; padding: 3px 10px; border-radius: 20px; font-size: 0.72em; font-weight: 700; text-transform: uppercase; letter-spacing: 0.4px; }}
.pill-red {{ background: #7f1d1d; color: #fecaca; }}
.pill-yellow {{ background: #713f12; color: #fde68a; }}
.pill-green {{ background: #14532d; color: #bbf7d0; }}
.pill-blue {{ background: #1e3a8a; color: #bfdbfe; }}
.item {{ background: #1e293b; border-radius: 10px; padding: 16px; margin: 12px 0; }}
.item h3 {{ color: #e2e8f0; font-size: 1.05em; margin-bottom: 8px; }}
.item p {{ color: #cbd5e1; margin-bottom: 10px; }}
.item ul {{ margin-left: 18px; }}
.item li {{ margin: 7px 0; color: #cbd5e1; }}
.tip {{ background: #0f172a; border-left: 3px solid #3b82f6; padding: 10px 14px; margin-top: 10px; border-radius: 0 8px 8px 0; font-size: 0.92em; color: #94a3b8; }}
.footer {{ text-align: center; color: #64748b; font-size: 0.82em; padding: 20px 0 40px; }}
</style>
</head>
<body>
<div class="wrap">

<div class="hero">
  <h1>🧬 Генетический план здоровья</h1>
  <p>Простыми словами: питание, спорт, лекарства и что проверить</p>
  <p>{len(snps):,} SNP · найдено ключевых маркеров: {sum(1 for rs in CURATED if rs in snps)} · {now}</p>
</div>

{render_section('💊 Лекарства — важно знать врачу', by_cat.get('meds', []))}
{render_section('🍽️ Питание и метаболизм', by_cat.get('food', []))}
{render_section('❤️ Сердце, сосуды и кровь', by_cat.get('heart', []))}
{render_section('🏃 Спорт', by_cat.get('sport', []))}
{render_section('🧠 Стресс, сон и нервная система', by_cat.get('stress', []))}
{render_section('☀️ Кожа и солнце', by_cat.get('skin', []))}
{render_section('🧬 Происхождение', by_cat.get('ancestry', []))}
{render_reference_section(ref_hits)}

<div class="card">
  <h2>📋 Что сделать сегодня</h2>
  <div class="item"><ul>{checklist}</ul></div>
</div>

<div class="footer">
  Этот отчёт создан в образовательных целях и не заменяет консультацию врача.<br>
  Генетика — предрасположенность, а не приговор. Consumer raw DNA может содержать strand/build/provider quirks.<br>
  Важные решения по лекарствам, беременности и лечению подтверждать клиническим тестом и врачом.
</div>

</div>
</body>
</html>
'''


def render_markdown(snps: dict[str, str], markers: list[Marker]) -> str:
    lines = ["# 🧬 Генетический план здоровья", "", f"SNP в файле: **{len(snps):,}**", f"Ключевых маркеров найдено: **{sum(1 for rs in CURATED if rs in snps)}**", "", "> Образовательный отчёт, не медицинское заключение.", ""]
    for m in markers:
        lines += [f"## {m.title} — {m.status}", "", f"- Ген/маркер: `{m.gene}` / `{m.rsid}`", f"- Генотип: `{m.genotype or 'не найден'}`", f"- Что значит: {m.meaning}", "", "Что делать:"]
        lines += [f"- {a}" for a in m.actions]
        lines.append("")
    lines += ["## Чек-лист", ""] + [f"- [ ] {x}" for x in build_checklist()]
    return "\n".join(lines) + "\n"


def cmd_stats(args: argparse.Namespace) -> int:
    snps = parse_raw(Path(args.input))
    print(f"SNPs: {len(snps):,}")
    print(f"Curated markers found: {sum(1 for rs in CURATED if rs in snps)} / {len(CURATED)}")
    return 0


def cmd_analyze(args: argparse.Namespace) -> int:
    src = Path(args.input)
    if not src.exists():
        print(f"error: input not found: {src}", file=sys.stderr)
        return 2
    snps = parse_raw(src)
    markers = build_markers(snps)
    ref_hits = None
    if args.refs_dir:
        refs_dir = Path(args.refs_dir)
        if refs_dir.exists():
            ref_hits = load_reference_hits(snps, refs_dir, max_hits=args.refs_max)
            print(f"Reference hits: {sum(len(v) for v in ref_hits.values())} (from {refs_dir})")
        else:
            print(f"warning: refs dir not found: {refs_dir}", file=sys.stderr)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render_html(snps, markers, ref_hits=ref_hits), encoding="utf-8")
    print(f"HTML report: {out}")
    if args.md:
        md = Path(args.md)
        md.parent.mkdir(parents=True, exist_ok=True)
        md.write_text(render_markdown(snps, markers), encoding="utf-8")
        print(f"Markdown report: {md}")
    return 0



def print_refs_info() -> None:
    print("Full-reference mode downloads public databases for a deeper local report.")
    print("Why: the built-in panel is practical; GWAS/ClinVar add broad research/clinical-reference matches.")
    print("Privacy: these are public references; your DNA file stays local and is not uploaded.")
    print()
    for key, meta in REF_SOURCES.items():
        print(f"{key}: {meta['name']}")
        print(f"  URL:  {meta['url']}")
        print(f"  File: {meta['file']}")
        print(f"  Size: {meta['size']}")
        print(f"  Use:  {meta['why']}")
        print()


def cmd_refs_info(args: argparse.Namespace) -> int:
    print_refs_info()
    print("Download everything:")
    print("  python3 bin/chip_atlas.py download-refs --all --dir data --yes")
    print("Then generate expanded report:")
    print("  python3 bin/chip_atlas.py analyze raw.txt --refs-dir data --out output/full_report.html --md output/full_report.md")
    return 0


def cmd_download_refs(args: argparse.Namespace) -> int:
    """Download optional public reference files for advanced/local experimentation."""
    data_dir = Path(args.dir)
    print_refs_info()
    jobs = []
    if args.gwas or args.all:
        jobs.append(("gwas", data_dir / REF_SOURCES["gwas"]["file"]))
    if args.clinvar or args.all:
        jobs.append(("clinvar", data_dir / REF_SOURCES["clinvar"]["file"]))
    if not jobs:
        print("Nothing selected. Use --all, --gwas, or --clinvar. Use refs-info to see concrete links.")
        return 2
    if not args.yes:
        print("This can download hundreds of MB or more. Continue? [y/N] ", end="")
        answer = input().strip().lower()
        if answer not in {"y", "yes", "д", "да"}:
            print("cancelled")
            return 1
    for key, dest in jobs:
        if dest.exists() and not args.force:
            print(f"skip existing: {dest}")
            continue
        print(f"Downloading {key}: {REF_URLS[key]}")
        download_url(REF_URLS[key], dest)
        print(f"saved: {dest}")
    print()
    print("Next step:")
    print(f"  python3 bin/chip_atlas.py analyze path/to/raw_dna.txt --refs-dir {data_dir} --out output/full_report.html --md output/full_report.md")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="chip_atlas.py", description="Generate a human-readable HTML report from Atlas/23andMe raw DNA TSV.")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("stats", help="Print raw file stats")
    s.add_argument("input")
    s.set_defaults(func=cmd_stats)
    a = sub.add_parser("analyze", help="Generate HTML report")
    a.add_argument("input")
    a.add_argument("--out", "-o", default="output/chip_atlas_report.html")
    a.add_argument("--md", help="Optional Markdown output path")
    a.add_argument("--refs-dir", default=None, help="Optional directory with downloaded GWAS/ClinVar references for expanded local annotation")
    a.add_argument("--refs-max", type=int, default=80, help="Maximum hits per reference source to include in the report")
    a.set_defaults(func=cmd_analyze)

    ri = sub.add_parser("refs-info", help="Explain optional large references, concrete links, sizes, and full-report workflow")
    ri.set_defaults(func=cmd_refs_info)

    r = sub.add_parser("download-refs", help="Download optional public GWAS/ClinVar references for expanded reports")
    r.add_argument("--dir", default="data", help="Destination directory")
    r.add_argument("--all", action="store_true", help="Download all optional references")
    r.add_argument("--gwas", action="store_true", help="Download GWAS Catalog TSV")
    r.add_argument("--clinvar", action="store_true", help="Download ClinVar VCF.gz")
    r.add_argument("--force", action="store_true", help="Overwrite existing files")
    r.add_argument("--yes", "-y", action="store_true", help="Do not ask for confirmation before large downloads")
    r.set_defaults(func=cmd_download_refs)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
