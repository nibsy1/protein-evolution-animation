---
name: sequence-scraper
description: Scrapes one homologous protein sequence per model plant species (16
  species across 5 lineages) from NCBI Protein and UniProt for a gene name (ECA3),
  AGI code (AT1G10130), or protein name. Progressive cascade: exact isoform lookup,
  all-family collection, k-mer similarity verification, term broadening to protein
  class, UniProt gap-fill. Outputs FASTA + TSV with quality flags.
version: 4
---
# Sequence Scraper (final)

Scrapes exactly one homologous protein sequence per model plant species
(16 species across 5 lineages) for a query: gene symbol (ECA3), AGI code
(AT1G10130), or protein name.

## Usage

```bash
python scripts/scrape_sequences.py ECA3 -o ECA3.fasta --table ECA3.tsv --email you@example.com
```

## Species list (16)

Green Algae: Chlamydomonas reinhardtii, Volvox carteri, Chlorella vulgaris
Bryophytes: Physcomitrium patens, Marchantia polymorpha
Ferns: Azolla filiculoides
Gymnosperms: Pinus taeda, Picea abies
Angiosperms: Arabidopsis thaliana, Oryza sativa, Zea mays,
Solanum lycopersicum, Populus trichocarpa, Phaseolus vulgaris,
Hordeum vulgare, Nicotiana benthamiana

(Salvinia cucullata and Ginkgo biloba were removed: their NCBI/UniProt
data is too sparse — Salvinia has only chloroplast records, Ginkgo almost
nothing. 100% coverage requires species with real proteome data.)

## How it works

1. **Query resolution** — gene symbols resolve via UniProt exact-gene
lookup in Arabidopsis (ECA3 -> Q9SY55/AT1G10130, never ECA1); AGI codes
resolve via NCBI. The resolved record pins the correct isoform and
provides the reference sequence/length.
2. **Per-species collection** — NCBI Protein search tries gene symbol,
AGI, family protein name (species then genus scope), collecting ALL
family isoforms.
3. **Closest-match selection** — candidates are verified by k-mer
similarity to the query sequence (threshold 0.008; distant algal
homologs pass at ~0.017, junk scores ~0.0); the closest one is kept.
4. **Term broadening (last resort)** — if nothing passed, the query
widens to the protein class (e.g. "ATPase"), fetches up to 100
records, and the similarity check identifies the true homolog among
unnamed/loosely-annotated records. This fills gaps like Picea abies
(unnamed protein product) or Marchantia (P-type Ca2+ transporter).
5. **UniProt gap-filling** — species still without a hit are filled
from UniProt REST, same verification.
6. **Fragment rescue** — if the best hit is <60% of query length, a
full-length close relative from the same genus replaces it.

## Output quality flags

query | exact | homolog | broad | uniprot | genus_fill | related/fallback | miss

## Validation

- ECA3: 16/16 species (~376 s; broadening adds time only for gaps)
- rbcL: 16/16 (fast path)

## Dependencies

Python 3.8+, `requests`.
