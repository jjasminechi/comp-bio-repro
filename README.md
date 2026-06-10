# Reproducing the CREMA Paper

Reproduces the figures from the paper across five species (yeast, *E. coli*, human, mouse, castor plant) using five search engines (Tide, Comet, MSGF+, MSFragger, MSAmanda), plus an extension comparing Percolator against CREMA's three peptide-level FDR methods.

---

## Setup

### Install CREMA

```bash
pip install crema-ms
```

### Install Crux (for Tide, Comet, Percolator)

Download Crux from http://crux.ms.

### Other tools

- **MSGF+**: download `MSGFPlus.zip` from https://github.com/MSGFPlus/msgfplus/releases and unzip into `datasets/MSGFPlus/`
- **MSFragger**: download from https://github.com/Nesvilab/MSFragger/wiki/Preparing-MSFragger#Downloading-MSFragger and place in `datasets/`
- **MSAmanda**: download from https://github.com/hgb-bin-proteomics/MSAmanda and place in `datasets/MSAmanda/`

---

## Datasets

All commands below assume mzML and FASTA files are placed in `datasets/`. Run all commands from `datasets/` unless noted otherwise.

| Species      | PRIDE ID  | Raw Filename                                    |
|--------------|-----------|-------------------------------------------------|
| Yeast        | PXD009420 | Tre1                                            |
| *E. coli*    | PXD011189 | 134 2018 ZBS6 Ecoli SP3 2                       |
| Human        | PXD011189 | 228 2018 ZBS6 HeLa FASP 3                       |
| Mouse        | PXD028550 | QHAGP181116_53                           |
| Castor plant | PXD007933 | Rcom 9 M4 AM R1 7Mar16 Samwise 15-08-55         |

FASTAs are downloaded from UniProt.

---

## Tide Searches (Figure 2)

Each species requires three separate Tide searches.

Replace `{sp}` with: `yeast`, `ecoli`, `human`, `mouse`, `castor`

### 1. Build index

```bash
crux tide-index \
  --output-dir {sp}-tide-index \
  --peptide-list T \
  {sp}.fasta \
  {sp}_index
```

### 2. p-value run (XCorr p-value + refactored XCorr)

```bash
crux tide-search \
  --output-dir {sp}-crux-output-pvalue \
  --exact-p-value T \
  --concat F \
  {sp}.mzML \
  {sp}_index
```

### 3. Tailor run (XCorr + Tailor score)

```bash
crux tide-search \
  --output-dir {sp}-crux-output-tailor \
  --use-tailor-calibration T \
  --concat F \
  {sp}.mzML \
  {sp}_index
```

### 4. Combined p-value run

```bash
crux tide-search \
  --output-dir {sp}-crux-output-combined \
  --exact-p-value T \
  --use-tailor-calibration T \
  --concat F \
  {sp}.mzML \
  {sp}_index
```

---

## Comet Searches (Figure 2)

```bash
crux comet \
  --output-dir {sp}-crux-output-comet-target \
  --decoy_search 0 \
  --output_txtfile 1 \
  --num_output_lines 1 \
  {sp}.mzML \
  {sp}.fasta

crux comet \
  --output-dir {sp}-crux-output-comet-decoy \
  --decoy_search 0 \
  --output_txtfile 1 \
  --num_output_lines 1 \
  {sp}.mzML \
  {sp}_decoy.fasta
```

---

## MSGF+ Searches (Figure 3)

```bash
java -Xmx8g -jar datasets/MSGFPlus/MSGFPlus.jar \
  -s {sp}.mzML \
  -d {sp}_msgf_db.fasta \
  -o {sp}_msgf.mzid \
  -t 20ppm \
  -tda 0 \
  -m 3 \
  -inst 3 \
  -thread 4

# Convert mzid to TSV
java -Xmx8g -jar datasets/MSGFPlus/MSGFPlus.jar \
  -MzIDToTsv \
  -i {sp}_msgf.mzid \
  -o {sp}_msgf.tsv \
  -showQValue 1 \
  -showDecoy 1
```

---

## MSFragger Searches (Figure 3)

```bash
java -Xmx8g -jar datasets/MSFragger-4.4.1.jar \
  closed_fragger.params \
  {sp}.mzML
```

Output: `{sp}.pepXML`

---

## MSAmanda Searches

Run from within the MSAmanda directory`:

```bash
cd datasets/MSAmanda
./MSAmanda \
  -s ../{sp}.mzML \
  -d ../{sp}.fasta \
  -e settings_crema.xml \
  -o ../{sp}_msamanda.csv
```

---

## Percolator (Extension)

Run from `datasets/`. Requires Tide searches to be completed first.

```bash
# Yeast
crux make-pin yeast-crux-output-tailor/tide-search.target.txt yeast-crux-output-tailor/tide-search.decoy.txt --output-dir yeast-percolator-tailor
crux percolator yeast-percolator-tailor/make-pin.pin --output-dir yeast-percolator-tailor
crux make-pin yeast-crux-output-combined/tide-search.target.txt yeast-crux-output-combined/tide-search.decoy.txt --output-dir yeast-percolator-combined
crux percolator yeast-percolator-combined/make-pin.pin --output-dir yeast-percolator-combined

# E. coli
crux make-pin ecoli-crux-output-tailor/tide-search.target.txt ecoli-crux-output-tailor/tide-search.decoy.txt --output-dir ecoli-percolator-tailor
crux percolator ecoli-percolator-tailor/make-pin.pin --output-dir ecoli-percolator-tailor
crux make-pin ecoli-crux-output-combined/tide-search.target.txt ecoli-crux-output-combined/tide-search.decoy.txt --output-dir ecoli-percolator-combined
crux percolator ecoli-percolator-combined/make-pin.pin --output-dir ecoli-percolator-combined

# Human
crux make-pin human-crux-output-tailor/tide-search.target.txt human-crux-output-tailor/tide-search.decoy.txt --output-dir human-percolator-tailor
crux percolator human-percolator-tailor/make-pin.pin --output-dir human-percolator-tailor
crux make-pin human-crux-output-combined/tide-search.target.txt human-crux-output-combined/tide-search.decoy.txt --output-dir human-percolator-combined
crux percolator human-percolator-combined/make-pin.pin --output-dir human-percolator-combined

# Mouse
crux make-pin mouse-crux-output-tailor/tide-search.target.txt mouse-crux-output-tailor/tide-search.decoy.txt --output-dir mouse-percolator-tailor
crux percolator mouse-percolator-tailor/make-pin.pin --output-dir mouse-percolator-tailor
crux make-pin mouse-crux-output-combined/tide-search.target.txt mouse-crux-output-combined/tide-search.decoy.txt --output-dir mouse-percolator-combined
crux percolator mouse-percolator-combined/make-pin.pin --output-dir mouse-percolator-combined

# Castor
crux make-pin castor-crux-output-tailor/tide-search.target.txt castor-crux-output-tailor/tide-search.decoy.txt --output-dir castor-percolator-tailor
crux percolator castor-percolator-tailor/make-pin.pin --output-dir castor-percolator-tailor
crux make-pin castor-crux-output-combined/tide-search.target.txt castor-crux-output-combined/tide-search.decoy.txt --output-dir castor-percolator-combined
crux percolator castor-percolator-combined/make-pin.pin --output-dir castor-percolator-combined
```

---

## Generate Figures

From the repo root:

```bash
python reproduce.py
```

This generates:
- `figure2.png` — 5×5 grid (Figure 2 reproduction + Percolator extension)
- `figure3.png` — 5×3 grid (Figure 3 reproduction)

And prints overall improvement averages:
```
OVERALL AVERAGE (all species, q-value <= 0.01)
  Avg psm-peptide vs psm-only:     17.73%  (paper: 17.14%)
  Avg psm-peptide vs peptide-only:  8.95%  (paper:  9.52%)

PERCOLATOR vs CREMA (all species, q-value <= 0.01)
  Avg percolator vs psm-only:      26.32%
  Avg percolator vs psm-peptide:    4.62%
```
