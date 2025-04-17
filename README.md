# DEM Benchmark — Code Verification Competition Submission

This repository reproduces **Figures 6, 9 and 11** of the paper
> *Comparing open‑source DEM frameworks for simulations of common bulk processes* (Computer Physics Communications 296 (2024) 109066).

---

## Repository layout
```
├── code/
│   ├── dem_benchmark_plots.m   # main MATLAB script (run me!)
│   └── linspecer.m             # 9‑colour palette helper
├── data/                       # raw benchmark data (Excel files)
│   ├── ResultsRev1DrumCOM.xlsx
│   ├── ResultsRev1Penetration25K.xlsx
│   └── ResultsRev1SiloLargeM1.xlsx
├── output/                     # figures are written here by the script
├── README.md                   # this file
└── LICENSE                     # BSD 3‑Clause licence
```

## Prerequisites
* **MATLAB R2018b** or newer (tested up to R2024a)
* The **Statistics & Machine Learning Toolbox** is required for `quantile`.

If you use GNU Octave you will need to replace `xlsread` with an Octave‑compatible spreadsheet reader.

---

## Quick start
```bash
# 1. Clone the repository
git clone <this‑repo‑url>
cd <repo>/code

# 2. Launch MATLAB and run the script
matlab -batch "run('dem_benchmark_plots.m')"
# or interactively:
%>> run('dem_benchmark_plots.m')
```
The script will
1. load three `.xlsx` files from `../data/`,
2. produce three figures, and
3. save them to `../output/figure6.png`, `figure9.png`, `figure11.png`.

Estimated runtime on a typical laptop: **< 60 s**.

---

## Data source
The Excel files inside `data/` are extracted from the Zenodo record accompanying the paper:
```
Dosta M. et al. (2023) “Comparing open‑source DEM frameworks …”.
Zenodo doi:10.5281/zenodo.8252892
```
Only the three sheets needed for these plots are included here to keep the repository lightweight.

---

## Reproducing the *full* benchmarks
Running the underlying DEM simulations (25 k – 100 k particles, GPU/CPU) is **not** part of this submission.  If you wish to replicate the raw data, follow the instructions in the Zenodo dataset for each solver (MercuryDPM, Kratos, …).

---

## Citing
If you use this script or the accompanying data, please cite the original paper:
```
M. Dosta et al., Computer Physics Communications 296 (2024) 109066.
```

---

## Contact
*Maintainer:* **T Weinhart** · `<t.weinhart@utwente.nl>`

Feedback and pull requests are welcome!
