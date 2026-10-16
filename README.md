# wgdlife_analysis

Figure code of *The life history of whole-genome doubling across human cancers* (Zhang et al., 2026): the panels of
Figs. 2–5 that visualise the analyses of the [wgdlife](https://github.com/pigudog/wgdlife) package. Each `figureN/figNx.R`
reads the result tables in `figureN/data/` and writes the panel to `results/figureN/figNx.pdf` and `.png`.

## Run

R 4.4.3. `environment/install.R` installs the package versions the panels were rendered with (only the ones that are
missing; they, and `remotes` if absent, come from CRAN). The panels use the Arial font, which Windows and macOS ship; on
Debian/Ubuntu it is installed by `sudo apt install ttf-mscorefonts-installer`. Without Arial the scripts fall back to
Liberation Sans (package `fonts-liberation`), which has the same glyph widths, so the layout is unchanged.

```bash
Rscript environment/install.R
bash run.sh              # every panel -> results/
bash run.sh figure2      # one figure
```

`run.sh` changes to the repository root first; a single panel is rendered from there too (`Rscript figure2/fig2b.R`).

## Tables

The tables in `figureN/data/` are the study's own results: the values each panel draws and the identifiers of its rows,
aggregated over loci, events, cancer types or cohorts. No third-party release is included. The one per-tumour table is
`figure3/data/fig3b_hrd_onset.tsv` (Fig. 3b draws one point per tumour): without identifiers, the study's own estimate of
the onset of the HRD-associated mutational process relative to the WGD time in each of the 95 PCAWG tumours.

## Licence

MIT. Please cite the article when reusing the tables.
