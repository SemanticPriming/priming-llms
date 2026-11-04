# priming-llms

Do LLMs show priming? An item-level comparison of language-model "priming" with human semantic priming from the Semantic Priming Across Many Languages project.

## Abstract (Psychonomic Society)

**Using the Semantic Priming Across Many Languages Project to Examine if Predictive Language Models Capture Human Semantic Priming.**
Erin Michelle Buchanan (Harrisburg University of Science and Technology), Addie Clark (Harrisburg University of Science and Technology), Maria Montefinese (University of Padova)

Semantic priming, faster target recognition following a related versus unrelated cue, is a key index of human semantic memory organization. As part of the Semantic Priming Across Many Languages Project, we conducted a pilot study to determine the best method for measuring item-level priming in large language models (LLMs) that aligns with human data. Multiple models and context sentences were tested on English cue-target pairs using surprisal-based priming deltas and embedding similarity differences as candidate indices. Surprisal-based priming deltas were small, but non-zero, indicating "priming" within LLMs. However, LLM priming was uncorrelated with human priming across all models. Embedding similarity differences were near zero for most models. Notably, one embedding similarity modestly correlated with human item-level priming, but this similarity did not translate into surprisal-based facilitation. Network analysis revealed five semantic communities with identifiable hub nodes, suggesting structured but human-dissimilar organization. These findings suggest that LLMs encode structured semantic representations that do not map onto human priming; implications for extending these comparisons to additional languages will be discussed.

> The abstract above is the text as submitted, with line-break hyphenation removed and one missing word restored ("near zero for most models"). It describes the pilot analysis (archived in `archive/540_pilot/`). The talk and manuscript report the follow-up analysis, whose results differ from the pilot: see `manuscript/` and `presentation/`.

## Layout

| Folder | Contents |
|---|---|
| `data/` | Human priming data (`priming_data/`), item data, matched stimuli, subtitle counts and semantic diversity; `data_creation.Rmd`, `diversity.rmd`, `phonemizer.rmd` build them. `data/spp/` holds item files for the original Semantic Priming Project (external validation). |
| `scripts/` | Analysis code (Python): `compute_priming.py` (surprisal and per-layer similarity), `analyze.py` (split, grid, confirmation, distributions), `prepare_spp.py` and `run_spp_validation.py` (external validation), `common.py`. `functions.r` holds the R helpers used by the pilot. |
| `outputs/` | Model outputs and result tables and figures (`surprisal/`, `discovery_grid.csv`, `layer_profile.png`, `distributions.png`, ...). |
| `manuscript/` | Draft paper (`manuscript.Rmd`, papaja) and references. |
| `presentation/` | Psychonomic Society slides (`slides.qmd`, Quarto reveal.js). |
| `archive/` | The original 540 class project and its outputs. |

## Running the analysis

From the project root:

```bash
python scripts/compute_priming.py --models gpt2 roberta-base --formats bare pair_begins
python scripts/analyze.py
```

`compute_priming.py` skips any model and format that already has an output file. `analyze.py` reads everything in `outputs/surprisal/`. The discovery/confirmation split is created once and saved to `outputs/split.csv`.

To render the slides:

```bash
quarto render presentation/slides.qmd
```

## Notes

- The English human data were collected at a single, long SOA. The original Semantic Priming Project (Hutchison et al., 2013) has both short and long SOA and is used as a second validation set.
- Few-shot examples are drawn only from the discovery half of the data.
- Masked models (BERT, RoBERTa) are scored with a sentence-final period after the target. Without it, BERT puts nearly all probability on punctuation and its surprisals are meaningless; outputs from before this fix are kept in `outputs/archive_no_period/`.
