# Aligned Translation Surprisal and Source-Text Eye Movements

Code for an MSc dissertation on whether neural machine translation (NMT)
signals predict English source-text fixation durations during English--Czech
sight translation. Eye-tracking data come from the EMMT corpus (Bhattacharya
et al., 2022).

Aligned NMT surprisal, $c_\mathrm{nmt}$, is the surprisal of the NMT model's
generated Czech tokens distributed back to English source words through soft
cross-attention alignment. Monolingual surprisal, $c_\mathrm{mono}$, is summed
GPT-2 subword surprisal for each English source word.

## Research questions

- **RQ1:** Does $c_\mathrm{nmt}$ predict source-text fixation duration during
  sight translation, and how does its held-out gain compare with
  $c_\mathrm{mono}$ and six NMT attention features?
- **RQ2 (exploratory):** Is the $c_\mathrm{nmt}$ slope larger during
  translation than during prior oral reading, after lexical and positional
  controls and $c_\mathrm{mono}$ are included jointly?
- **RQ3 (exploratory):** Which part of source-text processing carries the
  overall $c_\mathrm{nmt}$ association? Conditional re-reading time is the
  theory-driven primary outcome; first-fixation duration, gaze duration, and
  go-past time are secondary first-encounter contrasts. The probability of a
  subsequent revisit is not modelled.

Because reading always precedes translation, RQ2 identifies a
translation-related stage difference within this fixed-order paradigm. It
does not isolate a translation-specific causal effect from task, repetition,
and order.

## Headline results

> **Reanalysis in progress.** The numerical values in this section were
> obtained before the July 2026 corrections to EMMT timestamp parsing and
> punctuation-sensitive SUBTLEX matching, and before trial-level correction of
> vertical gaze drift. They are retained only as a record of the previous
> analysis. The corrected Bayesian models, including the new go-past outcome,
> are being regenerated on the CSF; these values must not be treated as final.

In RQ3, conditional RRT retains its theory-driven primary status. The result
file also reports Holm-adjusted values across the three secondary first-encounter
contrasts (FFD, GD, and go-past); the figure uses those adjusted values for
the secondary outcomes and the unadjusted primary RRT value.

Within each comparison or model family, predictive models use a shared
sentence-grouped 10-fold allocation, sentence-clustered standard errors, and
sentence-level sign-flip tests. S031 and S032 are retained as distinct
sentence IDs in the mixed models but are
kept in the same fold and treated as one cluster for predictive uncertainty,
because they form a near-minimal contrastive pair. The primary predictive
analyses therefore contain 200 sentence IDs and 199 inference clusters; the
leave-pair-out sensitivity analyses contain 198 of each.

| Comparison | Delta ELPD | Clustered SE | p |
|---|---:|---:|---:|
| $c_\mathrm{nmt}$ vs controls | +14.48 | 5.65 | .038 Holm |
| $c_\mathrm{mono}$ vs controls | +4.80 | 2.77 | .042 raw; .25 Holm |
| $M_\mathrm{nmt}-M_\mathrm{mono}$, paired directly | +9.67 | 5.59 | .047 one-sided; .095 two-sided |
| $c_\mathrm{nmt}$ beyond alignment mass | +9.15 | 4.99 | .031 |
| $c_\mathrm{nmt}$ with *stoplight* retained | +10.99 | 5.61 | .027 |

The direct NMT--monolingual contrast is statistically borderline, rather than
decisive. In the RQ2 joint model, the condition-by-$c_\mathrm{nmt}$ coefficient
was 0.053 (95% credible interval [0.013, 0.092]); the implied slopes were 0.013
[-0.017, 0.044] during reading and 0.066 [0.032, 0.101] during translation.
On the translation stage, $c_\mathrm{nmt}$ also improved prediction beyond a
baseline already containing $c_\mathrm{mono}$ (Delta ELPD = +11.36,
clustered SE = 5.09, p = .015).

For RQ3, gains were -1.22 for first-fixation duration, +0.59 for gaze duration,
and +8.46 for conditional re-reading time (clustered SE = 4.54, nominal
one-sided p = .032). This is tentative evidence about re-reading duration
among revisited words, not evidence that higher $c_\mathrm{nmt}$ makes a word
more likely to be revisited. Non-significant first-encounter comparisons are not
equivalence tests. Go-past time was added to the corrected analysis and has no
pre-correction result in this table.

## Data

The EMMT recordings are available from the
[UFAL EMMT repository](https://github.com/ufal/eyetracked-multi-modal-translation)
and are not redistributed here. SUBTLEX-US supplies the frequency norms.
The corpus contains 43 participants. The released gaze files for P38 contain
column headers only and no gaze samples. P10 and P14 contain some raw gaze
data, but none of their trial-stages provides sufficient fixation coverage
across the rendered sentence for a defensible word mapping. The final derived
eye-movement files therefore contain observations from 40 participants. This
exclusion is determined by the trial-level geometric quality checks, not by
the timestamp correction.
Expected analysis inputs are:

- `fixation_durations_word.csv`
- `eye_measures_word.csv`
- `fixation_durations_word_line_diagnostics.csv`
- `eye_measures_word_line_diagnostics.csv`
- `nmt_surprisal_soft_word.csv`
- `nmt_alignment_mass_word.csv`
- `monolingual_surprisal_word.csv`
- `attention_features_6_norm.csv`
- `subtlex_us.csv`

With the July 2026 line correction, the two eye-movement extractors produce
the same 19,857 word keys and identical TFD values: 10,751 READ rows and 9,106
TRANSLATE rows. The line diagnostics cover all 2,746 available stage files;
2,459 pass the geometric quality checks (1,233 READ and 1,226 TRANSLATE), of
which 190 TRANSLATE fits use the READ prior.

The main R scripts accept `--data-dir=PATH` and `--output-dir=PATH`. Model
caches and derived CSVs are deliberately not tracked.

## Reproducing feature extraction

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

python data-extraction/extract_fixation_duration.py \
  --read_dir /path/to/EMMT/preprocessed-data/gaze/Read \
  --translate_dir /path/to/EMMT/preprocessed-data/gaze/Translate \
  --sentences /path/to/EMMT/probes/Sentences.csv \
  --output /path/to/data/fixation_durations_word.csv

python data-extraction/extract_eye_measures.py \
  --read_dir /path/to/EMMT/preprocessed-data/gaze/Read \
  --translate_dir /path/to/EMMT/preprocessed-data/gaze/Translate \
  --sentences /path/to/EMMT/probes/Sentences.csv \
  --output /path/to/data/eye_measures_word.csv

python data-extraction/extract_nmt_surprisal_soft.py \
  --sentences /path/to/EMMT/probes/Sentences.csv \
  --output /path/to/data/nmt_surprisal_soft_word.csv

python data-extraction/extract_nmt_alignment_mass.py \
  --sentences /path/to/EMMT/probes/Sentences.csv \
  --output /path/to/data/nmt_alignment_mass_word.csv

python data-extraction/extract_attention_features_norm.py \
  --sentences /path/to/EMMT/probes/Sentences.csv \
  --output /path/to/data/attention_features_6_norm.csv

python data-extraction/extract_monolingual_surprisal.py \
  --sentences /path/to/EMMT/probes/Sentences.csv \
  --output /path/to/data/monolingual_surprisal_word.csv
```

The extractors default to the revisions used for the dissertation:

- `Helsinki-NLP/opus-mt-en-cs`:
  `2820c6a540ddc2b7c4ea4c95c39b3150bd3ac27e`
- `gpt2`: `607a30d783dfa663caf39e06633721c8d4cfcd7e`

The six source-side attention definitions and uniform-attention
normalisation follow `src_seq_att` in Lim et al.'s released code at commit
[`2265d5c`](https://github.com/ZhengWeiLim/pred-trans-difficulty-NMT/blob/2265d5cc875d3dc8d6ac24071d7016086a4084c8/src/attention_util.py).
This is an adaptation of their feature construction rather than an exact
replication of their full pipeline: the present study uses the pinned
English--Czech Marian model, whitespace-delimited source words, and the same
four-beam model output used for the other NMT-derived predictors. Layers and
heads are averaged only after each multi-subword word has been treated as one
segment and divided by its corresponding uniform-attention value.

The monolingual extractor aligns GPT-2 tokens through the leading-space BPE
marker (`Ġ`). This replaces the earlier offset-based implementation, which
could assign a token spanning a leading space to the preceding word.

The EMMT gaze extractors split timestamps on colons rather than fixed character
positions. This is required because hour, minute, and second fields are not
consistently zero-padded (for example, `9:18:51.2205` and `12:9:28.1865`).
Fixation bouts shorter than 20 ms are removed before word-level aggregation in
both eye-movement extractors.

The recorded vertical gaze coordinate is not treated as a fixed screen
location. Each READ trial is fitted with a robust, potentially tilted sentence
line from fixation bouts near the rendered sentence. A TRANSLATE trial uses an
independent fit when it contains a well-supported scanpath; a weak fit may use
the immediately preceding READ slope as a geometric prior. The prior-supported
fit keeps that slope fixed throughout refitting, must remain within 80 px of
the READ line, and is subjected to the same competing-band diagnostic. It is
otherwise rejected. Acceptance requires evidence distributed across the
sentence rather than a dense cluster at one screen location. An independent
fit must cover at least four word regions, except in the shortest four-word
sentences, where three regions are required; the other support, bin-coverage,
span, dispersion, and competing-band criteria are unchanged. No rejected
trial falls back to a fixed y coordinate.

Every run writes a companion `*_line_diagnostics.csv` containing the fitted
intercept, slope, band width, residual dispersion, horizontal coverage,
equal-bout and duration support, competing-mode score, the independent fit's
failure reason when a READ prior was needed, READ--TRANSLATE displacement,
WORD/OFFTEXT/UNKNOWN counts and durations, rejection reasons, and non-fatal
review flags. Each word row also records `line_fit_source`, allowing an
independent-fit-only sensitivity analysis. Review flags identify minimally
supported prior fits and the rare case in which an adaptive prior band cycles
at one boundary bout; the latter is resolved deterministically using the
intersection of the cycling memberships and their narrowest band.

Word regions are reconstructed with the experiment's Free Sans Bold font at
28 px. The repository vendors the GNU FreeFont 2012 release under
`assets/fonts/` (SHA-256 for `FreeSansBold.ttf`:
`982534a3731416a15e2756601721f26053f68bf4239011550f3dd23ce6308215`).
Fixations more than 30 px beyond the bounded sentence region are not forced
onto the first or last word. Inter-word spaces are divided at the midpoint of
the rendered gap rather than by distance to word centres.

The eye-measure extractor keeps every retained-duration bout in sequence as a
mapped word, `OFFTEXT`, or `UNKNOWN`. This prevents an intervening unmapped
fixation from joining two visits to the same word into one first-encounter
sequence. It also derives go-past time from the ordered bouts: it accumulates
fixation time from a word's first landing until, but not including, the first
later fixation to its right,
including intervening regressions to earlier words. It leaves the value missing
when no rightward crossing is observed, when the word was first encountered
through a regression from later text, or when an unmapped bout makes the path
unrecoverable; `go_past_status` records which case applies.
`first_encounter_status` separately records whether the word was first reached
progressively or only after a word to its right had already been visited.
`reread_occurrence` indicates a later return after the first visit. The legacy
`regress_in` column is retained as an alias for compatibility, but it should
not be interpreted as proving that the return came from the word's right.

SUBTLEX lookup keys and the word-length control use a common lexical form:
surrounding Unicode punctuation is removed and case is normalised, while
internal apostrophes and hyphens are retained. The original token remains
unchanged for word-region mapping and display. This prevents sentence-final
punctuation from being mistaken for an out-of-vocabulary word.

The shared contrastive pair S031/S032 is retained
in the primary derived files: these are recorded experimental trials, not the
four-item practice round reported separately in the EMMT paper. Because the
pair was included in every probe while the remaining sentences were
distributed across probes, the core models are also refitted
after excluding both sentences. This leave-pair-out analysis is an influence
check, not a different fixation-cleaning pipeline.

## Main analysis scripts

```text
RQ1/rq1_kfold_elpd.R                  eight predictors vs the shared baseline
RQ1/rq1_joint_surprisal_kfold.R        direct and bidirectional unique ELPD tests
RQ1/rq1_mass_stoplight_robustness.R   alignment-mass and stoplight checks
RQ1/rq_locus_kfold.R                  current/preceding/following c_mono check

RQ2/rq2_joint_maximal.R               joint stage-interaction model
RQ2/rq2_interaction_kfold.R            two-model interaction predictive check
RQ2/rq2_reading_cmono_validation.R     reading c_mono and position diagnostic
RQ2/rq2_stoplight_importance.R        sequential stoplight sensitivity check

RQ3/rq3_kfold_elpd.R                  total c_nmt phase-localisation CV
RQ3/rq3_gd_rrt.R                      RQ3 figure from saved model results
```

The complete directories also contain diagnostic, plotting, and alternative
specification scripts. The previous README's “known gap” no longer applies:
the exact sentence-grouped RQ1 and RQ3 scripts are now included. The scripts
listed above, rather than the retained historical alternatives, are the
authoritative implementations for the reported and CSF-submitted analyses.

## Computational environment

Feature extraction and local validation used:

- Python 3.13.3; PyTorch 2.11.0; Transformers 5.6.2; NumPy 2.4.3;
  pandas 3.0.1
- R 4.5.1; brms 2.23.0; loo 2.9.0; lme4 2.0.1; lmerTest 3.2.1;
  dplyr 1.2.1; posterior 1.7.0

The final Bayesian refits use the CSF3 R 4.4.1 module. The environment check
prints the exact installed package versions, while each batch log records its
R session and input hashes; the local R versions above should not be reported
as the final CSF model-fitting environment.

Random seed 42 fixes the primary grouped folds and sign-flip tests. All formal
sign-flip tests use 10,000 permutations and the finite-simulation correction
$(b+1)/(B+1)$.

The final sentence-grouped cross-validation fits use four chains, 4,000
iterations per chain (2,000 warm-up), `adapt_delta = 0.95`, and
`max_treedepth = 12`. The RQ2 joint maximal model uses four chains, 8,000
iterations per chain (4,000 warm-up), `adapt_delta = 0.99`, and
`max_treedepth = 15`; the longer run addresses slow mixing previously confined
to the population and participant intercepts. Sampler warnings are emitted
immediately beside the active model, and cache versions encode the longer
sampling schedule so earlier short-run caches cannot be reused silently.

### CSF3 / Slurm

The supplied jobs use the University of Manchester CSF3 R 4.4.1 module and
four CPU cores for each full analysis. Run the environment check before
submitting the long jobs:

The frozen July 2026 input archive is
`csf_analysis_inputs_linecorrected_lim6_20260721.tar.gz` (SHA-256
`d0c20dd51d072b22f09ccf5685fc1a04361c23caafd8ed74330a383c6329da38`).
Its `MANIFEST.sha256` hash is
`6baa2003930cdbea8bbe7930b29f01eca268d63bd0e3bb4bdbeffdc6f0306001`;
the six-feature attention CSV hash is
`eda8eb1d4f51ae1d208e974e1310a562445b82241bcb357e72d68cb71145f15b`.
Earlier archives contain the five-feature or pre-Lim attention extraction and
must not be used for the final run.

```bash
export DISSERTATION_DATA_DIR=/path/to/Dissertation_Data
export DISSERTATION_OUTPUT_DIR=/path/to/Dissertation_Data/results
export EXPECTED_GIT_COMMIT="$(git rev-parse HEAD)"
export EXPECTED_MANIFEST_SHA256=6baa2003930cdbea8bbe7930b29f01eca268d63bd0e3bb4bdbeffdc6f0306001

sbatch --export=ALL,DISSERTATION_REPO_DIR="$PWD" hpc/csf3_check.sbatch
squeue -u "$USER"
```

The check verifies that the corrected full inputs contain 200 sentence IDs,
including S031/S032, checks the lexical-normalised analysis counts and go-past
column, confirms that removing the pair yields the intended 198-sentence
sensitivity sample, and performs a small Stan fit.
An older input archive in which S031/S032 were removed during extraction will
fail this check and must be replaced with the regenerated full-data archive.
If it reports missing packages, install them on a compute node and rerun the
check:

```bash
sbatch hpc/csf3_install_packages.sbatch
```

An empty `squeue` result only means that the job has left the queue. Confirm
installation from Slurm accounting and the job log before running analyses:

```bash
sacct -j JOB_ID --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
tail -n 80 dissertation-r-packages-JOB_ID.out
tail -n 120 dissertation-r-packages-JOB_ID.err
```

Proceed only when the batch step reports `COMPLETED` with exit code `0:0` and
the output confirms that package installation completed.

After the check succeeds, submit the core coefficient and CV jobs:

```bash
bash hpc/submit_core_jobs.sh \
  "$DISSERTATION_DATA_DIR" "$DISSERTATION_OUTPUT_DIR"
```

The helper submits the primary analyses with S031/S032 retained and matched
leave-pair-out refits for the core RQ1--RQ3 models, the reading-stage
$c_\mathrm{mono}$ validation, and the neighbouring-word locus checks.
Sensitivity results are
written below `results/exclude_contrastive`, while model caches carry a distinct
`_exclude_contrastive` suffix. It also keeps dependent jobs in the correct
order: the three-model RQ1 comparison and the alignment-mass check reuse the
relevant RQ1 fold allocation and caches. The RQ2 predictive job fits only the
common-slope and stage-specific-slope models needed to test the interaction.
The RQ3 job treats the total $c_\mathrm{nmt}$ gain over lexical/positional
controls as its phase-locating contrast and saves coefficient models for
direction and effect size. Comparisons conditional on $c_\mathrm{mono}$ are
handled by RQ1/RQ2 rather than repeated as a requirement for RQ3. The helper
also schedules full sampling diagnostics
for the primary RQ1 coefficient, the RQ2 joint fit, and all primary and
leave-pair-out RQ3 coefficient fits.
The CV pairs use the same random-intercept structure on both sides of each
contrast, whereas the full-data coefficient models add focal random slopes.
The former estimate held-out gain and the latter estimate direction and effect
size; they should not be described as the same fitted specification.
Each job is pinned to the repository commit and input-manifest hash present at
submission; a queued job fails rather than silently running after either the
checkout or the data snapshot changes. The preflight verifies
`MANIFEST.sha256`, and each model log records the R session, whether the pair
was excluded,
and SHA-256 hashes of the main inputs. The core CV caches and the primary
coefficient/joint-model caches additionally store MD5 metadata for every input
they use; scripts that reuse those caches reject them if the files change,
even when the observation count and fold vector remain the same.

## References

- Bhattacharya, S., Kloudova, V., Zouhar, V., & Bojar, O. (2022). EMMT: A
  simultaneous eye-tracking, 4-electrode EEG and audio corpus for multi-modal
  reading and translation scenarios. *arXiv:2204.02905*.
- Carl, M. (2013). Dynamic programming for re-mapping noisy fixations in
  translation tasks. *Journal of Eye Movement Research, 6*(2), Article 5.
- Carr, J. W., Pescuma, V. N., Furlan, M., Ktori, M., & Crepaldi, D. (2022).
  Algorithms for the automated correction of vertical drift in eye-tracking
  data. *Behavior Research Methods, 54*, 287--310.
- Cohen, A. L. (2013). Software for the automatic correction of recorded eye
  fixation locations in reading experiments. *Behavior Research Methods, 45*,
  679--683.
- Lim, Z. W., Vylomova, E., Kemp, C., & Cohn, T. (2024). Predicting human
  translation difficulty with neural machine translation. *TACL*, 12,
  1479--1496.
- Lijewska, A., Chmiel, A., & Inhoff, A. W. (2022). Stages of sight
  translation: Evidence from eye movements. *Applied Psycholinguistics*, 43,
  997--1018.
- Wilcox, E. G. et al. (2023). Testing the predictions of surprisal theory in
  11 languages. *TACL*.
