# Custom SurVirus filtering

`CUSTOM_FILTER_SURVIRUS` runs after `RUN_SURVIRUS`, reading the complete
`results.remapped.txt` and the host/virus breakpoint FASTAs and sdust BEDs.
It never uses the stock filtered calls as input. It writes `accepted.tsv`,
`rejected.tsv`, `all_candidates.tsv`, and `summary.json` under
`<outputDir>/<sample>/survirus/<sample>_custom_filter/`.

## Parameters

| Nextflow parameter | Default | Meaning |
|---|---:|---|
| `filter_min_host_pbs` | 0.8 | Minimum remapped HOST_PBS |
| `filter_min_pairs` | 2 | Supporting pairs sufficient without split evidence |
| `filter_split_min_pairs` | 1 | Alternative minimum pairs when split evidence is present |
| `filter_min_split_reads` | 1 | Minimum split reads for that alternative |
| `filter_min_seq_length` | 30 | Minimum sequence length on each side |
| `filter_mask_mode` | reject | `reject` or `flag` high masked fractions |
| `filter_max_masked_fraction` | 0.8 | Flag/reject strictly above this fraction |
| `filter_pairing_rescue` | true | Rescue a rejected partner of a passing call |
| `filter_pair_min_distance` | -50 | Minimum signed reverse-host minus forward-host position |
| `filter_pair_max_distance` | 1000 | Maximum signed distance |
| `filter_deduplicate` | true | Remove later sequence duplicates |
| `filter_dedup_fraction` | 0.8 | Minimum aligned query span fraction on both sides |

Support means `pairs >= 2 OR (pairs >= 1 AND split_reads >= 1)` by default.
The two counts may overlap biologically; this does not prove two independent
DNA molecules. Coverage is reported but never used to reject or rescue.
NaN/infinite HOST_PBS is recorded as `INVALID_HOST_PBS`.

Both modes calculate masking to provide auditable flags. `flag` disables
mask-based rejection. BED intervals are merged and counted as zero-based,
half-open intervals (`end-start`). Unlike stock SurVirus, all intervals count,
without overlap double counting or an extra base. `REPEAT_FLAGS` identifies
sides above the configured fraction. These are sdust low-complexity flags,
not overlaps with a genomic RepeatMasker annotation. Missing BEDs are errors;
an existing empty BED means no masked bases.

Pairing requires the same host chromosome and virus contig, opposite host
strands and opposite virus strands, and the configured host distance interval.
Passing calls are paired first, then each still-unpaired original passing call
can rescue one nearest rejected call. Ties follow input order. Rescue can
override any initial failure, including masking; `STATUS=rescued` and
`INITIAL_FILTER_FAILURES` preserve that fact. Rescued calls cannot start rescue
chains. `PAIRED_WITH` records the pre-dedup partner, which may later be removed.

Deduplication uses SurVirus's striped Smith-Waterman scoring (1,4,6,1), with
reverse complement alignment when strands differ. As in the stock code, the
test is `(query_end-query_begin+1)/shorter_sequence_length >= threshold` on BOTH
sides, not percent identity. Candidates follow input order; each is compared
against retained earlier calls. `DUPLICATE_OF` identifies the retained call.
The tiny C++ wrapper compiles against the installed SurVirus `libs/` sources;
gcc and g++ must be available. Python 2.7 and Python 3 are supported.

`REJECTION_REASON` describes final rejection. `INITIAL_FILTER_FAILURES` also
remains populated for rescued or duplicate calls. Empty results still have
TSV headers. Raw SurVirus breakpoint strings retain their original coordinates.

## Run on existing N10 outputs

```bash
nextflow run main.nf \
  --filter_workdir /data/xrz/capint/nftide-survirus/work/e4/01e92d7a3f63912919b2a13a27b027/N10_survirus \
  --filter_sample N10
```

This executes only the custom filter, without reference indexing or alignment.
For repeat flagging, add `--filter_mask_mode flag`. All other parameters also
apply to normal pipeline runs. Use `-output-dir` to keep alternative runs apart.

## N10 verification

September 16 correction: inclusive SSW endpoints now count the full aligned
span; empty sequences receive short-sequence rejection rather than aborting.
Deduplication disabled no longer requires loading a shared library.
`HOST_HAS_MASKED_BASES` and `VIRUS_HAS_MASKED_BASES` flag any mask overlap,
while `REPEAT_FLAGS` continues to indicate fractions above the threshold.
`RESCUED_BY` retains the rescue sponsor ID. Removed calls have no active pair
link; a rescue whose sponsor was removed is rejected with
`RESCUE_PARTNER_REMOVED`. Initial failures remain available separately.
The earlier table below records the original implementation's validation;
corrected results are in `tests/n10_validation/fixed_reject/` and
`tests/n10_validation/fixed_flag/`.

Both configurations were executed through Nextflow 25.10.2 with the existing
SurVirus environment and actual N10 files (411 candidates).

| Mode | Initial passes | Final retained | Final rejected |
|---|---:|---:|---:|
| reject | 297 | 269 | 142 |
| flag | 341 | 297 | 114 |

No N10 calls were rescued. Results are retained under
`tests/n10_validation/published/N10/survirus/N10_custom_filter/` and
`tests/n10_validation/flag_published/N10/survirus/N10_custom_filter/`.
Separate small tests cover overlapping mask intervals, exact thresholds,
zero coverage, one-pair/one-split support, pairing rescue, disabling rescue,
deduplication, and rejection audit fields.
