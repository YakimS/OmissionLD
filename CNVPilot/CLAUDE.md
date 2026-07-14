# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Auditory **CNV (Contingent Negative Variation)** experiment. Each trial is a *couple*:
a warning stimulus **S1** followed, after a fixed S1-S2 interval (ISI), by an imperative
stimulus **S2**. The anticipatory CNV is a slow negative potential that builds during the
S1-S2 interval. Design is based on the mixed-couple paradigm of Gauthier et al. (1985)
(summarised in `../CNV_paradigm_summary_and_sleep_mapping.docx`). Built for EEG with
PsychToolbox audio and optional EGI NetStation triggers. Sibling of `../OmissionPilot`.

## Paradigm

**S1 is acoustically identical across both couples**, so the couple type is unknown to the
subject until S2. Two couples (interleaved within a block):

- `same` — S1 -> S1 (same tone repeated). **TARGET**: respond ASAP with an eye movement.
- `diff` — S1 -> S2 (a different tone). No response.

The anticipatory CNV in the S1-S2 window is, by construction, the same for both couples;
they diverge only at/after S2 (target detection vs frequency change).

**Response** is by **eye movement only**, analysed offline from EOG. Nothing is collected
online (no button/keyboard response). The subject watches a **silent nature movie played
externally** (this code does audio + triggers only; no PsychToolbox `Screen`). Instruct the
subject the eye movement is *reactive to S2* so the pre-S2 terminal CNV stays movement-free.

## Running the Experiment

**Phase 1 - Generate + validate configuration:**
```matlab
cd CNVPilot
config_CNV_exp                        % writes Config-CNVPilot.mat
analyze_CNV_config('Config-CNVPilot.mat')   % validates counts/timing/metadata (asserts)
```

**Phase 2 - Run:**
```matlab
cd CNVPilot
run_CNV_exp                           % loads config, plays audio, sends triggers
```
Set `isLabComputer = true` in `run_CNV_exp.m` for PsychToolbox + NetStation; `false` uses
MATLAB's built-in `sound()` for a desktop dry run (no NetStation, `input()` breaks).
Lab controls: SPACE (any key) advances a self-paced break; ESC during a rest aborts.

## Key Parameters (`config_CNV_exp.m`)

```matlab
exp_meta.STIM_DURATION     = 0.2;      % tone duration (S1 and S2), seconds
exp_meta.S1_FREQ           = 2000;     % S1 (warning) tone; also S2 of couple 'same'
exp_meta.S2_FREQ           = 500;      % S2 tone of couple 'diff'
exp_meta.SWAP_FREQS        = false;    % swap S1/S2 between subjects (counterbalance)
exp_meta.ISI_VALUES        = [0.8 1 1.2];% S1-S2 intervals (s); held CONSTANT within a block
exp_meta.N_BLOCKS_PER_ISI  = 2;        % blocks per ISI -> n_total_blocks = this x numel(ISI)
exp_meta.REST_MEAN         = 6.0;      % rest after S2 before next S1 (s)
exp_meta.REST_JITTER       = 3.0;      % -> rest uniform in [MEAN-JIT, MEAN+JIT] (unpredictable S1)
exp_meta.N_TRIALS          = 86;       % test trials PER ISI value (block type)
exp_meta.PERCENT_DIFF      = 20;       % percent of test trials that are 'diff' (rest = 'same'/target)
exp_meta.N_WARMUP          = 2;        % leading target trials per block (excluded)
```
`N_TRIALS` test trials are generated **per ISI value** (block type) and split as evenly as
possible across that ISI's `N_BLOCKS_PER_ISI` blocks, each `PERCENT_DIFF` 'diff' / rest
'same' (grand total = `N_TRIALS` x numel(`ISI_VALUES`)). The config prints trial counts and
a duration estimate; `analyze_CNV_config` validates them. **No two 'diff' trials are ever adjacent**
(they are the minority; 'same' may repeat). **Frequencies are fixed for the whole session**
(the target rule "S2 == S1" must be stable); only swap S1/S2 between subjects via `SWAP_FREQS`.

## Experimental Design

- **ISI is blocked** (constant within a block): a fixed, predictable S1-S2 rhythm is what
  lets the subject time S2 and build the late/expectancy CNV. Couples are randomly
  interleaved within a block.
- **`N_BLOCKS_PER_ISI` x numel(`ISI_VALUES`) blocks** (default 2 x 3 = 6); block order
  shuffled so same-ISI blocks are not adjacent (avoids time-on-task confounding the ISI
  contrast).
- Each block: `N_WARMUP` warmup `same` trials (re-establish the rhythm) + the test trials
  (`same`/`diff` interleaved, its share of `N_TRIALS`), with **no two `diff` trials adjacent**.

## Timing (onset-to-onset)

Per trial the onset-to-onset ITI is `ISI + STIM + REST`, identical for both couples. In the
events table this is encoded as: `wait_after(S1) = ISI - STIM`; `wait_after(S2) = REST`. The
runner (lab) also schedules S2 at an **absolute deadline** `S1_onset + ISI` (sample-accurate
`PsychPortAudio('Start', ..., when)`) so the critical S1->S2 interval never drifts.

## Data Flow

```
config_CNV_exp.m -> exp_meta struct -> Config-CNVPilot.mat -> run_CNV_exp.m -> audio + triggers
                                                    \-> analyze_CNV_config.m (validation)
```
`exp_meta` contains: waveforms (`beep_S1`, `beep_S2`), `events_table`, block ISI order,
all parameters, `params_string`, and duration estimates.

## Events Table

One row per event, **two rows per trial** (S1 then S2). Columns: `block, trial,
event_in_trial (1/2), couple, trial_kind (warmup/test), is_target, isi, rest, event_type
(tone), freq, s1_freq, s2_freq, wait_after, trigger, condition_name, display_str, ctyp`.

## NetStation Triggers

- `STRT` / `SEND` — experiment start/end (`SEND` with `abrt=1` if aborted).
- `BGIN` — block begin (`BNUM`, `CNAM`, `ISI`, `PRMS`).
- `WARN` — S1 (warning tone), sent on every trial at the true S1 onset.
- `S2SM` — S2 of couple `same` (target) at the scheduled S2 onset.
- `S2DF` — S2 of couple `diff` (different tone).

Every stimulus event carries metadata: `BNUM, TNUM, EVNT, COUP, TTYP, TRGT, ISI (s),
REST (s), S1FQ, S2FQ, CNAM, CTYP='cnv', WAIT, DISP`.

## Analysis note

CNV epochs are **S1-locked** and must span the S1-S2 interval plus post-S2: use roughly
`[-0.2, +2.0] s` around `WARN` with a `[-200, 0] ms` pre-S1 baseline — longer than the
OmissionPilot preprocessing epoch (`[-0.1, 1.2] s`), which must be extended before reusing
that pipeline. Sort the S1-locked CNV by upcoming couple via `COUP`; S2-locked ERPs
(P3 to the `same` target, N1 to the `diff` tone) time-lock to the `S2SM`/`S2DF` events.

## Dependencies

- MATLAB with PsychToolbox-3 (for lab mode)
- EGI NetStation (optional, for EEG triggers)
