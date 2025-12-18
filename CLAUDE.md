# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

MATLAB implementation of an auditory omission oddball paradigm based on Hughes et al. (2001) "Responses of Human Auditory Association Cortex to the Omission of an Expected Acoustic Event". Used for EEG studies with PsychToolbox audio and optional EGI NetStation integration.

## Project Structure

All active code is in `Hughes/`:
- `config_Hughes_exp.m` - Configuration generator (Phase 1)
- `run_Hughes_exp.m` - Experiment runner (Phase 2)
- `Config-*.mat` - Pre-generated experiment configurations

## Running the Experiment

**Phase 1 - Generate Configuration:**
```matlab
cd Hughes
config_Hughes_exp  % Generates .mat file with trial sequences and stimuli
```

**Phase 2 - Run Experiment:**
```matlab
cd Hughes
run_Hughes_exp  % Loads config and runs with audio playback
```

## Key Configuration Parameters

In `config_Hughes_exp.m` (lines 18-38):
- `OUTPUT_FILE` - Config filename or `false` for auto datetime naming
- `WITHIN_PAIR_ISI` - Array of ISI values; creates separate blocks for each
- `N_TRIALS`, `OMISSION_FREQ` - Trial count and omission proportion
- `STIM_DURATION`, `STIM_FREQ` - Tone parameters (default: 100ms, 1000Hz)

In `run_Hughes_exp.m` (lines 9, 17-19):
- `CONFIG_FILE` - Path to .mat from Phase 1
- `USE_NETSTATION` - Enable EEG triggers
- `NETSTATION_HOST`, `NETSTATION_PORT` - NetStation connection

## Experimental Design

Two procedures are randomly interleaved across blocks:

**Procedure 1 (Single Tones):** Individual tones with random omissions
- Trial: [tone or silence] → ISI

**Procedure 2 (Tone Pairs):** Paired tones where second may be omitted
- Trial: [tone1] → within_pair_ISI → [tone2 or silence] → remaining_ISI

Trial constraints: Initial N trials are always standards; minimum spacing enforced between omissions.

## Data Flow

```
config_Hughes_exp.m → exp_meta struct → .mat file → run_Hughes_exp.m → audio + triggers
```

The `exp_meta` struct contains: stimuli (`beep_sound`), trial sequences, ISI arrays, randomized block order, and all parameters.

## NetStation EEG Triggers

- `STRT`/`SEND` - Experiment start/end
- `BGIN` - Block begin (includes block#, procedure, within-pair ISI)
- `BEEP` - Tone onset (sent with precise PsychPortAudio timestamp)

All triggers include metadata: `BNUM` (block), `TNUM` (trial), `WISI` (within-pair ISI in ms).

## Dependencies

- MATLAB with PsychToolbox-3
- EGI NetStation (optional, for EEG)
