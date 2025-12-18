# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is an auditory omission oddball experiment implementation based on Hughes et al. (2001) "Responses of Human Auditory Association Cortex to the Omission of an Expected Acoustic Event". The experiment uses PsychToolbox for precise audio presentation and optionally integrates with EGI NetStation for EEG trigger delivery.

## Two-Phase Workflow

### Phase 1: Configuration Generation (`config_Hughes_exp.m`)

This script generates experiment parameters and saves them to a `.mat` file. It must be run before the experiment.

**Key parameters to configure (lines 18-38):**
- `OUTPUT_FILE`: Set filename (string) or `false` for auto datetime naming
- `WITHIN_PAIR_ISI`: Array of ISI values (e.g., `[0.1, 0.3, 0.5]`) - creates separate blocks for each
- `N_TRIALS`: Total trials per block
- `OMISSION_FREQ`: Proportion of omission trials (e.g., 0.1 = 10%)
- `STIM_DURATION`, `STIM_FREQ`: Tone parameters
- `ISI_MIN`, `ISI_MAX`: Range for randomized inter-stimulus intervals

**Output:**
- Displays complete experiment preview including:
  - Randomized block order (Procedure 1 and 2 for each ISI condition)
  - Duration estimates per block and total
  - Trial sequence statistics
- Saves `exp_meta` struct to `.mat` file containing all stimuli, sequences, and parameters

**Run:** Open MATLAB, run `config_Hughes_exp` from Hughes directory

### Phase 2: Experiment Execution (`run_Hughes_exp.m`)

Loads configuration and runs the actual experiment with audio playback.

**Configuration (lines 10-21):**
- `CONFIG_FILE`: Path to `.mat` file from Phase 1
- `USE_NETSTATION`: Enable/disable EEG triggers
- `NETSTATION_HOST`, `NETSTATION_PORT`: NetStation connection settings

**Run:** Open MATLAB, run `run_Hughes_exp` from Hughes directory

## Experimental Design

The experiment implements two procedures that are randomly interleaved:

**Procedure 1: Single Tones**
- Plays single tones with random omissions
- Each trial: [tone or silence] + ISI

**Procedure 2: Tone Pairs**
- Plays tone pairs where the second tone is sometimes omitted
- Each trial: [tone1] + within_pair_ISI + [tone2 or silence] + remaining_ISI

**Block structure:**
- For each ISI value in `WITHIN_PAIR_ISI`, both Procedure 1 and Procedure 2 are run
- Block order is randomized
- Each block has independent randomized trial sequences respecting constraints

**Trial constraints:**
- First N trials are always standards (no omissions)
- Minimum spacing between omissions
- Randomized ISI durations (uniform distribution between ISI_MIN and ISI_MAX)

## NetStation EEG Triggers

When `USE_NETSTATION = true`, the following triggers are sent:

- `STRT`: Experiment start
- `BGIN`: Block begin (includes block#, procedure type, within-pair ISI)
- `TONE`: Single tone in Procedure 1
- `TON1`: First tone in Procedure 2
- `TON2`: Second tone in Procedure 2
- `OMIS`: Omission event
- `SEND`: Experiment end

All triggers include metadata: `BNUM` (block number), `TNUM` (trial number), `PROC` (procedure type).

## Code Architecture

**Two-file structure:**
1. `config_Hughes_exp.m` - Configuration generator (script)
   - Contains `create_omission_sequence()` helper function for constraint-based trial generation
   - Generates audio stimulus with onset/offset ramps
   - Creates randomized conditions and trial sequences
   - Calculates accurate duration estimates

2. `run_Hughes_exp.m` - Experiment runner (script)
   - Contains `play_beep()` helper function wrapping PsychPortAudio calls
   - Loads exp_meta from .mat file
   - Handles NetStation connection and triggers
   - Executes trial loops for both procedures

**Data flow:**
```
config_Hughes_exp.m
  → generates exp_meta struct
  → saves to .mat file
  → run_Hughes_exp.m loads
  → executes experiment
```

## Duration Calculation Logic

Both procedures take identical time per trial despite different structure:
- **Procedure 1**: `STIM_DURATION + ISI`
- **Procedure 2**: Also `STIM_DURATION + ISI` (within-pair timing is absorbed into total ISI)

The config script correctly calculates this (see line 102 in config_Hughes_exp.m).

## Dependencies

- MATLAB (tested environment)
- PsychToolbox-3 (for audio playback via PsychPortAudio)
- EGI NetStation (optional, for EEG integration)

## File Types

- `.m` files: MATLAB scripts (configuration and experiment runner)
- `.mat` files: MATLAB binary data containing experiment configurations (exp_meta struct)
