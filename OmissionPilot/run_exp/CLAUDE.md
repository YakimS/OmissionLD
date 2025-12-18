# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

OmissionPilot experiment: auditory paradigm with single and double tone conditions, featuring omission and deviant tone deviants. Used for EEG studies with PsychToolbox audio and optional EGI NetStation integration.

## Running the Experiment

**Phase 1 - Generate Configuration:**
```matlab
cd OmissionPilot/run_exp
config_OmissionPilot  % Generates .mat file with trial sequences and stimuli
```

**Phase 2 - Run Experiment:**
```matlab
cd OmissionPilot/run_exp
run_OmissionPilot  % Loads config and runs with audio playback
```

Set `isLabComputer = true` in run_OmissionPilot.m for PsychToolbox + NetStation; `false` uses MATLAB's built-in `sound()` for testing.

## Experiment Structure

**Conditions** (per ITI × condition type):
- For each ITI value:
  - 1 Single condition: single tone at t=0
  - N Double conditions: tone pairs with varying ISI

**Deviants:**
- Omission: tone is not played
- Deviant tone: different frequency played
- In single conditions: the single tone is manipulated
- In double conditions: only the second tone is manipulated

**Block Structure:**
- Each condition has BLOCKS_PER_CONDITION blocks
- Blocks are shuffled randomly across all conditions
- Per block, which frequency is "standard" vs "deviant" is randomly assigned

## Key Parameters (config_OmissionPilot.m)

```matlab
% Frequencies (lines 26-27)
exp_meta.FREQ_1 = 1000;           % First frequency
exp_meta.FREQ_2 = 1200;           % Second frequency

% Timing (lines 30-34)
exp_meta.ITI_VALUES = [1.2];      % Inter-trial intervals
exp_meta.ISI_VALUES = [0.1, 0.25]; % ISI for double conditions
exp_meta.ITI_JITTER_MIN/MAX       % Jitter range applied to ITI

% Trials (lines 37-42)
exp_meta.N_TRIALS = 600;          % Trials per condition
exp_meta.BLOCKS_PER_CONDITION = 6;
exp_meta.OMISSION_RATE = 0.1;
exp_meta.DEVIANT_RATE = 0.1;
exp_meta.N_STANDARD_BEGIN = 10;
exp_meta.MIN_STANDARD_BETWEEN = 3; % Between ANY deviants
```

## Data Flow

```
config_OmissionPilot.m → exp_meta struct → .mat file → run_OmissionPilot.m → audio + triggers
```

The `exp_meta` struct contains: stimuli (`beep_freq1`, `beep_freq2`), events table, block info, frequency assignments, and all parameters.

## Events Table

The config generates a complete events table with columns:
- `block`, `trial`, `event_in_trial`
- `condition_type` ('single'/'double'), `condition_name`
- `iti`, `isi`
- `trial_type` ('standard'/'omission'/'deviant')
- `event_type` ('tone'/'omission')
- `freq` (tone frequency, NaN for omission)
- `wait_after` (seconds to wait after this event)
- `trigger` (NetStation trigger code)
- `display_str` (console output string)

Double condition trials have 2 rows (one per tone).

## NetStation Triggers

- `STRT`/`SEND`: Experiment start/end
- `BGIN`: Block begin (includes BNUM, COND, PRMS)
- `STND`: Standard tone (single condition)
- `TON1`: First tone in double condition
- `TON2`: Second tone in double condition (standard)
- `DEVT`: Deviant tone
- `OMIS`: Omission

All event triggers include full metadata: BNUM, TNUM, EVNT, CTYP, CNAM, ITI, ISI, TTYP, ETYP, FREQ, WAIT, DISP.

## Dependencies

- MATLAB with PsychToolbox-3 (for lab mode)
- EGI NetStation (optional, for EEG triggers)
