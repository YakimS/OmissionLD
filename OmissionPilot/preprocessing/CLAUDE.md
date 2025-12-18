# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

EEGLAB-based preprocessing pipeline for OmissionPilot EEG data. Processes EGI .mff files and outputs EEGLAB .set files.

## Running Preprocessing

```matlab
cd OmissionPilot/preprocessing
preprocess_OmissionPilot  % Processes all .mff files in RAW_DATA_DIR
```

## Data Paths

- EEGLAB: `D:\matlab_libs\eeglab2024.2`
- Raw MFF files: `D:\omissionPilot\raw_data`
- Output: `D:\omissionPilot\preprocessed`

## Pipeline Steps

1. Load MFF file (`pop_mffimport`)
2. High-pass filter (default: 0.1 Hz)
3. Low-pass filter (default: 40 Hz)
4. Re-reference to average
5. Epoch around events [-0.2, 0.8] sec
6. Baseline correction [-200, 0] ms
7. Artifact rejection (±100 µV threshold)
8. Save as .set file

## Key Parameters (preprocess_OmissionPilot.m)

```matlab
% Filtering (lines 18-19)
HIGHPASS_FREQ = 0.1;
LOWPASS_FREQ = 40;

% Epoching (lines 22-23)
EPOCH_START = -0.2;
EPOCH_END = 0.8;

% Baseline (lines 26-27)
BASELINE_START = -200;  % ms
BASELINE_END = 0;

% Events to epoch (line 30)
EVENTS_OF_INTEREST = {'STND', 'TON1', 'TON2', 'DEVT', 'OMIS'};

% Artifact threshold (line 33)
ARTIFACT_THRESHOLD = 100;  % µV
```

## Event Types from Experiment

- `STND`: Standard tone (single condition)
- `TON1`: First tone in double condition
- `TON2`: Second tone in double condition
- `DEVT`: Deviant tone
- `OMIS`: Omission (no tone)
- `BGIN`: Block begin (not epoched)
- `STRT`/`SEND`: Experiment start/end (not epoched)

## Dependencies

- MATLAB
- EEGLAB 2024.2 (with MFFMatlabIO plugin for .mff import)
