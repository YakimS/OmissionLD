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

- EEGLAB: `D:\matlab_libs\eeglab2025.1.0`
- Raw MFF files: `D:\omissionPilot\raw_data\gel` (128-channel) and `D:\omissionPilot\raw_data\saline` (256-channel)
- Output: `D:\omissionPilot\preprocessed\gel` and `D:\omissionPilot\preprocessed\saline`

## Data Type Selection

All scripts support a `DATA_TYPE` flag: `'gel'`, `'saline'`, or `'both'`

Electrode sets per data type:
- **Gel (128-channel)**:
  - Central: E106, E7, E80, E55, E31
  - Frontal: E4, E5, E10, E11, E12, E16, E18, E19
- **Saline (256-channel)**:
  - Central: E9, E186, E45, E81, E132, E257
  - Frontal: E6, E7, E14, E15, E16, E22, E23

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
