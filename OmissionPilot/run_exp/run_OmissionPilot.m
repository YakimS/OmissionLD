% run_OmissionPilot.m
% Runs OmissionPilot experiment from pre-generated configuration
%
% Usage:
%   1. First run config_OmissionPilot.m to generate configuration
%   2. Set CONFIG_FILE below to the generated .mat file
%   3. Run this script

clear all;
close all;
clc;

isLabComputer = true;
if isLabComputer; Priority(1); end

CONFIG_FILE = 'Config-DoubleOmission_iti-1200.mat';

%% Load experiment configuration
loaded = load(CONFIG_FILE);
exp_meta = loaded.exp_meta;

% NetStation Configuration
USE_NETSTATION = isLabComputer & true;
NETSTATION_HOST = '10.10.10.42';
NETSTATION_PORT = 55513;

%% Connect to NetStation
if isLabComputer && USE_NETSTATION
    fprintf('Connecting to NetStation at %s:%d...\n', NETSTATION_HOST, NETSTATION_PORT);
    [nsstatus, nserror] = NetStation('Connect', NETSTATION_HOST, NETSTATION_PORT);
    if nsstatus ~= 0
        error('Could not connect to NetStation: %s', nserror);
    end
    NetStation('GetNTPSynchronize', '10.10.10.51');
    WaitSecs(1);
    fprintf('NetStation connected and NTPSynchronized.\n');
end

%% Sound setup
pahandleplay = [];

if isLabComputer
    InitializePsychSound(1);

    devices = PsychPortAudio('GetDevices');
    for ii = 1:length(devices)
        if strcmp(devices(ii).DeviceName, 'Headphones (Realtek(R) Audio)') && strcmp(devices(ii).HostAudioAPIName, 'MME')
            OutputDevice = devices(ii).DeviceIndex;
        end
    end

    pamaster = PsychPortAudio('Open', OutputDevice, 1+8, 3, exp_meta.SAMPLE_RATE, 1);
    PsychPortAudio('Start', pamaster, 0, 0, 1);
    pahandleplay = PsychPortAudio('OpenSlave', pamaster);
    PsychPortAudio('Volume', pahandleplay, 1);
end

%% Precompute audio buffers for both frequencies
% Actually just map both frequencies to their sounds
beep_sounds = containers.Map('KeyType', 'double', 'ValueType', 'any');
beep_sounds(exp_meta.FREQ_1) = exp_meta.beep_freq1;
beep_sounds(exp_meta.FREQ_2) = exp_meta.beep_freq2;

%% Run Experiment
fprintf('\n===========================================\n');
fprintf('   OMISSION PILOT EXPERIMENT\n');
fprintf('===========================================\n\n');

if USE_NETSTATION
    NetStation('StartRecording');
end

% Send experiment start trigger
if USE_NETSTATION
    NetStation('Event', 'STRT', GetSecs, 0.001, 'expt', 1);
    WaitSecs(2);
end

events = exp_meta.events_table;
current_block = 0;

for i = 1:height(events)
    event = events(i, :);

    % Check if new block
    if event.block ~= current_block
        % Block transition
        if current_block > 0
            fprintf('\nBlock %d complete.\n', current_block);
            if current_block < exp_meta.n_total_blocks && exp_meta.WAIT_BETWEEN_BLOCKS > 0
                do_wait(exp_meta.WAIT_BETWEEN_BLOCKS, isLabComputer);
            end
        end

        current_block = event.block;
        cond_name = event.condition_name{1};

        fprintf('\n===========================================\n');
        fprintf('Block %d/%d: %s\n', current_block, exp_meta.n_total_blocks, cond_name);
        fprintf('===========================================\n');

        % Send block start trigger
        if USE_NETSTATION
            NetStation('GetNTPSynchronize', '10.10.10.51');
            NetStation('Event', 'BGIN', GetSecs, 0.001, ...
                       'BNUM', current_block, ...
                       'COND', cond_name, ...
                       'PRMS', exp_meta.params_string);
        end

        if exp_meta.WAIT_BETWEEN_BLOCKS > 0
            do_wait(exp_meta.WAIT_BETWEEN_BLOCKS, isLabComputer);
        end
    end

    % Display trial info (only on first event of trial)
    if exp_meta.SHOW_TRIAL_OUTPUT && event.event_in_trial == 1
        fprintf('Trial %d: ', event.trial);
    end

    % Execute event
    if strcmp(event.event_type{1}, 'tone')
        % Play tone
        freq = event.freq;
        sound_data = beep_sounds(freq);

        if exp_meta.SHOW_TRIAL_OUTPUT
            fprintf('%s ', event.display_str{1});
        end

        startTime = play_beep(pahandleplay, sound_data, isLabComputer, exp_meta.SAMPLE_RATE);

        if USE_NETSTATION
            NetStation('Event', event.trigger{1}, startTime, 0.001, ...
                       'BNUM', event.block, ...
                       'TNUM', event.trial, ...
                       'EVNT', event.event_in_trial, ...
                       'CTYP', event.condition_type{1}, ...
                       'CNAM', event.condition_name{1}, ...
                       'ITI', event.iti, ...
                       'ISI', event.isi, ...
                       'TTYP', event.trial_type{1}, ...
                       'ETYP', event.event_type{1}, ...
                       'FREQ', freq, ...
                       'WAIT', event.wait_after, ...
                       'DISP', event.display_str{1});
        end

    else
        % Omission - no sound, just trigger
        if exp_meta.SHOW_TRIAL_OUTPUT
            fprintf('%s ', event.display_str{1});
        end

        if isLabComputer
            omis_time = GetSecs;
        else
            omis_time = posixtime(datetime('now'));
        end

        if USE_NETSTATION
            NetStation('Event', event.trigger{1}, omis_time, 0.001, ...
                       'BNUM', event.block, ...
                       'TNUM', event.trial, ...
                       'EVNT', event.event_in_trial, ...
                       'CTYP', event.condition_type{1}, ...
                       'CNAM', event.condition_name{1}, ...
                       'ITI', event.iti, ...
                       'ISI', event.isi, ...
                       'TTYP', event.trial_type{1}, ...
                       'ETYP', event.event_type{1}, ...
                       'FREQ', event.freq, ...
                       'WAIT', event.wait_after, ...
                       'DISP', event.display_str{1});
        end
    end

    % Wait after event
    do_wait(event.wait_after, isLabComputer);

    % Newline after last event in trial
    if exp_meta.SHOW_TRIAL_OUTPUT
        % Check if this is the last event in this trial
        if i == height(events) || events.trial(i+1) ~= event.trial || events.block(i+1) ~= event.block
            fprintf('\n');
        end
    end
end

fprintf('\nBlock %d complete.\n', current_block);
fprintf('\n===========================================\n');
fprintf('EXPERIMENT COMPLETE\n');
fprintf('===========================================\n');

% Send experiment end trigger
if USE_NETSTATION
    WaitSecs(2);
    NetStation('Event', 'SEND', GetSecs, 0.001, 'end', 1);
    NetStation('StopRecording');
end

% Cleanup
if isLabComputer
    PsychPortAudio('Close', pahandleplay);
end

%% Helper functions
function startTime = play_beep(pahandle, beep_sound, isLabComputer, sampleRate)
    if isLabComputer
        PsychPortAudio('FillBuffer', pahandle, beep_sound);
        startTime = PsychPortAudio('Start', pahandle, 1, 0, 1);
        PsychPortAudio('Stop', pahandle, 1);
    else
        startTime = posixtime(datetime('now'));
        sound(beep_sound, sampleRate);
    end
end

function do_wait(duration, isLabComputer)
    if isLabComputer
        WaitSecs(duration);
    else
        pause(duration);
    end
end
