% run_Hughes_exp.m
% Omission oddball experiment based on Hughes et al. (2001)



clear all;
close all;
clc;

isLabComputer = false;
if isLabComputer Priority(1); end

CONFIG_FILE = 'Config-test_timingTest.mat'; % "Config-Hughes_params.mat" 'Config-test.mat'


%% Load experiment configuration
loaded = load(CONFIG_FILE);
exp_meta = loaded.exp_meta;

% NetStation Configuration
USE_NETSTATION = isLabComputer & true;  % Set to true to enable NetStation triggers
NETSTATION_HOST = '10.10.10.42';  % NetStation host IP
NETSTATION_PORT = 55513;             % NetStation port


%% Connect to NetStation
if isLabComputer & USE_NETSTATION
    fprintf('Connecting to NetStation at %s:%d...\n', NETSTATION_HOST, NETSTATION_PORT);
    [nsstatus, nserror] = NetStation('Connect', NETSTATION_HOST, NETSTATION_PORT);
    if nsstatus ~= 0
        error('Could not connect to NetStation: %s', nserror);
    end
    NetStation('GetNTPSynchronize','10.10.10.51');
    WaitSecs(1);
    fprintf('NetStation connected and NTPSynchronized.\n');
end

%% sound stuff
pahandleplay = [];  % Initialize for non-lab computer case

if isLabComputer
    InitializePsychSound(1);

    devices = PsychPortAudio('GetDevices');
    for ii=1:length(devices)
        if strcmp(devices(ii).DeviceName, 'Headphones (Realtek(R) Audio)') && strcmp(devices(ii).HostAudioAPIName, 'MME')
            OutputDevice=devices(ii).DeviceIndex;
        end
    end

    pamaster = PsychPortAudio('Open', OutputDevice, 1+8, 3, exp_meta.SAMPLE_RATE, 1); %mastimingtester sound device
    PsychPortAudio('Start', pamaster, 0, 0, 1);
    pahandleplay = PsychPortAudio('OpenSlave', pamaster);
    PsychPortAudio('Volume', pahandleplay,1);
end

%% Run Experiment
fprintf('\n===========================================\n');
fprintf('   OMISSION ODDBALL EXPERIMENT\n');
fprintf('===========================================\n\n');

NetStation('StartRecording');

% Send experiment start trigger
if USE_NETSTATION
    NetStation('Event', 'STRT', GetSecs, 0.001, 'expt', 1);
    WaitSecs(2);
end

for block = 1:exp_meta.n_total_blocks
    current_proc = exp_meta.conditions(block, 1);
    current_isi = exp_meta.conditions(block, 2);
    trial_seq = exp_meta.trial_sequences{block};
    isi_array = exp_meta.isi_arrays{block};

    % Determine procedure name
    if current_proc == 1
        proc_name = 'Procedure 1 (Single Tones)';
    else
        proc_name = 'Procedure 2 (Tone Pairs)';
    end

    fprintf('===========================================\n');
    fprintf('Block %d/%d: %s\n', block, exp_meta.n_total_blocks, proc_name);
    fprintf('Within-pair ISI: %.3f s\n', current_isi);
    fprintf('Trials: %d (%.0f%% omissions)\n', exp_meta.N_TRIALS, exp_meta.OMISSION_FREQ*100);

    % Send block start trigger
    if USE_NETSTATION
        NetStation('GetNTPSynchronize','10.10.10.51');
        NetStation('Event', 'BGIN', GetSecs, 0.001, 'BNUM', block, 'PROC', current_proc, 'WISI', current_isi*1000);
    end

    if exp_meta.WAIT_BETWEEN_BLOCKS
        WaitSecs(exp_meta.WAIT_BETWEEN_BLOCKS);
    end
    fprintf('===========================================\n');

    if current_proc == 1
        % Procedure 1: Single tones with omissions
        fprintf('Starting Procedure 1...\n');
        for trial = 1:exp_meta.N_TRIALS
            if exp_meta.SHOW_TRIAL_OUTPUT; fprintf('Trial %d/%d: ', trial, exp_meta.N_TRIALS); end

            if trial_seq(trial) == 1
                % Standard trial: play single tone
                if exp_meta.SHOW_TRIAL_OUTPUT; fprintf('TONE\n'); end
                startTime = play_beep(pahandleplay, exp_meta.beep_sound, isLabComputer, exp_meta.SAMPLE_RATE);
                NetStation('Event', 'BEEP', startTime, 0.001, 'BNUM', block, 'TNUM', trial, 'WISI', current_isi*1000);
                if USE_NETSTATION
                    % NetStation('Event', 'P1ST', GetSecs, 0.001, 'BNUM', block, 'TNUM', trial, 'WISI', current_isi*1000);
                end
                WaitSecs(isi_array(trial));
            else
                % Omission trial: silence
                if exp_meta.SHOW_TRIAL_OUTPUT; fprintf('OMISSION\n'); end
                if USE_NETSTATION
                    % NetStation('Event', 'P1DV', GetSecs, 0.001, 'BNUM', block, 'TNUM', trial, 'WISI', current_isi*1000);
                end
                WaitSecs(exp_meta.STIM_DURATION + isi_array(trial));
            end
        end

    else
        % Procedure 2: Tone pairs with second tone omissions
        fprintf('Starting Procedure 2...\n');
        for trial = 1:exp_meta.N_TRIALS
            if exp_meta.SHOW_TRIAL_OUTPUT; fprintf('Trial %d/%d: ', trial, exp_meta.N_TRIALS); end

            % First tone always plays
            startTime = play_beep(pahandleplay, exp_meta.beep_sound, isLabComputer, exp_meta.SAMPLE_RATE);
            NetStation('Event', 'BEEP', startTime, 0.001, 'BNUM', block, 'TNUM', trial, 'WISI', current_isi*1000);

            % Send trigger at trial start indicating standard or deviant
            if USE_NETSTATION
                if trial_seq(trial) == 1
                    % NetStation('Event', 'P2ST', GetSecs, 0.001, 'BNUM', block, 'TNUM', trial, 'WISI', current_isi*1000);
                else
                    % NetStation('Event', 'P2DV', GetSecs, 0.001, 'BNUM', block, 'TNUM', trial, 'WISI', current_isi*1000);
                end
            end

            WaitSecs(current_isi);

            if trial_seq(trial) == 1
                % Standard trial: play second tone
                if exp_meta.SHOW_TRIAL_OUTPUT; fprintf('TONE PAIR\n'); end
                startTime = play_beep(pahandleplay, exp_meta.beep_sound, isLabComputer, exp_meta.SAMPLE_RATE);
                NetStation('Event', 'BEEP', startTime, 0.001, 'BNUM', block, 'TNUM', trial, 'WISI', current_isi*1000);
                WaitSecs(isi_array(trial) - current_isi - exp_meta.STIM_DURATION);
            else
                % Omission trial: omit second tone
                if exp_meta.SHOW_TRIAL_OUTPUT; fprintf('TONE-OMISSION\n'); end
                WaitSecs(exp_meta.STIM_DURATION);  % Silent period where second tone would be
                WaitSecs(isi_array(trial) - current_isi - exp_meta.STIM_DURATION);
            end
        end
    end

    fprintf('\nBlock %d/%d complete.\n\n', block, exp_meta.n_total_blocks);

    if block < exp_meta.n_total_blocks && exp_meta.WAIT_BETWEEN_BLOCKS
        WaitSecs(exp_meta.WAIT_BETWEEN_BLOCKS);
    end
end
fprintf('EXPERIMENT COMPLETE\n');

% Send experiment end trigger
if USE_NETSTATION
    WaitSecs(2);   NetStation('Event', 'SEND', GetSecs, 0.001, 'end', 1);
end

if isLabComputer
    PsychPortAudio('Close', pahandleplay);
end


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