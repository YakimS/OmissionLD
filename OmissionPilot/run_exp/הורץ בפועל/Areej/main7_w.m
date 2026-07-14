close all
clear
clc
PsychDefaultSetup(2);
InitializeMatlabOpenGL;
InitializePsychSound

try
    parallelPortAddress = hex2dec('EEFC'); %LPT3
    ioObj = io64;
    status = io64(ioObj);
    io64(ioObj,parallelPortAddress,0);
catch
    parallelPortAddress=NaN;
end

% Directories
exp_folderpath = 'C:\Users\Owner\Documents\Experiments\Areej\REFLECT\Experiment_Scripts\Wakefulness\stim_2';
baseDir = 'C:\Users\Owner\Documents\Experiments\Areej\REFLECT\Data\RawData\exp_output';

% Enter participant ID and Language 
[sub, Lang] = getParticipantID();


% create folders in file explorer
participantID = upper(sub);
session = 'W';
folderPath = fullfile(baseDir, participantID, session);
mkdir(folderPath)

%create file for automatic save of exp_output
filePath = fullfile(folderPath, 'stimuli_info.csv');

% initialize variables
mode = 1; %flag for choosing trial number 
number_of_trials = 60; % around 160 for sleep session
number_of_conditions = 3; % 2 pleasant, 2 unpleasant and 1 blank
odor_value= zeros(1,number_of_trials); 

% ensures the function produces a new sequence of stimuli order
rng("shuffle")

% generate a randomized vector of odor values (1-5)
perm_ind = 1:number_of_conditions:number_of_trials;
for i = perm_ind
    odor_value(1,i:i+(number_of_conditions-1))= randperm(number_of_conditions,number_of_conditions);
end

% aromashooter serial numbers
lookup1 = {"All_unpls","All_pls","ASN3A01270"}; % add serial number of blank aromashooter

% ports numbers
port = {[1,2,3],[4,5,6],[1,2,3,4,5,6]}; % 1+2+3 Urine , 1+2+3 Jasmine, 4+5+6 Smoke, 4+5+6 Shampoo

%names of odors
O_name = {"Urine","Shampoo","blank"};

% pleasantness of odors
O_pls = {"unpls","pls","blank"};

% arrange aromashooter serial numbers that correspond to the odor value 
Olfacto_name = lookup1(odor_value); 

% arrange port numbers for activation that correspond to the odor value 
port_num = port(odor_value); 

%
odor_name = O_name(odor_value);

%
odor_pls = O_pls(odor_value);


%% connect to speakers
devices = PsychPortAudio('GetDevices');

for ii=1:length(devices)
    if strcmp(devices(ii).DeviceName, 'Headphones (Realtek(R) Audio)') && strcmp(devices(ii).HostAudioAPIName, 'MME')
        OutputDevice=devices(ii).DeviceIndex;
    end
end

% open psychotoolbox audi device and designate it as a master audio device
% that controls timing and settings for audio

pamaster = PsychPortAudio('Open', OutputDevice, 1+8, 3 ,[], 1); %master sound device
PsychPortAudio('Start', pamaster, 0, 0, 1); % start the master audio device
pahandleplay = PsychPortAudio('OpenSlave', pamaster); % inherits master device srate and clock and it is the device that actually plays the sound

% read audio files 

soundsfolderPath =[exp_folderpath '\sounds2']; %folder path for sounds %%[running_folder,'\\]
soundfiles = dir(fullfile(soundsfolderPath, '*.wav'));
for numSounds = 1:length(soundfiles) %goes through each sound file, makes it ready to be played through the psychport audio
            chosenSound{numSounds} = soundfiles(numSounds).name;
            soundPath = [soundsfolderPath '\' chosenSound{numSounds}];
            [audiodata, freq] = audioread(soundPath);
            numSamples = size(audiodata, 1);
            numChannels = size(audiodata, 2);

            % 2. Define fade duration in seconds
            fadeDuration = 1.5; 
            fadeSamples = round(fadeDuration * freq);

            % 3. Create linear fade envelope
            fadeIn = linspace(0, 1, fadeSamples)';
            fadeOut = linspace(1, 0, fadeSamples)';

            % 4. Apply fade-in to the beginning
            if numSamples > fadeSamples * 2
                audiodata(1:fadeSamples, :) = audiodata(1:fadeSamples, :) .* repmat(fadeIn, 1, numChannels);
                % 5. Apply fade-out to the end
                audiodata(end-fadeSamples+1:end, :) = audiodata(end-fadeSamples+1:end, :) .* repmat(fadeOut, 1, numChannels);
            end
            %audiodata=audiodata*soundvol(numSounds); chnage volume
            sndlength(numSounds)=length(audiodata)/freq;
            Buffer(numSounds) = PsychPortAudio('CreateBuffer', pahandleplay, mean(audiodata,2)');
end

% ensures the function produces a new sequence of stimuli order
rng("shuffle")
                                                      
% generate a randomized vector of sound values (1-5)
sound_value= zeros(1,number_of_trials);
for i = perm_ind
    sound_value(1,i:i+(number_of_conditions-1))= randperm(number_of_conditions,number_of_conditions);
end
       
% audio variables names
lookup2 = {chosenSound{1}, chosenSound{2}, "blank"};

% arrange audio variables in an order that correspond to the sound value 
sound_name = lookup2(sound_value);                                                                                                                                                                                                           

%
S_pls = {"pls","unpls","blank"};
%
sound_pls = S_pls(sound_value);

% stimuli codes for netstation 
STIM1 = ["OD01","OD02","OD03"]; % odors
STIM2 =["SO01","SO02","SO03"]; % sounds

%
odor_EEG_code = STIM1(odor_value);
sound_EEG_code = STIM2(sound_value);

% table for stimuli order
% Define size (0 rows, 3 columns)
sz = [number_of_trials 11];

% Define variable types
varTypes = ["double","string","cell","string","double","string","string","string","double","string","string"];

% Define variable names
varNames = ["Trial","serial_num","ports","odor_name","odor_value","odor_pls","odor_EEG_code","sound_name","sound_value","sound_pls","sound_EEG_code"];
% Create empty table
stimuli_info_table = table('Size', sz,'VariableTypes', varTypes, 'VariableNames', varNames);
disp(stimuli_info_table)

resp = table('Size', [3 1],'VariableTypes',["string"],'VariableNames',["Q"]);

%% Xtrodes
fprintf('============================================================\n');
fprintf('Smart Multi-Annotation Trigger System - Main Script\n');
fprintf('============================================================\n\n');

% Initialize the annotation system
try
    config = AnnotationTriggerSystem.Init('annotation_options.txt');
catch ME
    fprintf('Initialization failed: %s\n', ME.message);
    % pause;
    % config = AnnotationTriggerSystem.Init('annotation_options.txt'); %% in case disabling the allow connection window failed
    return;
end


%% Instructions
KbName('UnifyKeyNames');
keySpace = KbName('space'); 
keyEsc   = KbName('ESCAPE');

Screen('Preference', 'VisualDebugLevel', 3);
Screen('Preference', 'SkipSyncTests', 2); 
Screen('Preference', 'TextRenderer', 0);

screenNumber = 0;
black = BlackIndex(screenNumber);
 
screenRect = Screen('Rect', screenNumber);
[win, winRect] = Screen('OpenWindow', screenNumber, black, [0 0 1920 1080]); %% try 


instructions = fileread(['instruction_texts\Start_inst_' Lang '.txt']);
title = fileread(['instruction_texts\Start_inst_title_' Lang '.txt']);

% DrawFormattedText(win, instructions, 'center', 'center', 255, 120);
[screenX, screenY] = Screen('WindowSize', win);

margin = 100;

% Title
Screen('TextSize', win, 42);
Screen('TextStyle', win, 1);
if ismember(Lang,"ARB")
    DrawFormattedText(win,[ 'Reflect ' title],...
  'center', screenY * 0.15, 255, ...
    screenX - margin);

elseif ismember(Lang,"HEB")
    DrawFormattedText(win, ['Reflect ' title],'center', screenY * 0.15, 255);
else
    DrawFormattedText(win, title,'center', screenY * 0.15, 255);
end

% Body
Screen('TextSize', win, 24);
Screen('TextStyle', win, 0);


DrawFormattedText(win,instructions,...
'center', screenY * 0.25, 255, ...
screenX - margin,[],[],2);


Screen('Flip', win);

while true
    [keyIsDown, ~, keyCode] = KbCheck;
    if keyIsDown
        if keyCode(keySpace)
            KbReleaseWait;
            break;
        elseif keyCode(keyEsc) 
            KbReleaseWait;
            sca;
            return;
        end
    end
end

sca;

% start movie
videoFile = [exp_folderpath '\video\The Wonder of Whales_' Lang '.mp4'];
vlcExe = 'C:\Program Files\VideoLAN\VLC\vlc.exe';

rcHost = '127.0.0.1';
rcPort = 4212;
StartVLC(videoFile, vlcExe, rcHost, rcPort);


% open black background 
screenSize = get(0, 'ScreenSize'); % [x y width height]

f = figure( ...
    'Color', 'k', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'NumberTitle', 'off', ...
    'Name', '', ...
    'Resize', 'off', ...
    'Units', 'pixels', ...
    'Position', screenSize);

ax = axes('Parent', f, ...
          'Position', [0 0 1 1], ...
          'Color', 'k', ...
          'XColor', 'k', ...
          'YColor', 'k');

% give vlc a moment 
pause(2.5);

vlc = tcpclient(rcHost , rcPort, 'Timeout', 5);
flush(vlc);

SendVLC(vlc, 'volume 0');


% clear greeting text from vlc
pause(0.2)
if vlc.NumBytesAvailable > 0
    read(vlc, vlc.NumBytesAvailable, 'char');
end

%%
% Connection to net station
%------------------------------------------------------------------------
%-------------------------------------------------------------------------

NTPhost = '10.10.10.42'; %net station host IP address
NTPport = 55513;
NetStation('Connect', NTPhost,NTPport);

disp('Connected to Net Station.');

% % Synchronize clocks (the value '10' is an example time in milliseconds)
NetStation('Synchronize', 10); 
WaitSecs(0.5); % Wait briefly after synchronization

NetStation('GetNTPSynchronize', '10.10.10.51');
NetStation('StartRecording');


% in case of code crash mid experiment, we can restart the experiment by
% running this section and choosing the trial from which we want to
% continue
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

if mode==0
    Trial_number = getTrial();
else
    Trial_number = 1;
    mode = 0;
end




%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% start experiment
breakIndex = 1;
% start of session  
disp ("Start sesison")
 
for iTrial = Trial_number:number_of_trials
        
        
        % activate aromashooter
        %-----------------------------------------------------------------
        %-----------------------------------------------------------------
        if odor_value(iTrial) == 3
            [ShooterTiming3(iTrial)] = aromaShooter4(Olfacto_name{iTrial}, port_num{iTrial});

        else 
            [ShooterTiming3(iTrial)] = aromaShooter3(port_num{iTrial});

        end
        % xtrodes_trigger
        try
            success = AnnotationTriggerSystem.SendTrigger(config, odor_pls{iTrial});
            if ~success
                fprintf('Failed to send trigger for ''%s''\n', odor_pls{iTrial});
            end
        catch ME
            fprintf('Error sending trigger: %s\n', ME.message);
        end


        % send trigger to netstation
        NetStation('Event',STIM1(odor_value(iTrial)), ShooterTiming3(iTrial)+0.006, 0.001, num2str(iTrial), odor_value(iTrial)); 

        % 
        disp(['Sent trigger for trial: ' num2str(iTrial)]);

        pause(10)       
        aromaShooter_disconnect();

        % send end of trigger to xtrodes
        try
            success = AnnotationTriggerSystem.SendTrigger(config, "end_smell");
            if ~success
                fprintf('Failed to send trigger for ''%s''\n', "end_smell");
            end
        catch ME
            fprintf('Error sending trigger: %s\n', ME.message);
        end
        
        % inter-stimulus interval
        % WaitSecs(18+rand*6-3)
        WaitSecs(18+randi([1,4]))

      
         
        
        % play audio
        %-----------------------------------------------------------------
        %-----------------------------------------------------------------
        if(sound_name{iTrial}~="blank")

            PsychPortAudio ('Stop', pahandleplay);
            PsychPortAudio ('UseSchedule',pahandleplay,2); %<>% it then resets the presentation schedule
            PsychPortAudio ('AddToSchedule',pahandleplay,Buffer(sound_value(iTrial))); %<>% it then adds the sound to the schedule
            soundstart= GetSecs() + 0.23; %0.73
            PsychPortAudio('Start',pahandleplay,1,soundstart); %<>%it finally start playing the schedule or playlist add(planned_soundstart)
        else
            soundstart= GetSecs() + 0.225;%0.725
        end

        
        % send trigger to netstation
        NetStation('Event',STIM2(sound_value(iTrial)), soundstart+0.0713, 0.001, num2str(iTrial), sound_value(iTrial)); 

        % send trigger to xtrodes
        try
            success = AnnotationTriggerSystem.SendTrigger(config, sound_pls{iTrial});
            if ~success
                fprintf('Failed to send trigger for ''%s''\n', sound_pls{iTrial});
            end
        catch ME
            fprintf('Error sending trigger: %s\n', ME.message);
        end
        
        
        disp(['Sent trigger for trial: ' num2str(iTrial)]);
        
        pause(5)   
        % send end of trigger to xtrodes
         try
            success = AnnotationTriggerSystem.SendTrigger(config, "end_sound");
            if ~success
                fprintf('Failed to send trigger for ''%s''\n', "end_sound");
            end
        catch ME
            fprintf('Error sending trigger: %s\n', ME.message);
        end
    
        % saving variables
        %-----------------------------------------------------------------
        %-----------------------------------------------------------------
        stimuli_info_table.Trial(iTrial) = iTrial;
        stimuli_info_table.serial_num(iTrial) = Olfacto_name{iTrial};
        stimuli_info_table.ports(iTrial)= port_num(iTrial);
        stimuli_info_table.odor_name(iTrial) = odor_name{iTrial};
        stimuli_info_table.odor_value(iTrial) = odor_value(iTrial);
        stimuli_info_table.odor_pls(iTrial) = odor_pls{iTrial};
        stimuli_info_table.odor_EEG_code(iTrial) = odor_EEG_code{iTrial};
        stimuli_info_table.sound_name(iTrial) = sound_name{iTrial};
        stimuli_info_table.sound_value(iTrial) = sound_value(iTrial);
        stimuli_info_table.sound_pls(iTrial) = sound_pls{iTrial};
        stimuli_info_table.sound_EEG_code(iTrial) = sound_EEG_code{iTrial};
        

        % automatic data saving in csv
        if iTrial == 1
            % first write → creates file + writes header
            writetable(stimuli_info_table(iTrial,:), filePath);
        else
            writetable(stimuli_info_table(iTrial, :), filePath, 'WriteMode', 'append');
        end
        
        % inter-trial interval
        % WaitSecs(18+rand*6-3)
        WaitSecs(18+randi([1,4]))

        % break or experiment termination
        %-----------------------------------------------------------------
        %------------------------------------------------------------------
        if iTrial == 60
            resp.Q(breakIndex) = DoBreak2(vlc, breakIndex, black,Lang);
            StopVLC(vlc, keyEsc,black,Lang)
            close(f)
        elseif iTrial == 20 || iTrial ==40

            resp.Q(breakIndex) = DoBreak(vlc, breakIndex, keySpace, keyEsc, black,Lang);
            breakIndex = breakIndex +1;
        end

   
end



% file save 
% -------------------------------------------------------------------------
%--------------------------------------------------------------------------
disp('End of Experiment')
save([folderPath '\stimuli_info.mat'], 'stimuli_info_table');
save([folderPath '\response_to_questions.mat'], 'resp');
disp("Events were saved successfully")



% Stop the recording (if started via MATLAB)
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
NetStation('StopRecording'); 

% Disconnect from Net Station
NetStation('Disconnect');
disp('Disconnected from Net Station.');