function generate_tones()
% GENERATE_TONES Generate pure tone WAV files for oddball experiment
%
% Generates three pure tones:
%   - Low:    650 Hz
%   - Medium: 1428 Hz
%   - High:   3137 Hz
%
% All tones are:
%   - Duration: 50 ms
%   - Amplitude: Normalized to [-1, 1] range
%   - Sample rate: 44100 Hz
%   - NOTE: Amplitude should be adjusted manually per subject's preference
%
% Output files:
%   - tone_low_650Hz.wav
%   - tone_medium_1428Hz.wav
%   - tone_high_3137Hz.wav

    % Define parameters
    sample_rate = 44100;  % Hz
    duration = 0.050;      % 50 ms

    % Frequency definitions
    frequencies = [650, 1428, 3137];  % Hz
    tone_names = {'low_650Hz', 'medium_1428Hz', 'high_3137Hz'};

    % Time vector
    t = 0:1/sample_rate:duration;

    % Apply 5ms cosine-squared ramp to avoid clicks
    ramp_duration = 0.005;  % 5 ms
    ramp_samples = round(ramp_duration * sample_rate);

    % Create output directory if it doesn't exist
    output_dir = 'tones';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    fprintf('Generating pure tones...\n');
    fprintf('=====================================\n');

    % Generate each tone
    for i = 1:length(frequencies)
        freq = frequencies(i);

        % Generate pure tone
        tone = sin(2 * pi * freq * t);

        % Apply onset/offset ramps (cosine-squared envelope)
        envelope = ones(size(tone));

        % Onset ramp
        for j = 1:ramp_samples
            envelope(j) = sin(pi * j / (2 * ramp_samples))^2;
        end

        % Offset ramp
        for j = 1:ramp_samples
            envelope(end - j + 1) = sin(pi * j / (2 * ramp_samples))^2;
        end

        % Apply envelope
        tone = tone .* envelope;

        % Normalize to [-1, 1] range to prevent clipping
        % Amplitude will be adjusted manually per subject's preference during experiment
        tone = tone / max(abs(tone));

        % Save to WAV file
        filename = fullfile(output_dir, sprintf('tone_%s.wav', tone_names{i}));
        audiowrite(filename, tone, sample_rate);

        fprintf('Generated: %s (%d Hz, 50 ms)\n', filename, freq);
    end

    fprintf('=====================================\n');
    fprintf('All tones generated successfully!\n');
    fprintf('Output directory: %s\n', output_dir);
end
