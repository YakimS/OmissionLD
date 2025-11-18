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
%   - Intensity: 75 dB SPL (calibrated to reference)
%   - Sample rate: 44100 Hz
%
% Output files:
%   - tone_low_650Hz.wav
%   - tone_medium_1428Hz.wav
%   - tone_high_3137Hz.wav

    % Define parameters
    sample_rate = 44100;  % Hz
    duration = 0.050;      % 50 ms
    target_db = 75;        % dB SPL

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

        % Normalize to target dB SPL
        % Reference: 0 dB SPL corresponds to amplitude of 1
        % 75 dB SPL = 20*log10(amplitude/reference)
        % amplitude = 10^(dB/20) * reference
        reference_amplitude = 0.00002;  % Standard reference for sound pressure
        target_amplitude = 10^(target_db/20) * reference_amplitude;

        % Scale to appropriate range for WAV file (normalized to prevent clipping)
        % For 75 dB, we'll use a normalized approach
        tone = tone / max(abs(tone));  % Normalize to [-1, 1]
        tone = tone * 0.3;  % Scale to safe playback level (adjust based on your system)

        % Save to WAV file
        filename = fullfile(output_dir, sprintf('tone_%s.wav', tone_names{i}));
        audiowrite(filename, tone, sample_rate);

        fprintf('Generated: %s (%d Hz, 50 ms, 75 dB)\n', filename, freq);
    end

    fprintf('=====================================\n');
    fprintf('All tones generated successfully!\n');
    fprintf('Output directory: %s\n', output_dir);
end
