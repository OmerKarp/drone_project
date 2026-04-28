% +---------------------+
% | Author - Omer Karp  |
% |  Desc  - SDR Test   |
% +---------------------+

clc
clear
close all

fc = 2.4145e9;
fs = 1e6;

connectedRadios = findsdru();
disp(connectedRadios);

%% Rx

% Create SDRu Receiver System object
rx = comm.SDRuReceiver(...
    'Platform', 'B210', ...
    'SerialNum', '34D62B0', ...
    'CenterFrequency', fc, ...
    'Gain', 40, ...
    'MasterClockRate', 56e6, ...
    'DecimationFactor', 256);

% Setup signal logging
rxLog = dsp.SignalSink;

% Capture loop
for counter = 1:20
    data = rx(); % Receive data
    rxLog(data); % Log data
end

% Release hardware
release(rx);
% release(rxLog)

%% Received data is in rxLog.Buffer

% --- USRP Spectrum Analyzer Setup ---
% 1. Create the receiver object
rx = comm.SDRuReceiver('Platform', 'B210', ... % Change to your platform (e.g., N210, X310)
    'SerialNum', '34D62B0', ...               % Change to your radio serial number
    'CenterFrequency', 2483.5e6, ...             % Example: 2.4 GHz
    'Gain', 30, ...                           % Set gain
    'SampleRate', 1e6);                       % 1 MHz bandwidth

% 2. Initialize the Spectrum Analyzer
spectrumScope = spectrumAnalyzer('SampleRate', rx.SampleRate, ...
    'ViewType', 'spectrum', ...
    'Title', 'USRP Spectrum Analyzer');

% 3. Main Loop: Capture and Visualize
numFrames = 1000;
for i = 1:numFrames
    % Receive data
    data = rx();
    
    % Send data to Spectrum Analyzer
    spectrumScope(data);
end

% 4. Cleanup
release(rx);
release(spectrumScope);

%% das
rx = comm.SDRuReceiver('Platform', 'B210', ... % Change to your platform (e.g., N210, X310)
    'SerialNum', '34D62B0', ...               % Change to your radio serial number
    'CenterFrequency', fc, ...             % Example: 2.4 GHz
    'Gain', 30, ...                           % Set gain
    'SampleRate', fs);                       % 1 MHz bandwidth

hSpectrum = spectrumAnalyzer(...
    'Name',             'Passband Spectrum',...
    'Title',            'Passband Spectrum', ...
    'Method',           'Welch', ...
    'SpectrumType',     'Power density', ...
    'FrequencySpan',    'Full', ...
    'SampleRate',       fc, ...
    'FrequencyOffset',  fc, ...
    'YLimits',          [-120 10], ...
    'YLabel',           'Magnitude-squared, dB', ...
    'Position',         figposition([50 30 30 40]));

% Initialize radio time
radioTime = 0;

% Main loop
number_of_frames = 10000;
frames = 0;
while frames < number_of_frames
    rcv = rx();
  
    rcv = double(rcv) - mean(abs(double(rcv)));  % Remove DC component.
    step(hSpectrum, rcv);

    frames = frames + 1;
end

% Release all System objects
% release(sigSrc);
release(hSpectrum);
release(rx);

%% F

function set_new_fc()
    posibble_freqs = [5.7965e9, 5.7765e9, 5.7565e9, 2.3995e9, 2.4145e9, 2.4295e9, 2.4445e9, 2.4595e9];

end