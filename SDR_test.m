
clc
clear
close all
%%
fc = 2.4295e9;
fs = 15.36e6;

%% RX
mcr = 15.36e6; 
decim = round(mcr/fs);    

rx1 = comm.SDRuReceiver('Platform', 'B210', ...
    'SerialNum', '356E654', ...
    'CenterFrequency', fc, ...
    'Gain', 20, ...
    'OutputDataType', 'double', ...
    'MasterClockRate', mcr, ...  
    'SamplesPerFrame',2^10, ...
    'DecimationFactor', decim); 

rx1.ChannelMapping = 2; 
rx1.ReceiveAntennaPort = 'RX2';

% hSpectrum = spectrumAnalyzer(...
%     'Name',             'Passband Spectrum',...
%     'Title',            'Passband Spectrum', ...
%     'Method',           'Welch', ...
%     'SpectrumType',     'Power density', ...
%     'FrequencySpan',    'Full', ...
%     'SampleRate',       fs, ...
%     'FrequencyOffset',  fc, ...
%     'YLimits',          [-140 10], ...
%     'YLabel',           'Magnitude-squared, dB', ...
%     'Position',         figposition([50 30 30 40]));
% LO shift -> moving the high frequency peak at 0 Hz

% Main loop
number_of_frames = 100000;
frames = 0;
ofdm_sym_len = 1024;
SCS = 15e3; %sub carrier spacing
cp_len = 72;
cp_len_extended = 80;
packet_length = 9960;
samp_diff_sym_4 = (2 * (ofdm_sym_len + cp_len) + ofdm_sym_len + cp_len_extended); % sample difference between symbols

buffer = zeros(20*packet_length,1);
saved_last_buffer = zeros(samp_diff_sym_4,1);
papr_thresh = 300;
%samples_test = capture(rx1, 20000*packet_length);

%%
%physical_layer_demod(samples_test,papr_thresh,fs);

%%
while frames < number_of_frames
    %rcv = rx1();
    %step(hSpectrum, rcv);
    %saved_last_buffer = buffer(end-packet_length+1:end);
    buffer = capture(rx1, 20*packet_length);
    %buffer = [saved_last_buffer; buffer];
    raw_hex = physical_layer_demod(buffer,papr_thresh,fs);

    disp(raw_hex);
    frames = frames + 1;
    %pause(0.5);
end

% Release all System objects
% release(sigSrc);
%release(hSpectrum);
release(rx1);



%% 
function raw_hex = physical_layer_demod(signal,papr_thresh,fs)
    % Params
    ofdm_sym_len = 1024;
    SCS = 15e3; %sub carrier spacing
    fs_spc_1 = ofdm_sym_len * SCS;

    % Resample the raw_samples if fs<fs_spc_1
    if fs<fs_spc_1
        signal = resample(signal, fs_spc_1, fs);
        fs = fs_spc_1; %update fs
    end
    % Creating the ZaddOff-Chu sequences
    [zc_symbol_4, zc_symbol_6] = create_zc_symbols();

    % Extracting the CP
    [cp_symbol3,zc_4_peak_samp,is_threshold] = extract_symbol3_cp(zc_symbol_4,papr_thresh, signal);
    if ~is_threshold
        raw_hex = 'No signal Found';
        return;
    end
    %if there are a lot of packets, filter by max value of zc_4 corr peak-
    %use this packet, the code is designed for one packet at a time.
    overhead = 1000; %overhead in samples
    signal = initial_zc_very_coarse_time_sync(signal,zc_4_peak_samp,overhead);

    % Find freq_offset coarse
    freq_offset_coarse = find_freq_offset_coarse(signal, cp_symbol3, fs);

    % Correction with freq_coarse
    signal = fix_freq(signal, freq_offset_coarse, fs);
    
    % Frequency soft fix
    freq_offset_soft = find_freq_soft(signal, zc_symbol_4, zc_symbol_6, fs);
    
    % Correction with freq_soft
    signal = fix_freq(signal, freq_offset_soft, fs);

    % Time AND phase sync
    [samp_offset, phase_offset] = find_samp_and_phase_offset(signal, zc_symbol_4, zc_symbol_6);

    % Sync the signal with phase and time
    packet_length = 9880;
    if(samp_offset<1  || samp_offset + packet_length -1 >length(signal))

    end
    synced_samples = signal(samp_offset : samp_offset + packet_length - 1) .* exp(-1j * phase_offset);

    % Take the synced_samples -> raw_hex for the data_layer (also channel estimation)
    raw_hex = demod_synced_samples(synced_samples, zc_symbol_4);
end

function signal = initial_zc_very_coarse_time_sync(signal,zc_4_peak_samp,overhead)
    ofdm_sym_len = 1024;
    cp_len = 72;
    cp_len_ends = 80;
    packet_length = 9980;
    
    samp_diff_sym_4 = (2* (ofdm_sym_len+cp_len) + ofdm_sym_len+cp_len_ends ); %sample difference between symbols
    first_index = max(zc_4_peak_samp-samp_diff_sym_4-overhead,1);
    last_index = min(zc_4_peak_samp-samp_diff_sym_4+packet_length+overhead,length(signal)-1);
    signal = signal(first_index:last_index);

end
% TESTED
function [samp_offset, phase_offset] = find_samp_and_phase_offset(signal, zc_symbol_4, zc_symbol_6)
    ofdm_sym_len = 1024;
    cp_len = 72;
    cp_len_extended = 80;

    samp_diff_sym_4 = (2 * (ofdm_sym_len + cp_len) + ofdm_sym_len + cp_len_extended); % sample difference between symbols
    signal2corr = [zc_symbol_4 ; zeros(ofdm_sym_len + cp_len, 1) ; zc_symbol_6];
    
    corr_4_6 = do_fast_corr(signal, signal2corr);

    [~,samp_offset_signal2corr] = max(abs(corr_4_6));
    samp_offset = samp_offset_signal2corr - samp_diff_sym_4 - length(signal2corr) + 1;
    phase_offset = angle(corr_4_6(samp_offset_signal2corr));
end

% TESTED
function freq_offset_soft = find_freq_soft(signal, zc_symbol_4, zc_symbol_6, fs)
    ofdm_sym_len = 1024;
    cp_len = 72;

    corr1 = do_fast_corr(signal, zc_symbol_4);
    corr2 = do_fast_corr(signal, zc_symbol_6);
    
    [~,samp1] = max(abs(corr1));
    [~,samp2] = max(abs(corr2));
    phase_diff = angle(corr2(samp2) .* conj(corr1(samp1)));      % we can detect freq changes from -3500 to 3500 hz (1/2*t_diff4_6)
    time_diff_zc_sym = 2 * (ofdm_sym_len + cp_len - 1) ./ fs;    % 2 ofdm_symbols
    freq_offset_soft = phase_diff ./ (2 * pi * time_diff_zc_sym);
end

% TESTED
function correct_freq_signal = fix_freq(raw_samples, freq_offset, fs)
    t_val = ((0:length(raw_samples)-1)./fs ).';
    correct_freq_signal = raw_samples.*exp(-1j * 2 * pi * freq_offset * t_val);
end

% TESTED
function freq_offset_coarse = find_freq_offset_coarse(raw_samples, data_only_cp, fs)
    ofdm_sym_len = 1024;

    corr_with_data_cp = do_fast_corr(raw_samples, data_only_cp);

    sig2corr = [1 ; zeros(1023,1) ; 1];
    corr_2_peaks = do_fast_corr(abs(corr_with_data_cp), sig2corr);

    [~,first_peak_samp] = max(abs(corr_2_peaks));
    first_peak_samp  = first_peak_samp - 1024;
    second_peak_samp = first_peak_samp + ofdm_sym_len;

    freq_offset_coarse = angle(conj(corr_with_data_cp(first_peak_samp))*corr_with_data_cp(second_peak_samp))./(2*pi*(ofdm_sym_len-1)/fs);
end

% TESTED
function [zc_symbol_4, zc_symbol_6] = create_zc_symbols()
    %create raw zc sequences.
    seq1 = zadoffChuSeq(600,601);
    seq2 = zadoffChuSeq(147,601);

    % pad them with zeros at the start and end
    seq1 = [zeros(212,1) ; seq1 ; zeros(1024-813,1)];
    seq2 = [zeros(212,1) ; seq2 ; zeros(1024-813,1)];
    
    % ofdm mod (do ifft)
    seq1 = ifft(ifftshift(seq1));
    seq2 = ifft(ifftshift(seq2));
    
    zc_symbol_4 = [seq1(end-72+1:end) ; seq1];
    zc_symbol_6 = [seq2(end-72+1:end) ; seq2];
end

function [data_only_cp,zc_4_peak_samp,is_threshold] = extract_symbol3_cp(zc_symbol_4,papr_thresh, raw_samples)
    cp_len = 72;
    ofdm_sym_len = 1024;

    corr1 = do_fast_corr(raw_samples, zc_symbol_4);
    [max_corr_value,zc_4_peak_samp] = max(abs(corr1));
    
    first_index = max(zc_4_peak_samp-ofdm_sym_len-cp_len,1);
    last_index = min(zc_4_peak_samp+ofdm_sym_len+cp_len,length(corr1)-1);

    corr1_peak_envirnoment = corr1(first_index:last_index);
    papr = max_corr_value.^2/mean(abs(corr1_peak_envirnoment).^2);
    %disp(papr);
    is_threshold = papr > papr_thresh;
    if ~is_threshold
        data_only_cp=[];
        zc_4_peak_samp = 0;
        return;
    end
    data_cp_samp = zc_4_peak_samp - length(zc_symbol_4) - cp_len + 1;
    data_only_cp = raw_samples(data_cp_samp : data_cp_samp + cp_len - 1);
    
    %normalize samples to regular samples.
    zc_4_peak_samp = zc_4_peak_samp-length(zc_symbol_4)+1;
end

% TESTED
function corr = do_fast_corr(sig1, sig2)
    % make that sig1 is the longer one, else switch
    if length(sig2) > length(sig1)
        temp = sig1;
        sig1 = sig2;
        sig2 = temp;
    end

    % padding with 0's so they match lengths
    sig2 = [conj(flip(sig2)) ; zeros(length(sig1) - length(sig2), 1)];

    ffted_sig1 = fft(sig1);
    ffted_sig2 = fft(sig2);

    corr = ifft(ffted_sig1 .* ffted_sig2);
end

% synced_samples -> raw_hex, using ofdm and qpsk demodulations.
% TESTED
function raw_hex = demod_synced_samples(synced_samples,zc_seq_sym_4)
    demoded_ofdm = ofdm_demod_with_ch_est(synced_samples,zc_seq_sym_4);
    raw_bit_vec = qpskdemod(demoded_ofdm);
    raw_hex = binaryVectorToHex(raw_bit_vec);
end

% TESTED
function demod_channel_est = ofdm_demod_with_ch_est(sig,zc_seq_sym_4)
    cp_len = 72;
    cp_len_extended = 80;
    ofdm_sym_len = 1024;

    %we know symbol 1 is not relevant.
    zc_4_sample = cp_len_extended+ofdm_sym_len+2*(cp_len+ofdm_sym_len)+1;
    zc_4_cand = sig(zc_4_sample:zc_4_sample+cp_len+ofdm_sym_len-1); %zc_4 candidate
    H = fft(zc_4_cand(cp_len+1:end))./fft(zc_seq_sym_4(cp_len+1:end)); %estimate channel
    
    %ofdm symbols 2 - 8
    sig2_8 = sig(ofdm_sym_len+cp_len_extended+1:end-ofdm_sym_len-cp_len_extended); %symbols 2 - 8
    
    %ofdm symbol 9
    sig9 = sig(end-ofdm_sym_len+1:end);                                %symbols 9
    
    demod2_8 = reshape(sig2_8,ofdm_sym_len+cp_len,[]);
    demod2_8 = demod2_8(cp_len+1:end,:);%cp removal
    
    demod2_8 = fftshift(fft(demod2_8)./H,1);
    demod2_8 = demod2_8(:);
    
    demod2_8_no_zc = [demod2_8(1:2*ofdm_sym_len) ; demod2_8(3*ofdm_sym_len+1:4*ofdm_sym_len) ; demod2_8(5*ofdm_sym_len+1:end)];
    demod9 = fftshift(fft(sig9)./H,1);
    
    demod_channel_est = [demod2_8_no_zc ; demod9];
    demod_channel_est = reshape(demod_channel_est,ofdm_sym_len,[]);
    demod_channel_est = [demod_channel_est(213:512,:); demod_channel_est(514:813,:)];
    
    %flatten matrix
    demod_channel_est = demod_channel_est(:);

end

% TESTED
function raw_bit_vec = qpskdemod(sig)
    symbols = zeros(1, length(sig));

%   We demod the data with this demod constalation:
%
%                    |
%         Ⅱ - 01     |   Ⅰ - 00
%                    |
%        ————————————+————————————
%                    |
%         Ⅲ - 11    |   Ⅳ - 10
%                    |

    symbols(real(sig) > 0 & imag(sig) > 0) = 0; % Ⅰ - 00
    symbols(real(sig) < 0 & imag(sig) > 0) = 2; % Ⅱ - 01
    symbols(real(sig) < 0 & imag(sig) < 0) = 3; % Ⅲ - 11 
    symbols(real(sig) > 0 & imag(sig) < 0) = 1; % Ⅳ - 10

    raw_bit_mat_columns = de2bi(symbols, 2, 'left-msb');
    raw_bit_mat_rows = raw_bit_mat_columns.';
    raw_bit_vec = reshape(raw_bit_mat_rows, 1, []);
end