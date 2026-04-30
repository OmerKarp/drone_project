% ┌—————————————————————————————————————————————┑
% │ Authors - Omer Karp & Ori Harel & Elad Mor  │
% │  Desc   - Full data extraction of DJI drone │
% │  Date   -   The 23 until 30 of April 2026   │
% ┕—————————————————————————————————————————————┙

% ┌———————————————————————————————————————————————————————————————————————————————————————————————┑
% │  The goal of this project is to extract the data from the given drone.                        │
% │  Steps:                                                                                       │
% │                                                                                               │
% │   1) Use the B210 SDR to capture real-time samples                                            │
% │   2) Those samples go through the Physical Layer:                                             │
% │    2.1) Time/Phase/Frequency synchronization                                                  │
% │    2.2) OFDM and QPSK demod                                                                   │
% │  3) The bits (in hex form) go to the Data Layer:                                              │
% │    3.1) They expirence XOR masks, Cyclical buffer manipulations                               │
% │    3.2) Then we apply a Turbo decoder for ECC, using S/P1/P2 with hard viterby                │
% │    3.3) More data and bits manipulations such as interleavers/padding/etc..                   │
% │    3.4) CRC24 check on the bits, packet data extraction and another CRC16 checks              │
% │  4) Then we use the (lat, lon, alt) = (Latitude, Longitude, Altitude) [meters],               │
% │     to define a 3D point on Earth.                                                            │
% │  5) We then use a series of points show a GUI animation in 2D and 3D of                       │
% │     the path the drone took. (We also record this path as a video and sent it via telegram)   │
% │                                                                                               │
% │                      GitHub repo: https://github.com/OmerKarp/drone_project                   │
% ┕———————————————————————————————————————————————————————————————————————————————————————————————┙

%% <================== SDR Layer ==================>

clc
clear
close all

is_recording_path = false;
is_showing_path = false;

% Params
fc = 2.4295e9;
fs = 15.36e6;
packet_length = 9960;
packets_per_buffer = 12;
packet_overlap = 1;
buffer = zeros((packets_per_buffer+packet_overlap)*packet_length,1);
last_buffer = zeros(packet_overlap*packet_length,1);

papr_thresh = 300; % Threshold for detecting a packet, found empirically

% RX
mcr = 15.36e6;
decim = round(mcr/fs);    

rx1 = comm.SDRuReceiver('Platform', 'B210', ...
    'SerialNum', '34D6290', ...
    'CenterFrequency', fc, ...
    'Gain', 20, ...
    'OutputDataType', 'double', ...
    'MasterClockRate', mcr, ...  
    'SamplesPerFrame',2^12, ...
    'DecimationFactor', decim); 

rx1.ChannelMapping = 1;
rx1.ReceiveAntennaPort = 'RX2';

% Main loop
number_of_frames = 10000;
frames = 0;
positions = zeros(3, 100);
position_index = 1;
freq_index = 1;

while frames < number_of_frames
    buffer = [last_buffer ;capture(rx1, packets_per_buffer*packet_length)];
    raw_hex = physical_layer_demod(buffer, papr_thresh, fs);
    last_buffer = buffer(end-packet_overlap*packet_length+1:end);

    if isempty(raw_hex)
        %disp("(-) No signal found")
        [rx1.CenterFrequency, freq_index] = get_new_freq(freq_index); % get the new fc

    elseif position_index < length(positions) % if there is room left
        disp(raw_hex);
        [packet_obj, is_data_corrupted] = bits_layer(raw_hex);

        if is_data_corrupted
            disp("(-) Packet is corrupted!")
            % positions(:, position_index) = positions(:, max(position_index - 1, 1)); % Stay at the last known position
        end

        positions(:, position_index) = [packet_obj.AppLatitude, packet_obj.AppLongitude, packet_obj.Altitude];

        position_index = position_index + 1;
    end

    frames = frames + 1;
end

% Release all System objects
release(rx1);

if is_showing_path
    % 1. Sample Data [Lon, Lat, Alt(meters)]
    lat = positions(1, :);
    lon = positions(2, :);
    alt = positions(3, :);
    
    reply_path(lat, lon, alt, is_recording_path)
end

function [new_freq, new_freq_index]= get_new_freq(freq_index)
    freq_list = [2.3995e9, 2.4145e9, 2.4295e9, 2.4445e9, 2.4595e9];
    new_freq_index = mod(freq_index, 5) + 1;
    new_freq = freq_list(new_freq_index);
end

%% <================== Physical Layer ==================>

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
        raw_hex = [];
        return;
    end
    %if there are a lot of packets, filter by max value of zc_4 corr peak-
    %use this packet, the code is designed for one packet at a time.
    overhead = 1000; %overhead in samples
    signal = initial_zc_very_coarse_time_sync(signal,zc_4_peak_samp,overhead);

    % Find freq_offset coarse
    [freq_offset_coarse,is_boundaries] = find_freq_offset_coarse(signal, cp_symbol3, fs);
    if ~is_boundaries
        raw_hex = [];
        return;
    end
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

    %whole packet is out of boundaries, pass it to the next buffer
    if(samp_offset<1  || samp_offset + packet_length -1 >length(signal))
        raw_hex = [];
        return
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
function [freq_offset_coarse,is_boundaries_valid] = find_freq_offset_coarse(raw_samples, data_only_cp, fs)
    ofdm_sym_len = 1024;
    is_boundaries_valid = 1;
    corr_with_data_cp = do_fast_corr(raw_samples, data_only_cp);

    sig2corr = [1 ; zeros(1023,1) ; 1];
    corr_2_peaks = do_fast_corr(abs(corr_with_data_cp), sig2corr);

    [~,first_peak_samp] = max(abs(corr_2_peaks));
    first_peak_samp  = first_peak_samp - ofdm_sym_len;
    second_peak_samp = first_peak_samp + ofdm_sym_len;
    
    if(first_peak_samp < 1  || second_peak_samp>length(corr_with_data_cp))
        is_boundaries_valid = 0;
        freq_offset_coarse = 0;
        return
    end
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

    % If sample is out of boundries, move to the next buffer
    if data_cp_samp < 1 || data_cp_samp + cp_len - 1 > length(raw_samples)-1
        data_only_cp=[];
        zc_4_peak_samp = 0;
        return;
    end
    data_only_cp = raw_samples(data_cp_samp : data_cp_samp + cp_len - 1);
    
    %normalize samples to regular samples.
    zc_4_peak_samp = zc_4_peak_samp-length(zc_symbol_4)+1;
end

% TESTED
%does corr with fft - last sample should be the peak.
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

%% <================== Bits Layer ==================>

% Author - Omer Karp
% Desc - returns an obj of the packet, and also if the data is corrupt
function [packet_obj, is_data_corrupted] = bits_layer(packet_in_hex)
    % First we transform the hex raw data to bits raw data
    raw_bits = large_hex_data_to_bits(packet_in_hex);

    % Section 3.2 - XOR mask
    XOR_mask = gen_gold_code();
    data_unmasked = xor(raw_bits, XOR_mask.');

    % Section 3.3 - CB (cyclic buffer)
    data_with_no_CB = remove_CB_from_data(data_unmasked);

    % Section 3.4 - Error correction code & Interleaver
    data_no_ECC = remove_error_correction_code_bits(data_with_no_CB);

    % Section 3.5 - CRC
    % We are now left with 1408 bits, the last 24 are CRC, a quick check in
    % Wikipedia showed that this is the CRC-24-Radix-64 version.
    % Source: https://en.wikipedia.org/wiki/Cyclic_redundancy_check#Polynomial_representations
    % The polynomial:
    CRC24_poly = 'z^24 + z^23 + z^18 + z^17 + z^14 + z^11 + z^10 + z^7 + z^6 + z^5 + z^4 + z^3 + z + 1';
    
    [data_no_CRC24, is_data_corrupted] = check_CRC(data_no_ECC, CRC24_poly);
    if is_data_corrupted
        disp('(-) CRC24 Check Failed, Data is corrupted.');
    end

    % Section 3.6 - Packet Structure & Data Validation
    packet_obj = get_data_from_raw_packet(data_no_CRC24);
    if ~packet_obj.is_CRC_valid
        disp('(-) CRC16 Check Failed, Data is corrupted.');
        is_data_corrupted = true;
        return
    end
end

%% <================== Helper Functions - Bits Layer ==================>

% Author - Omer Karp
% Section - Start
function raw_bits = large_hex_data_to_bits(raw_hex)
    % Map each hex digit to its 4 bits binary form
    hex_map = dec2bin(0:15, 4) - '0';
    
    % Convert hex characters to indices (0-15)
    idx = hex2dec(raw_hex(:)) + 1;
    
    % Retrieve bits and reshape into a single row
    raw_bits = reshape(hex_map(idx, :).', 1, []);
end

% Author - Omer Karp
% Section 3.3
function data_with_no_CB = remove_CB_from_data(data)
    % the data starts at index 4148 and ends at index 4149, looping back to
    % the start and ending back at the starting index.

    index_start_of_data = 4149;
    index_last_of_data = 4236;

    data_with_no_CB = [data(index_start_of_data : index_last_of_data) data(1 : index_start_of_data - 1)];
end

% Author - Omer Karp
% Section 3.4.12
function [interleaved_P1, interleaved_P2] = undo_P1_P2_interlace(data)
    interleaved_P1 = data(1:2:end);
    interleaved_P2 = data(2:2:end);
end

% Author - Omer Karp
% Section 3.4.10
function out_data = undo_rectangular_interleaver(in_data)
    out_data_length = 1408;
    dummy_values_length = 28;

    rev_map = [1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31, 2, 18, 10, 26, 6, 22, 14, 30, 4, 20, 12, 28, 8, 24, 16, 32];
    
    data_mat_temp = -1 * ones(45, 32); % intializing an empty matrix
    current_index = 1;
    for column_replacment_index = rev_map
        if column_replacment_index < 29
            data_mat_temp(2:end, column_replacment_index) = in_data(current_index : current_index + 43);
            current_index = current_index + 44;
        else
            data_mat_temp(1:end, column_replacment_index) = in_data(current_index : current_index + 44);
            current_index = current_index + 45;
        end
    end

    reshaped_mat = data_mat_temp.';
    out_data_with_prefix_dummy_value = reshape(reshaped_mat, 1, []);
    out_data = out_data_with_prefix_dummy_value(1 + dummy_values_length : dummy_values_length + out_data_length);
end

% Author - Omer Karp
% Section 3.4
function data_no_ECC = remove_error_correction_code_bits(data)
    % We get "data" as a 4,236 bits array, it is made out of 3 parts, S, P1
    % and P2, interleaved and manipulated inside of it, each part is 4,236 / 3 = 1,412 bits.

    S_interleaved = data(1 : length(data)/3);
    [P1_interleaved, P2_interleaved] = undo_P1_P2_interlace(data(length(data)/3 + 1 : end));

    S = undo_rectangular_interleaver(S_interleaved);
    P1 = undo_rectangular_interleaver(P1_interleaved);
    P2 = undo_rectangular_interleaver(P2_interleaved);
    
    number_of_iterations = 1;
    for iteration_index = 1 : number_of_iterations
        % for each iteration, we interleave S and P1, then use viterby on
        % the product, map it using pi_interleaver and do the same with P2,
        % then do reverse pi and this will be the next S that we can use
        % for the next iteration, while using HARD viterbi, we only do 1
        % itration, when it is a SOFT viterbi, each iteration will improve
        % on the last one, therefor we do ~8 of them.
        
        S_P1_interleaved = do_interleaver(S, P1);
        S_improved = viterbi_using_matlab_functions(S_P1_interleaved);

        S_improved_pi = apply_pi_interleaver(S_improved);
        S_improved_pi_P2_interleaved = do_interleaver(S_improved_pi, P2);
        S_improved_twice_pi = viterbi_using_matlab_functions(S_improved_pi_P2_interleaved);
        S = undo_pi_interleaver(S_improved_twice_pi); % new S for the next iteration
    end
    
    data_no_ECC = S;
end

% Author - Omer Karp
% Section 3.4 (Helper)
function interleaved_vecs = do_interleaver(vec1, vec2)    
    % Stack them vertically: [vec1 ; bvec2] and this gives a 2xN matrix
    % then "linearize" by (:) reads down columns, then across rows
    
    vec1 = vec1(:).';
    vec2 = vec2(:).';
    vecs_mat = [vec1 ; vec2];
    interleaved_vecs = vecs_mat(:).';
end

% Author - Omer Karp
% Section 3.4 (Helper for reversing)
function u_tag = apply_pi_interleaver(u)
    pi_interleaver = mod((43 * (1 : 1408) + 88 * (1 : 1408).^2), 1408);
    u_tag = u(pi_interleaver + 1);
end

% Author - Omer Karp
% Section 3.4
function u = undo_pi_interleaver(u_tag)
    pi_mapping = apply_pi_interleaver(1:1408);
    [~, reverse_mapping] = sort(pi_mapping);
    u = u_tag(reverse_mapping);
end

% Author - Omer Karp
% Section 3.4
function decoded_data = viterbi_using_matlab_functions(data)
    % Constraint length = 4 (Memory 3 + 1)
    constraintLength = 4;

    feedbackPoly = 13;
    gen1 = 13;
    gen2 = 15;

    % Traceback depth (tblen) is typically 5*constraintLength
    tblen = 10 * constraintLength; 
    
    % Generate Trellis Structure
    trellis = poly2trellis(constraintLength, [gen1 gen2], feedbackPoly);
    
    decoded_data = vitdec(data, trellis, tblen, 'trunc', 'hard');
end

% Author - Omer Karp
% Section 3.5
function [crcBits, is_data_corrupted] = check_CRC(data, poly, initial_state, direct_method)    
    arguments
        data
        poly
        initial_state = 0 % default (0 when not used) state
        direct_method logical = false
    end

    config = crcConfig(Polynomial=poly, InitialConditions=initial_state, DirectMethod=direct_method);
    [crcBits, is_data_corrupted] = crcDetect(data.', config);
end

% Author - Omer Karp
% Section 3.6
function packet_obj = get_data_from_raw_packet(data)
    % Only the first 91 bytes have meaningful data.
    % So we cut only the first 91 * 8 = 728 bits.
    packet_data = data(1 : 728);
    packet_obj.RawBits = data(1 : 728);
    
    % In the docs it says "...should be divided by 174533"
    DIVIDING_FACTOR = 174533;

    % We then use the helper function the extract the relevant data to an object.
    packet_obj.PayloadLength = get_value_from_data(packet_data, 1, 0, 'uint');
    packet_obj.Unknown1 = get_value_from_data(packet_data, 1, 1, 'uint');
    packet_obj.Version = get_value_from_data(packet_data, 1, 2, 'uint');
    packet_obj.SequenceNumber = get_value_from_data(packet_data, 2, 3, 'uint');
    packet_obj.StateInfo = get_value_from_data(packet_data, 2, 5, 'info');
    packet_obj.SerialNumber = get_value_from_data(packet_data, 16, 7, 'ASCII', false);
    packet_obj.Longitude = get_value_from_data(packet_data, 4, 23, 'int');
    packet_obj.Latitude = get_value_from_data(packet_data, 4, 27, 'int');
    packet_obj.Altitude = get_value_from_data(packet_data, 2, 31, 'int');
    packet_obj.Height = get_value_from_data(packet_data, 2, 33, 'int');
    packet_obj.VelocityNorth = get_value_from_data(packet_data, 2, 35, 'int');
    packet_obj.VelocityEast = get_value_from_data(packet_data, 2, 37, 'int');
    packet_obj.VelocityUp = get_value_from_data(packet_data, 2, 39, 'int');
    packet_obj.Unknown2 = get_value_from_data(packet_data, 2, 41, 'int');
    packet_obj.Time = get_value_from_data(packet_data, 8, 43, 'time');
    packet_obj.AppLatitude = get_value_from_data(packet_data, 4, 51, 'int') / DIVIDING_FACTOR;
    packet_obj.AppLongitude = get_value_from_data(packet_data, 4, 55, 'int') / DIVIDING_FACTOR;
    packet_obj.HomeLongitude = get_value_from_data(packet_data, 4, 59, 'int');
    packet_obj.HomeLatitude = get_value_from_data(packet_data, 4, 63, 'int');
    packet_obj.DeviceType = get_value_from_data(packet_data, 1, 67, 'int');
    packet_obj.UUIDLength = get_value_from_data(packet_data, 1, 68, 'int');
    packet_obj.UUID = get_value_from_data(packet_data, 20, 69, 'ASCII', false);
    packet_obj.CRC16 = dec2hex(get_value_from_data(packet_data, 2, 89, 'int'));

    % Do the CRC check
    % The polynomial & Initial State are:
    % CRC16_poly = 'z^16 + z^11 + z^4 + 1';
    CRC16_poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];
    CRC16_initial_state = [0 1 0 0 1 0 0 1 0 1 1 0 1 1 0 0];

    bits_mat = reshape(packet_data, 8, []).';
    switched_bits_order_mat = fliplr(bits_mat).';
    bits = reshape(switched_bits_order_mat, 1, []);

    is_direct_method = true;
    [~, is_data_corrupted] = check_CRC(bits, CRC16_poly, CRC16_initial_state, is_direct_method);

    packet_obj.is_CRC_valid = ~is_data_corrupted;
end

% Author - Omer Karp
% Section 3.6 (Helper)
function value = get_value_from_data(packet_data, value_length, offset, value_type, is_little_endian)
    arguments
        packet_data (1,728)
        value_length (1,1)
        offset (1,1)
        value_type string
        is_little_endian logical = true
    end

    byte_length = 8;
    bits = packet_data(1 + offset * byte_length : offset * byte_length + value_length * byte_length);
    
    if is_little_endian
        bits_mat = reshape(bits, 8, []).';
        switched_bits_order_mat = flipud(bits_mat).';
        bits = reshape(switched_bits_order_mat, 1, []);
    end

    if strcmp(value_type, 'ASCII')
        % Convert bits to '0'/'1' characters
        bin_str = char(bits + '0'); 

        % Reshape to 8 bit rows
        bits_mat = reshape(bin_str, 8, []).';

        % Convert each row to 1 ASCII char
        value = char(bin2dec(bits_mat)).';
        
        % We want to remove all the 0x00 bits
        value = deblank(value);

    elseif strcmp(value_type, 'int')
        is_int_signed = true;
        value = bit2int(bits', value_length * byte_length, IsSigned = is_int_signed);

    elseif strcmp(value_type, 'uint')
        is_int_signed = false;
        value = bit2int(bits', value_length * byte_length, IsSigned = is_int_signed);
    
    elseif strcmp(value_type, 'time')
        MILLISECOND_FACTOR = 1000;

        time_decimal = uint64(bi2de(bits, 'left-msb'));
        value = datetime(time_decimal / MILLISECOND_FACTOR, 'ConvertFrom', 'posixtime');
    
    elseif strcmp(value_type, 'info')
        value = logical(bits); % Turning from 1/0 ⇒ true/false
    
    end
end

%% <================== Given Function - Gold Code Generator ==================>

% Author - None (was given)
function seq = gen_gold_code(Nc, L, seed)
    arguments
        Nc (1,1) = 1600
        L (1,1) = 7200
        seed (1,1) = 0x12345678
    end

    x1 = zeros(Nc + L + 31, 1);
    x2 = zeros(Nc + L + 31, 1);
    x2(1:32) = flip(fdec2bin(seed, 32));
    x1(1) = 1;

    for n = 1:(Nc + L)
        x1(n + 31) = xor(x1(n + 3), x1(n));
        x2(n + 31) = xor(xor(x2(n + 3), x2(n + 2)), xor(x2(n+1), x2(n)));
    end

    seq = xor(x1(Nc + 1:Nc+L), x2(Nc + 1:Nc+L));
end

% Author - Omer Karp
% Description - Helper function for gen_gold_code(), returns a list instead
% of a string, example: dec2bin(10) = '1010', fdec2bin(10) = [1 0 1 0]
function bin_result = fdec2bin(dec_munber, min_digits)
    bin_result = dec2bin(dec_munber, min_digits) - '0';
end

% Author - Omer Karp
% Description - Show the path on 2D and 3D, plus reply it as an anomation.
function reply_path(lat, lon, alt, is_record_video)
    mean_Latitude = mean(lat);
    mean_Longitude = mean(lon);
    mean_Altitude = mean(alt);

    figpos = [1000 500 800 400];
    uif = uifigure(Position=figpos, WindowState='maximized');
    ug = uigridlayout(uif,[1,2]);
    p1 = uipanel(ug);
    p2 = uipanel(ug);
    gx = geoaxes(p1,Basemap="satellite"); 
    gg = geoglobe(p2); 
    gx.InnerPosition = gx.OuterPosition;
    gg.Position = [0 0 1 1];
    
    % 2D sync
    heightAboveTerrain = 100;
    gx.MapCenter = [mean_Latitude mean_Longitude];
    zoomLevel = heightToZoomLevel(heightAboveTerrain, mean_Latitude);
    gx.ZoomLevel = zoomLevel;
    
    % 3D sync
    N = egm96geoid(mean_Latitude, mean_Longitude);
    mean_height = mean_Altitude + N;
    ellipsoidalHeight = mean_height + heightAboveTerrain;
    campos(gg, mean_Latitude, mean_Longitude, ellipsoidalHeight)
    drawnow
    
    % Calculate Flight Headings
    wgs84 = wgs84Ellipsoid;
    theading = azimuth(lat(1:end-1), lon(1:end-1), lat(2:end), lon(2:end), wgs84);
    theading = [theading(1); theading(:)];
    
    % Calculate 3-D Distances
    N = egm96geoid(lat, lon);
    h = alt + N;
    
    % Calculate distance offsets
    lat1 = lat(1:end-1);
    lat2 = lat(2:end);
    lon1 = lon(1:end-1);
    lon2 = lon(2:end);
    h1 = h(1:end-1);
    h2 = h(2:end);
    [dx,dy,dz] = ecefOffset(wgs84,lat1,lon1,h1,lat2,lon2,h2);
    
    % Calculate the Euclidean distance between each pair of adjacent points using the hypot function (in meters)
    distanceIncrementIn3D = hypot(hypot(dx, dy), dz);
    
    % Calculate cumulative distance in 3-D and the total distance
    cumulativeDistanceIn3D = cumsum(distanceIncrementIn3D);
    totalDistanceIn3D = sum(distanceIncrementIn3D);
    disp("——————————————————————————————————————————————————————")
    fprintf("(+) Total track distance of the drone is %f meters.\n", totalDistanceIn3D)
    disp("——————————————————————————————————————————————————————")
    tdist = [0 cumulativeDistanceIn3D];
    
    % Ploting Flight Line
    geoplot3(gg, lat, lon, alt,"c", LineWidth=2, HeightReference="geoid")
    ptrack = geoplot(gx, lat, lon, "c", LineWidth=2);
    
    % Center the camera again
    [clat, clon, cheight] = campos(gg);
    gx.MapCenter = [clat clon];
    % gx.ZoomLevel = heightToZoomLevel(cheight,clat);
    drawnow
    
    % Set Initial View
    campos(gg, lat(1), lon(1))
    camheight(gg, alt(1) + 75)
    campitch(gg, -90)
    camheading(gg, theading(3))
    
    % show the start and end locations of the flight track
    hold(gx,"on")
    mstart = geoplot(gx, lat(1), lon(1), "ow", MarkerSize=10, MarkerFaceColor="magenta");
    mend = geoplot(gx, lat(end), lon(end), "ow", MarkerSize=10, MarkerFaceColor="blue");
    icon = geoiconchart(gx, lat(1), lon(1), "drone_image.png", SizeData=30);
    
    icon.DisplayName = "Current Location";
    mstart.DisplayName = "Start Location";
    mend.DisplayName = "End Location";
    ptrack.DisplayName = "Drone Track";
    legend(gx)
    
    % url = "https://basemap.nationalmap.gov/ArcGIS/rest/services/USGSTopo/MapServer/tile/${z}/${y}/${x}/png";
    % addCustomBasemap("usgstopo",url,Attribution="USGS The National Map")
    
    % gx.Basemap = "usgstopo";
    % gx.ZoomLevel = 11;
    
    dt = datatip(ptrack,DataIndex=1,Location="northwest");
    
    ptrack.DataTipTemplate.DataTipRows(1).Format = "%.3f";
    ptrack.DataTipTemplate.DataTipRows(2).Format = "%.3f";
    
    dtrow = dataTipTextRow("Distance",tdist,"%.2f");
    dtrow(end+1) = dataTipTextRow("Altitude",alt,"%.2f");
    dtrow(end+1) = dataTipTextRow("Heading",theading,"%.2f");
    ptrack.DataTipTemplate.DataTipRows(end+1:end+3) = dtrow;
    
    % Reply the fly animation
    pitch = -2.7689;
    campitch(gg,pitch)

    % Start the sound
    [y, fs] = audioread('drone_sound_cut.mp3');
    player = audioplayer(y, fs);
    play(player); % Starts playback

    if is_record_video
        % Create a figure and the writer object
        v = VideoWriter('Drone_path_video.mp4', 'MPEG-4');
        v.FrameRate = 2; % Change FPS here
        open(v);
    end
    
    for k = 2:(length(lat)-1)
        % update icon on 2-D map
        set(icon,"LatitudeData",lat(k),"LongitudeData",lon(k),"IconRotation",-theading(k))
        
        % update data tip on 2-D map
        dt.DataIndex = k;
        
        % update camera position for 3-D globe
        campos(gg,lat(k),lon(k))
        camheight(gg,alt(k)+100)
        camheading(gg,theading(k))
        
        if is_record_video
            frame = getframe(uif);
            writeVideo(v, frame);
        end

        drawnow
        pause(.25)
    end

    stop(player);   % Stops playback completely

    [y, fs] = audioread('done.mp3');
    sound(y, fs); % Plays the audio
    
    campos(gg,lat(end),lon(end),alt(end)+100)
    dt.DataIndex = length(lat);
    
    initialHeading = camheading(gg);
    if is_record_video
        increment = 20;
    else
        increment = 5;
    end
    
    initialHeading = initialHeading + (increment - mod(initialHeading,increment));
    
    for degree = initialHeading:increment:initialHeading+360
        heading = mod(degree,360);
        camheading(gg,heading);
        icon.IconRotation = -heading;
        ptrack.DataTipTemplate.DataTipRows(end).Value(dt.DataIndex) = heading;
    
        if is_record_video
            frame = getframe(uif);
            writeVideo(v, frame);
        end

        drawnow
    end
    
    if is_record_video
        close(v);
        close(uif);
    end
end


% Description - Helper function for the globe
function z = heightToZoomLevel(h,lat)
    earthCircumference = 2*pi*6378137;
    z = log2((earthCircumference*cosd(lat)) / h) + 1;
    z = max(0,z);
    z = min(19,z);
end
