clc
clear
close all

known_raw_hex = 'D2C3DE7944A59EBC9B162D46F997D696676D5DB62EE36220D2E1D99DA7BE9F9CD8CC716B5AB2F0428634B86246DB3F2D5B48E622C65E8E1B8988488769F6780827A382E3A702508CD1149D49AA6233677F430D03668B6184FEB10B688305F63F882F500AF60D020618BD42CBDA2F69581552E24FD3D3ADB41D085102586301C79BD9E8AC76633F1ECBC074E4C175969C43B33F92E58F6FDDA1B6870C3BE04590DF825F71D2E5F424AEB1E5745A565B17330F496F4C83E7FBCBD19DAEFC7B37E002B2B9751458F567B164F316E4F2C9FAE7786E9FBD06A20ECEDDC9372048D0C8598150828ED071D91DA90CDB4CC3C04BFBC52BE843C24651C71C10D63683885658485B3376C374FA471B94EC5D924C1D9337C1EA0F8BCC83814E5C184C98C9F12F14E6D18C0535F2550B88EC6B91EA8668C1FE720687402F7A340AE3F5F3C04F6735712191A9349ED3A0BA02916FBDA6C340DA0AA91EF747A859DAB180E5BF8000C2415375E85F3C2264546021CF5BEBFA7B78D480A55E8FC3A9F50213302FAE17708CA2127A23A19C87C66B43A972583C6F05DEFDF070E072EAF5DD867C83468B0E9FD6C4472D05B1DC8248E24CE4548E81F03A7ED416916F516B3FA3A6311CE2F7080B40FE2CE0DB590FCAB21C38FE2E1F6FA16DE33E9D9B83F1F95F2B8C7FFEB4B1C3A67C1D4F6DFCDB2516B8DC99C14D5F4AA4C1264CB61068087676FCCA5224D4D16960F8C68733C3DD160BEDF4EDECB1D2B08D286538C00D6CC6B79E631F3B886942DFE29D6328D0C6A8673F2E534B5987B4E4B8E4ADAEA9EE3F0A2356894273DB1D18CE67D68712F7789D04E0BD27C386DB7AF413A8209EF7CBA9B776990A0DB50616323DDC70C14DB46AB1D922E4325115EF11DE17FB55765550BB67174BC182EF5DC41E846991B3C4F6CBCCE8B81B81C91734DC48223E8587ABD8F7EEA02DD344989F32B95E901B86E14254A9033C7E07DEA7904DEE94A1CF72089836BA9FA7B416740C61C1B0DDC90A6B6C1FBF2CC21CA5FA1CB9BDADC1CCB2386D6F3030112E74D0E8BB32098AD30A5A21C4EB89BF35201279A7411869F2253FBC9A2DE60576FE4D3C1F39CA208890D7CE98DEE35939E6C23C7AE2A07FEC2BB4AFC25E68169C6B60C45A38B3F634FA120C99CFC853F76546CDA19A17464027A84E4F2CE220CB1AD7163DD6205F906EB0AC818EC3892F8BD832A058EB45B7B2B5424AB6C6174E85EADCF49858FF37AC700FCEC6C417E7713A47BDB757C687FF3C21C91E52B9';

%% open signal
raw_sample_base = '.\Provided_files\raw_samples\';
synced_sample_base = '.\Provided_files\synced_samples\';

signal_path = [raw_sample_base,'raw_samples_12_Feb_2026_09_30_39_442_fs_10MHz.32fc'];
signal_sync_path = [synced_sample_base,'synced_samples_12_Feb_2026_09_30_39_442_fs_15.36MHz.32fc'];

%import raw signal
fid = fopen(signal_path, 'rb'); % 'rb' for read-binary
data = fread(fid, [inf], 'float32');
fclose(fid);
signal = (data(1:2:end) + 1j * data(2:2:end)) ;

%import synced signal
fid = fopen(signal_sync_path, 'rb'); % 'rb' for read-binary
data = fread(fid, [inf], 'float32');
fclose(fid);
known_synced_signal = (data(1:2:end) + 1j * data(2:2:end)) ;

%% <================== Physical Layer ==================>

raw_hex = physical_layer_demod(signal);

% disp(1800 - nnz(raw_hex == known_raw_hex))

% a = demod_synced_samples(known_synced_signal)
b = hexToBinaryVector(raw_hex);
known_bits = hexToBinaryVector(known_raw_hex);

find(b ~= known_bits)

%% <================== Physical Layer ==================>

% Handle the raw samples and turn to synced samples
function raw_hex = physical_layer_demod(raw_samples)
    % Params
    fs = 10e6;
    ofdm_sym_len = 1024;
    SCS = 15e3; %sub carrier spacing
    fs_spc_1 = ofdm_sym_len .* SCS;

    % Resample the raw_samples
    raw_samples = resample(raw_samples, fs_spc_1, fs);
    fs = fs_spc_1; %update fs
    
    % Creating the ZaddOff-Chu sequences
    [zc_symbol_4, zc_symbol_6] = create_zc_symbols();

    % Extracting the CP
    data_only_cp = extract_data_cp(zc_symbol_4, raw_samples);

    % Find freq_offset coarse
    freq_offset_coarse = find_freq_offset_coarse(raw_samples, data_only_cp, fs);

    % Correction with freq_coarse
    signal = fix_freq(raw_samples, freq_offset_coarse, fs);
    
    % Frequency soft fix
    freq_offset_soft = find_freq_soft(signal, zc_symbol_4, zc_symbol_6, fs);
    
    % Correction with freq_soft
    signal = fix_freq(signal, freq_offset_soft, fs);

    % Time AND phase sync
    [samp_offset, phase_offset] = find_samp_and_phase_offset(signal, zc_symbol_4, zc_symbol_6);

    % Sync the signal with phase and time
    packet_length = 9880;
    synced_samples = signal(samp_offset : samp_offset + packet_length - 1) .* exp(-1j * phase_offset);

    % Channel estimation
    % synced_samples_with_channel_estimation = fix_channel_estimation(synced_samples)

    % Take the synced_samples -> raw_hex for the data_layer
    raw_hex = demod_synced_samples(synced_samples,zc_symbol_4);
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

% TESTED
function data_only_cp = extract_data_cp(zc_symbol_4, raw_samples)
    cp_len = 72;

    corr1 = do_fast_corr(raw_samples, zc_symbol_4);

    [~,zc_4_peak_samp] = max(abs(corr1));
    data_cp_samp = zc_4_peak_samp - length(zc_symbol_4) - cp_len + 1;
    data_only_cp = raw_samples(data_cp_samp : data_cp_samp + cp_len - 1);
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
    
    sig2_8 = sig(ofdm_sym_len+cp_len_extended+1:end-ofdm_sym_len-cp_len_extended); %symbols 2 - 8
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

%% <================== Helper Functions - Plotting ==================>

function plot_IQ(z, show_radius, plot_title)
    arguments
        z
        show_radius (1,1) logical = false   % optinal param to see the radiuses
        plot_title string = ""              % optinal suffix name to the graph
    end

    % getting the radiuses
    radiuses = abs(z);
    radiuses = round(radiuses, 2);
    radiuses = unique(radiuses);

    figure;
    hold on;
    xline(0, 'k-');
    yline(0, 'k-');

    plot(real(z), imag(z), 'o', 'LineWidth', 10)
    if show_radius
        theta = linspace(0, 2*pi, 100);

        for raduis = radiuses
            z = raduis * exp(1j * theta);
            plot(real(z), imag(z), 'Color', 'r', 'LineStyle', '--');
        end
    end
    
    title('IQ Graph ' + plot_title);
    xlabel('I (real)');
    ylabel('Q (Imaginary)');
    grid on;
end

function plot_two_sided_fft(y, fs, plot_title)
    arguments
        y 
        fs 
        plot_title string = "" % optinal suffix for the graph
    end

    N = length(y);
    ffted_y = fft(y);
    P2 = abs(ffted_y/N); % 2 sided ffted data
    P2_shifted = fftshift(P2);
    freq_axis = fs * (-N/2 : N/2 - 1) / N;

    figure;
    plot(freq_axis, P2_shifted);
    title("Two Sided FFT " + plot_title);
    xlabel("Frequency (Hz)");
    ylabel("Amplitude");
    grid on;
end