clear
close all
clc

%% open signal
raw_sample_base = 'C:\Users\eladm\Desktop\drone sigint\drone_project\raw_samples-20260423T111818Z-3-001\raw_samples\';
synced_sample_base = 'C:\Users\eladm\Desktop\drone sigint\drone_project\synced_samples\';
signal_path = [raw_sample_base,'raw_samples_12_Feb_2026_13_29_59_531_fs_10MHz.32fc'];
signal_ref_path = [synced_sample_base,'synced_samples_12_Feb_2026_13_29_59_531_fs_15.36MHz.32fc'];

%import raw signal
fid = fopen(signal_path, 'rb'); % 'rb' for read-binary
data = fread(fid, [inf], 'float32');
fclose(fid);
signal = (data(1:2:end) + 1j * data(2:2:end)) ;

%import synced signal
fid = fopen(signal_ref_path, 'rb'); % 'rb' for read-binary
data = fread(fid, [inf], 'float32');
fclose(fid);
signal_ref = (data(1:2:end) + 1j * data(2:2:end)) ;

%% get to know the signal. plots
fs = 10e6;
ofdm_sym_len = 1024;
SCS = 15e3; %sub carrier spacing
cp_len = 72;
cp_len_ends = 80;
fs_spc_1 = ofdm_sym_len.*SCS;
M = 32;
g = hann(M);
L = 0.5*M;
Ndft = 1024;

signal = resample(signal,fs_spc_1,fs);
fs = fs_spc_1;%update fs

%cut the signal to ofdm symbol.
% samp_min = floor(0.059.*fs);
% samp_max = floor(0.061.*fs);
% signal = signal(samp_min:samp_max);

% [spec,f,t] = spectrogram(signal,g,L,Ndft,fs,"centered");
% 
% figure;
% mesh(t,f,abs(spec).^2)
% title("spectrogram")
% xlabel("t[s]");
% ylabel("f[Hz]");
% view(2), axis tight


%% sync signal.

%create zc sequences.
seq1 = zadoffChuSeq(600,601);
seq2 = zadoffChuSeq(147,601);

% seq1(301) = 0;
% seq2(301) = 0;

seq1 = [zeros(212,1);seq1;zeros(1024-813,1)];
seq2 = [zeros(212,1);seq2;zeros(1024-813,1)];

seq1 = (ifft(ifftshift(seq1)));
seq2 = (ifft(ifftshift(seq2)));

seq1 = [seq1(end-72+1:end) ; seq1];
seq2 = [seq2(end-72+1:end) ; seq2];

% fs = 10e6; %smaller because the ends are zeros already - more efficient.
% fs_spc_1 = 1024.*15e3;
% 
% zc_seq_sym_4 = resample(seq1,fs,fs_spc_1);
% zc_seq_sym_6 = resample(seq2,fs,fs_spc_1);

zc_seq_sym_4 = seq1;
zc_seq_sym_6 = seq2;

%%

zc_4_cp = seq1(end-72+1:end);
zc_6_cp = seq2(end-72+1:end);

corr1 = ifft(fft((signal)).*fft([conj(flip((zc_4_cp))); zeros(length(signal)-length(zc_4_cp),1)]));
figure; plot(abs(corr1)); %we get 2 peaks for cp, spaced 1024 samples from each other(T_ofdm_sym);
sig2corr = [1;zeros(1023,1);1];%1024 samples overall
corr11 = ifft(fft(abs(corr1)).*fft([conj(flip((sig2corr))); zeros(length(corr1)-length(sig2corr),1)]));
figure; plot(abs(corr11));
[~,first_peak_samp] = max(abs(corr11));
first_peak_samp = first_peak_samp-1024;
second_peak_samp= first_peak_samp+ofdm_sym_len;
freq_offset_4 = angle(conj(corr1(first_peak_samp))*corr1(second_peak_samp))./(2*pi*ofdm_sym_len/fs);


corr1 = ifft(fft((signal)).*fft([conj(flip((zc_6_cp))); zeros(length(signal)-length(zc_6_cp),1)]));
figure; plot(abs(corr1)); %we get 2 peaks for cp, spaced 1024 samples from each other(T_ofdm_sym);
sig2corr = [1;zeros(1022,1);1];
corr11 = ifft(fft(abs(corr1)).*fft([conj(flip((sig2corr))); zeros(length(corr1)-length(sig2corr),1)]));
figure; plot(abs(corr11));
[~,first_peak_samp] = max(abs(corr11));
first_peak_samp = first_peak_samp-1023;
second_peak_samp= first_peak_samp+ofdm_sym_len;
freq_offset_6 = angle(conj(corr1(first_peak_samp))*corr1(second_peak_samp))./(2*pi*ofdm_sym_len/fs);

%freq_offset = 0.5*(freq_offset_6+freq_offset_4);
freq_offset = freq_offset_4;

t_val = ((0:length(signal)-1)./fs ).';
signal = signal.*exp(-1j*2*pi*freq_offset*t_val);
[freq_offset,phase_offset,phase_ref_samp] = find_freq_phase_offset(signal,zc_seq_sym_4,zc_seq_sym_6,fs);
freq_phase_synced_signal = freq_phase_sync_signal(signal,phase_offset,freq_offset,phase_ref_samp,fs);
synced_signal = samp_sync_signal(freq_phase_synced_signal,1);
synced_signal_time_only = samp_sync_signal(signal,1);

sym_demod = ofdm_demod(synced_signal,cp_len,cp_len_ends,ofdm_sym_len);

%zc_4_cand = signal(phase_ref_samp:phase_ref_samp+ofdm_sym_len+cp_len);


%% ref signal
demod_symbols_ref = ofdm_demod(signal_ref,cp_len,cp_len_ends,ofdm_sym_len);
[freq_offset2,phase_offset,phase_ref_samp] = find_freq_phase_offset(signal_ref,zc_seq_sym_4,zc_seq_sym_6,fs);
freq_phase_synced_signal_ref = freq_phase_sync_signal(signal_ref,phase_offset,freq_offset2,phase_ref_samp,fs);
sym_ref = ofdm_demod(freq_phase_synced_signal_ref,cp_len,cp_len_ends,ofdm_sym_len);


%%
figure;
hold on;
scatterplot(sym_demod);
scatterplot(sym_ref);
hold off;

%% functions
function demod = ofdm_demod(sig,cp_len,cp_len_ends,ofdm_len)
%we know symbol 1 is not relevant.
sig2_8 = sig(ofdm_len+cp_len_ends+1:end-ofdm_len-cp_len_ends); %symbols 2 - 8
sig9 = sig(end-ofdm_len+1:end);                                %symbols 9

demod2_8 = reshape(sig2_8,ofdm_len+cp_len,[]);
demod2_8 = demod2_8(cp_len+1:end,:);%cp removal

demod2_8 = fftshift(fft((demod2_8)),1);
demod2_8 = demod2_8(:);

demod2_8_no_zc = [demod2_8(1:2*ofdm_len) ; demod2_8(3*ofdm_len+1:4*ofdm_len) ; demod2_8(5*ofdm_len+1:end)];
demod9 = fftshift(fft((sig9)));

demod = [demod2_8_no_zc ; demod9];
demod = reshape(demod,ofdm_len,[]);
demod = [demod(213:513,:); demod(515:813,:)];

demod = demod(:);
end

function [freq_offset,phase_offset,phase_ref_samp] = find_freq_phase_offset(signal,zc_seq_sym_4,zc_seq_sym_6,fs)
%zadoff+cyclic_prefix
corr1 = ifft(fft((signal)).*fft([conj(flip((zc_seq_sym_4))); zeros(length(signal)-length(zc_seq_sym_4),1)]));
corr2 = ifft(fft((signal)).*fft([conj(flip((zc_seq_sym_6))); zeros(length(signal)-length(zc_seq_sym_6),1)]));

% figure; plot(abs(corr1));
% figure; plot(abs(corr2));

[~,samp1] = max(abs(corr1));
[~,samp2] = max(abs(corr2));

phase_offset = angle(corr1(samp1));
phase_diff = angle(corr2(samp2).*conj(corr1(samp1))); % we can detect freq changes from -3500 to 3500 hz (1/2*t_diff4_6)

ofdm_sym_len = 1024;
cp_len = 72;
time_diff_zc_sym = 2* (ofdm_sym_len+cp_len)./fs; %2 ofdm_symbols
freq_offset = phase_diff./(2*pi*time_diff_zc_sym);
phase_ref_samp = samp1 - length(zc_seq_sym_4)+1;


end


function samp_offset = find_samp_offset(freq_phase_synced_signal,zc_seq_sym_4,zc_seq_sym_6)
ofdm_sym_len = 1024;
cp_len = 72;
cp_len_ends = 80;
samp_diff_sym_4 = (2* (ofdm_sym_len+cp_len) + ofdm_sym_len+cp_len_ends ); %sample difference between symbols

signal2corr = [zc_seq_sym_4 ;zeros(ofdm_sym_len+cp_len,1); zc_seq_sym_6];

[~,samp_offset_signal2corr] = max(abs(ifft(fft(freq_phase_synced_signal).*fft([conj(flip(signal2corr));zeros(length(freq_phase_synced_signal)-length(signal2corr),1)]))));
samp_offset = samp_offset_signal2corr - samp_diff_sym_4 - length(signal2corr) + 1;

end

function synced_signal = freq_phase_sync_signal(signal,phase_offset,freq_offset,phase_ref_sample,fs)
N = length(signal);
t_val = ((0:N-1) - phase_ref_sample).' ./fs;

fix_sig_freq_fix = exp(1j*(-2*pi*freq_offset.*t_val - phase_offset)); %signal offset
synced_signal = signal.*fix_sig_freq_fix;
end

function synced_signal = samp_sync_signal(freq_phase_synced_signal,samp_offset)
ofdm_sym_len = 1024;
cp_len = 72;
cp_len_ends = 80;
ofdm_packet_len = 7*(ofdm_sym_len+cp_len)+2*(ofdm_sym_len + cp_len_ends);
synced_signal = freq_phase_synced_signal(samp_offset:samp_offset+ofdm_packet_len-1);
end

function zc_seq = create_ofdm_zc(fs,zc_root,zc_len,ofdm_len,SCS,cp_len,usable_carriers)
zc_seq = zadoffChuSeq(zc_root,zc_len);
zc_seq(floor(zc_len/2)+1) = 0;
zc_seq = [zeros(usable_carriers(1)-1,1);zc_seq;zeros(ofdm_len-usable_carriers(2),1)];
zc_seq = (ifft(ifftshift(zc_seq)));
zc_seq = [zc_seq(end-cp_len+1:end) ; zc_seq];

%fs should be smaller because the ends are zeros already - more efficient.
fs_spc_1 = ofdm_len.*SCS;
zc_seq = resample(zc_seq,fs,fs_spc_1);

end

