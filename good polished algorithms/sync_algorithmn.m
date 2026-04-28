clear
close all
clc

%% open signal
raw_sample_base = 'C:\Users\eladm\Desktop\drone sigint\drone_project\raw_samples-20260423T111818Z-3-001\raw_samples\';
synced_sample_base = 'C:\Users\eladm\Desktop\drone sigint\drone_project\synced_samples\';
signal_path = [raw_sample_base,'raw_samples_12_Feb_2026_13_34_20_590_fs_10MHz.32fc'];
signal_ref_path = [synced_sample_base,'synced_samples_12_Feb_2026_13_34_20_590_fs_15.36MHz.32fc'];

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
signal = signal.*exp(1j*2*pi*0*((0:length(signal)-1)./fs).');

%create zc sequences.
seq1 = zadoffChuSeq(600,601);
seq2 = zadoffChuSeq(147,601);

% seq1(301) = 0;
% seq2(301) = 0;

seq1 = [zeros(212,1);seq1;zeros(ofdm_sym_len-813,1)];
seq2 = [zeros(212,1);seq2;zeros(ofdm_sym_len-813,1)];

seq1 = (ifft(ifftshift(seq1)));
seq2 = (ifft(ifftshift(seq2)));

seq1 = [seq1(end-cp_len+1:end) ; seq1];
seq2 = [seq2(end-cp_len+1:end) ; seq2];
zc_seq_sym_4 = seq1;
zc_seq_sym_6 = seq2;

zc_4_cp = seq1(end-cp_len+1:end);
zc_6_cp = seq2(end-cp_len+1:end);

corr1 = ifft(fft((signal)).*fft([conj(flip((seq1))); zeros(length(signal)-length(seq1),1)]));

figure; plot(abs(corr1)); 
[~,zc_4_peak_samp] = max(abs(corr1));
data_cp_samp = zc_4_peak_samp-length(seq1)-cp_len+2;
data_cp = signal(data_cp_samp:data_cp_samp+cp_len-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

corr1 = ifft(fft((signal)).*fft([conj(flip((data_cp))); zeros(length(signal)-length(data_cp),1)]));
figure; plot(abs(corr1)); %we get 2 peaks for cp, spaced 1024 samples from each other(T_ofdm_sym);
sig2corr = [1;zeros(1023,1);1];
corr11 = ifft(fft(abs(corr1)).*fft([conj(flip((sig2corr))); zeros(length(corr1)-length(sig2corr),1)]));
figure; plot(abs(corr11));
[~,first_peak_samp] = max(abs(corr11));
first_peak_samp = first_peak_samp-1024;
second_peak_samp= first_peak_samp+ofdm_sym_len;
freq_offset_corase = angle(conj(corr1(first_peak_samp))*corr1(second_peak_samp))./(2*pi*(ofdm_sym_len-1)/fs);

t_val = ((0:length(signal)-1)./fs ).';
signal = signal.*exp(-1j*2*pi*freq_offset_corase*t_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corr1 = ifft(fft((signal)).*fft([conj(flip((zc_seq_sym_4))); zeros(length(signal)-length(zc_seq_sym_4),1)]));
corr2 = ifft(fft((signal)).*fft([conj(flip((zc_seq_sym_6))); zeros(length(signal)-length(zc_seq_sym_6),1)]));

% figure; plot(abs(corr1));
% figure; plot(abs(corr2));

[~,samp1] = max(abs(corr1));
[~,samp2] = max(abs(corr2));
phase_diff = angle(corr2(samp2).*conj(corr1(samp1))); % we can detect freq changes from -3500 to 3500 hz (1/2*t_diff4_6)
time_diff_zc_sym = 2* (ofdm_sym_len+cp_len-1)./fs; %2 ofdm_symbols
freq_offset_soft = phase_diff./(2*pi*time_diff_zc_sym)-25;
t_val = ((0:length(signal)-1)./fs ).';
signal = signal.*exp(-1j*2*pi*freq_offset_soft*t_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


samp_diff_sym_4 = (2* (ofdm_sym_len+cp_len) + ofdm_sym_len+cp_len_ends ); %sample difference between symbols
signal2corr = [zc_seq_sym_4 ;zeros(ofdm_sym_len+cp_len,1); zc_seq_sym_6];

corr_4_6 = ifft(fft(signal).*fft([conj(flip(signal2corr));zeros(length(signal)-length(signal2corr),1)]));
[~,samp_offset_signal2corr] = max(abs(corr_4_6));
figure; plot(abs(corr_4_6));
title("corr_4_6")
samp_offset = samp_offset_signal2corr - samp_diff_sym_4 - length(signal2corr) + 1;
phase_offset_exp = corr_4_6(samp_offset_signal2corr)./abs(corr_4_6(samp_offset_signal2corr));



%2D frequency-time search would also work here with signal2corr because
%frequency shifts affect different zc sequences(chirps) differently!!

synced_sig_coarse = signal(samp_offset:samp_offset+9880-1)*conj(phase_offset_exp);
[spec,f,t] = spectrogram(synced_sig_coarse,g,L,Ndft,fs,"centered");
figure;
mesh(t,f,abs(spec).^2)
title("spectrogram")
xlabel("t[s]");
ylabel("f[Hz]");
view(2), axis tight

[spec,f,t] = spectrogram(signal_ref,g,L,Ndft,fs,"centered");
figure;
mesh(t,f,abs(spec).^2)
title("spectrogram")
xlabel("t[s]");
ylabel("f[Hz]");
view(2), axis tight

figure;
hold on;
plot(abs(synced_sig_coarse));
plot(abs(signal_ref));
hold off;

figure;
hold on;
plot(angle(synced_sig_coarse));
plot(angle(signal_ref));
hold off;

signal_ref = signal_ref.*exp(1j*2*pi*10*((0:length(signal_ref)-1)./fs).');
sym_ref = ofdm_demod(signal_ref,cp_len,cp_len_ends,ofdm_sym_len);
sym_demod = ofdm_demod(synced_sig_coarse,cp_len,cp_len_ends,ofdm_sym_len);

figure;
hold on;
scatterplot(sym_demod);
scatterplot(sym_ref);
hold off;

function demod = ofdm_demod(sig,cp_len,cp_len_ends,ofdm_sym_len)
%we know symbol 1 is not relevant.
sig2_8 = sig(ofdm_sym_len+cp_len_ends+1:end-ofdm_sym_len-cp_len_ends); %symbols 2 - 8
sig9 = sig(end-ofdm_sym_len+1:end);                                %symbols 9

demod2_8 = reshape(sig2_8,ofdm_sym_len+cp_len,[]);
demod2_8 = demod2_8(cp_len+1:end,:);%cp removal

demod2_8 = fftshift(fft((demod2_8)),1);
demod2_8 = demod2_8(:);

demod2_8_no_zc = [demod2_8(1:2*ofdm_sym_len) ; demod2_8(3*ofdm_sym_len+1:4*ofdm_sym_len) ; demod2_8(5*ofdm_sym_len+1:end)];
demod9 = fftshift(fft((sig9)));

demod = [demod2_8_no_zc ; demod9];
demod = reshape(demod,ofdm_sym_len,[]);
demod = [demod(213:513,:); demod(515:813,:)];

demod = demod(:);
end

