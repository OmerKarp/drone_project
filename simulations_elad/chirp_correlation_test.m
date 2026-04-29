clear
close all
clc


N = 100000;
fs = 10e7;
ofdm_sym_len = 1024;
SCS = 15e3; %sub carrier spacing
cp_len = 72;
cp_len_ends = 80;
fs_spc_1 = ofdm_sym_len.*SCS;
M = 64;
g = hann(M);
L = 0.5*M;
Ndft = 1024;

t_val = ((0:N-1)./fs ).';
T = (N-1)./fs;

twofourGcarriers = [ .3995e9, .4145e9, .4295e9, .4445e9, .4595e9];
twofourGcarriers = twofourGcarriers-mean(twofourGcarriers);

chirp_bw = SCS.*600;

freq_slope = chirp_bw./(2*T);
chirp = exp(1j*2*pi*t_val.^2.*freq_slope);
freq_offsets = exp(1j*2*pi*t_val*twofourGcarriers);
chirp_freq_offsets = chirp.*freq_offsets;

%chirp = chirp.*exp(1j*2*pi*t_val*(-fs/2));

correlation_peaks = ifft(fft(conj(flip(chirp))) .*fft(chirp_freq_offsets));
figure; plot(max(abs(correlation_peaks)));

[spec,f,t] = spectrogram(chirp,g,L,Ndft,fs,"centered");
figure;
mesh(t,f,abs(spec).^2)
title("spectrogram")
xlabel("t[s]");
ylabel("f[Hz]");
view(2), axis tight