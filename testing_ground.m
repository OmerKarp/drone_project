clear
close all
clc
%%
seq1 = zadoffChuSeq(600,601);
seq1_1 = seq1;
seq2 = zadoffChuSeq(147,601);

seq1(301) = 0;
seq2(301) = 0;

seq1 = [zeros(212,1);seq1;zeros(1024-813,1)];
seq2 = [zeros(212,1);seq2;zeros(1024-813,1)];

seq1 = (ifft(ifftshift(seq1)));
seq2 = (ifft(ifftshift(seq2)));

seq1 = [seq1(end-72+1:end) ; seq1];
seq2 = [seq2(end-72+1:end) ; seq2];

fs = 10e6; %smaller because the ends are zeros already - more efficient.
fs_spc_1 = 1024.*15e3;

seq1 = resample(seq1,fs,fs_spc_1);
seq2 = resample(seq2,fs,fs_spc_1);

M = 64;
g = hann(M,"periodic");
L = 0.5*M;
Ndft = 1024;

spectrogram(seq1,g,L,Ndft,fs,"centered","yaxis");

