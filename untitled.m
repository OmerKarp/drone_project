%% open signal
raw_smaple_base = 'C:\Users\eladm\Desktop\drone sigint\drone_project\raw_samples-20260423T111818Z-3-001\raw_samples\';
synced_sample_base = 'C:\Users\eladm\Desktop\drone sigint\drone_project\synced_samples\';
signal_path = [synced_sample_base,'synced_samples_12_Feb_2026_09_30_39_442_fs_15.36MHz.32fc'];


% 1. Open the file
fid = fopen(signal_path, 'rb'); % 'rb' for read-binary

% 2. Read the data

data = fread(fid, [inf], 'float32');

% 3. Close the file
fclose(fid);

% 4. Convert to complex MATLAB data
% sig = complex(data(1,:), data(2,:)).';
sig = data(1:2:end) + 1j * data(2:2:end);

ofdm_len = 1024;
SCS = 15e3; %sub carrier spacing
cp_len = 72;
cp_len_ends = 80;
fs_spc_1 = ofdm_len.*SCS;
Ndft = 1024;




fs = fs_spc_1;

sig2_8 = sig(ofdm_len+cp_len_ends+1:end-ofdm_len-cp_len_ends); %symbols 2 - 8
sig9 = sig(end-ofdm_len+1:end);                                  %symbols 9

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
demod2 = ofdm_demod(sig,cp_len,cp_len_ends,ofdm_len);
% demod_builtin = ofdmdemod(sig2_8,1024,72);
% demod_builtin = demod_builtin(:);



function demod = ofdm_demod(sig,cp_len,cp_len_ends,ofdm_len)
%we know symbol 1 is not relevant.
sig2_8 = sig(ofdm_len+cp_len_ends+1:end-ofdm_len-cp_len_ends); %symbols 2 - 8
sig9 = sig(end-ofdm_len+1:end);                                  %symbols 9

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