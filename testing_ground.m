clear
close all
clc

%% viterbi


constraintLength = 4;
data = [0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];

% 1 + D + D^2 + D^3 ⇒ 1111 (base 2) ⇒ 15
% 1 + D^2 + D^3 ⇒ 1101 (base 2) ⇒ 13
generators = [13 15];
trellis = poly2trellis(constraintLength, generators, 13);

% 2. Generate/Encode Data
encodedData = convenc(data, trellis);
% 3. Viterbi Decoding
tbdepth = 34; % Traceback depth
decodedData = vitdec(encodedData, trellis, tbdepth, 'trunc', 'hard');
% Verify
if isequal(data, decodedData)
    disp('Decoding successful');
end

