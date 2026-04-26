% Example: Simple CRC Implementation
% msg = [1 1 0 1 0 1 1 0 1 1]; % Message bits
% poly = [1 0 1 1];            % Generator Polynomial (x^3 + x + 1)


data = [1 0 1 0];
poly = [1 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 1];
initial_state = [0 1 0 0 1 0 0 1 0 1 1 0 1 1 0 0];

[~, msgLen] = size(data);
[~, polyLen] = size(poly);

% 1. Append zeros
paddedMsg = [data initial_state];
rem = paddedMsg;

% 2. Division loop
for i = 1:msgLen
    if rem(i) == 1
        % Perform XOR
        rem(i:i+polyLen-1) = mod(rem(i:i+polyLen-1) + poly, 2);
    end
end

% 3. Extract remainder
crc = rem(end-polyLen+2:end);
disp('Generated CRC:');
disp(binaryVectorToHex(crc));
