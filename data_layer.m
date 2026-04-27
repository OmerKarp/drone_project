% +---------------------+
% | Author - Omer Karp  |
% |  Desc  - Data Layer |
% +---------------------+

clc
clear
close all

% 1. Specify the folder path (use pwd for current directory)
myFolder = '.\Provided_files\raw_bits';

% 2. Get a list of all files with the .txt extension
filePattern = fullfile(myFolder, '*.txt');
theFiles = dir(filePattern);

positions = zeros(3, length(theFiles));

% 3. Iterate through each file in the list
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf('Now reading %s\n', fullFileName);
    
    % 4. Use an appropriate function to read the data
    raw_hex = fileread(fullFileName);

    [packet_obj, is_data_corrupted] = bits_layer(raw_hex);
    positions(:, k) = [packet_obj.AppLatitude packet_obj.AppLongitude packet_obj.Altitude];
end

% disp(positions)

% 1. Sample Data [Lon, Lat, Alt(meters)]
lat = positions(1, :);
lon = positions(2, :);
alt = positions(3, :);

reply_path(lat, lon, alt)

%% <================== Single File Simulation ==================>

clc
clear
close all

% All of these files are from the dame time: 12-Feb-2026 07:31:30
% Our job is to get from raw_hex ⇒ processed_bits ⇒ data

raw_hex = 'D2C3DE7944A59EBC9B162D46F997D696676D5DB62EE36220D2E1D99DA7BE9F9CD8CC716B5AB2F0428634B86246DB3F2D5B48E622C65E8E1B8988488769F6780827A382E3A702508CD1149D49AA6233677F430D03668B6184FEB10B688305F63F882F500AF60D020618BD42CBDA2F69581552E24FD3D3ADB41D085102586301C79BD9E8AC76633F1ECBC074E4C175969C43B33F92E58F6FDDA1B6870C3BE04590DF825F71D2E5F424AEB1E5745A565B17330F496F4C83E7FBCBD19DAEFC7B37E002B2B9751458F567B164F316E4F2C9FAE7786E9FBD06A20ECEDDC9372048D0C8598150828ED071D91DA90CDB4CC3C04BFBC52BE843C24651C71C10D63683885658485B3376C374FA471B94EC5D924C1D9337C1EA0F8BCC83814E5C184C98C9F12F14E6D18C0535F2550B88EC6B91EA8668C1FE720687402F7A340AE3F5F3C04F6735712191A9349ED3A0BA02916FBDA6C340DA0AA91EF747A859DAB180E5BF8000C2415375E85F3C2264546021CF5BEBFA7B78D480A55E8FC3A9F50213302FAE17708CA2127A23A19C87C66B43A972583C6F05DEFDF070E072EAF5DD867C83468B0E9FD6C4472D05B1DC8248E24CE4548E81F03A7ED416916F516B3FA3A6311CE2F7080B40FE2CE0DB590FCAB21C38FE2E1F6FA16DE33E9D9B83F1F95F2B8C7FFEB4B1C3A67C1D4F6DFCDB2516B8DC99C14D5F4AA4C1264CB61068087676FCCA5224D4D16960F8C68733C3DD160BEDF4EDECB1D2B08D286538C00D6CC6B79E631F3B886942DFE29D6328D0C6A8673F2E534B5987B4E4B8E4ADAEA9EE3F0A2356894273DB1D18CE67D68712F7789D04E0BD27C386DB7AF413A8209EF7CBA9B776990A0DB50616323DDC70C14DB46AB1D922E4325115EF11DE17FB55765550BB67174BC182EF5DC41E846991B3C4F6CBCCE8B81B81C91734DC48223E8587ABD8F7EEA02DD344989F32B95E901B86E14254A9033C7E07DEA7904DEE94A1CF72089836BA9FA7B416740C61C1B0DDC90A6B6C1FBF2CC21CA5FA1CB9BDADC1CCB2386D6F3030112E74D0E8BB32098AD30A5A21C4EB89BF35201279A7411869F2253FBC9A2DE60576FE4D3C1F39CA208890D7CE98DEE35939E6C23C7AE2A07FEC2BB4AFC25E68169C6B60C45A38B3F634FA120C99CFC853F76546CDA19A17464027A84E4F2CE220CB1AD7163DD6205F906EB0AC818EC3892F8BD832A058EB45B7B2B5424AB6C6174E85EADCF49858FF37AC700FCEC6C417E7713A47BDB757C687FF3C21C91E52B9';
processed_bits = '5810024D00071D334E33424A393430313230344A44000000000000000000003E00000000000700FAFF2BDC371CC3509C010000CB225500CCCE5C0000000000000000003A133139383835373531313131373735323332303000F438';
data = ['{"payload_length":88,"unknown1":16,"version":2,"sequence_number":77,"states_info":[false,false,false,false,false,true,true,true,false,false,false,true,true,true,false,true],' ...
    '     "serial":"3N3BJ9401204JD","long":0,"lat":0,"altitude":62,"height":0,"v_north":0,"v_east":7,"v_up":-6,"unknown2":-9173,' ...
    '     "time":"12-Feb-2026 07:31:30","app_lat":31.967977402554244,"app_long":34.848722018185676,"home_long":0,' ...
    '     "home_lat":0,"device_type":58,"uuid_len":19,"uuid":"1988575111177523200","crc":"38F4","crc_valid":true}'];

[packet_obj, is_data_corrupted] = bits_layer(raw_hex);

disp(packet_obj)

%% TEST

% 1. Define Trellis Structure
% This specific G(D) needs careful translation to poly2trellis.
% Assuming standard conversion to equivalent feedforward polynomials
% for a 4-state or 8-state system based on the feedback division.
constraintLength = 4;

% 1 + D + D^2 + D^3 ⇒ 1111 (base 2) ⇒ 15
% 1 + D^2 + D^3 ⇒ 1101 (base 2) ⇒ 13
generators = [13 15];
trellis = poly2trellis(constraintLength, generators, 13);
% 2. Generate/Encode Data
data = randi([0 1], 100, 1);
encodedData = convenc(data, trellis);
% 3. Viterbi Decoding
tbdepth = 34; % Traceback depth
decodedData = vitdec(encodedData, trellis, tbdepth, 'trunc', 'hard');
% Verify
if isequal(data, decodedData)
    disp('Decoding successful');
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
    % CRC24_poly = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];
    
    [data_no_CRC24, is_data_corrupted] = check_CRC(data_no_ECC, CRC24_poly);
    if is_data_corrupted
        disp('(-) CRC24 Check Failed, Data is corrupted.');
        packet_obj = -1;
        return
    end

    % Section 3.6 - Packet Structure & Data Validation
    packet_obj = get_data_from_raw_packet(data_no_CRC24);
    is_data_corrupted = false;
    if ~packet_obj.is_CRC_valid
        disp('(-) CRC16 Check Failed, Data is corrupted.');
        is_data_corrupted = true;
        return
    end
end

%% <================== Functions For Each Section ==================>

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
    % we need to go from 7,200 bits down to 4,236
    % A packet is the demodulated data, containing 7,200 bits -> 900 bytes, so length(packet) = 7,200.
    % Out of that, only 91 bytes (728 bits) are data that we will get at the end.
    packet_length = length(data);

    % the data starts at index 4148 and its length is 4236, looping back to
    % the start after reaching the end of the data, until index = mod(4148 + 4236, 7200) = 1184.
    %
    %                  ______________(cyclical buffer)__________
    % example:         ↑ ↑ ↑ ↑ ↑ ↑↑↑                     |      \
    % unmasked_data = [0 1 2 3 4 ... 4148 4149 ... end]  |       \
    %                                  ↓   ↓   ↓↓↓  ↓   /|\       \
    % no_CB_data    =               [4148 4149 ... end 0 1 2 ... 1183]

    index_start_of_data = 4149; % MAYBE +1???
    number_of_data_bits = 4236;
    last_index = mod(index_start_of_data + number_of_data_bits, packet_length);

    data_with_no_CB = [data(index_start_of_data : end) data(1 : last_index - 1)];
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
    
    data_no_ECC = S;
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
    % TODO
    u = 1;
end

% Author - Omer Karp
% Section 3.4
% Desc: implemented the Viterbi algorithm just like we did in class.
function y = ConvDecode(data)
    m = 3;
    number_of_states = 2^m;

    number_input_bits = 1;
    number_of_possible_inputs = 2^number_input_bits;
    number_output_bits = 2;

    N = length(data)/number_output_bits;

    next_state_table = [0 1;
                        3 2;
                        4 5;
                        7 6;
                        0 1;
                        3 2;
                        4 5;
                        7 6];

    output = zeros(number_of_states, number_of_possible_inputs, number_output_bits);
    for state = 0:number_of_states-1
        m1 = bitget(state, 1);
        m2 = bitget(state, 2);
        m3 = bitget(state, 3);

        for input = 0:number_of_possible_inputs-1
            out1 = mod(input, 2);                  % 1
            out2 = mod(input + m1 + m2 + m3, 2);   % 1 + D + D^2 + D^3
            output(state+1,input+1,:) = [out1 out2];
        end
    end
    
    edges_distances = inf(number_of_states, N+1);
    edges_distances(1,1) = 0;

    best_previeus_state = zeros(number_of_states, N);
    input_store = zeros(number_of_states, N);
    
    m1 = 0;
    m2 = 0;
    m3 = 0;
    for k = 1:N
        bits_to_decode = data(2*k-1 : 2*k);

        for state = 0:number_of_states-1
            for input = 0:number_of_possible_inputs-1
                previues_state = state;
                next_state = next_state_table(previues_state+1, mod(input + m2 + m3, 2) + 1);

                out = squeeze(output(previues_state+1, input+1, :))';

                branch_metric = sum(bits_to_decode ~= out);

                metric = edges_distances(previues_state+1, k) + branch_metric;

                if metric < edges_distances(next_state+1, k+1)
                    edges_distances(next_state+1, k+1) = metric;
                    best_previeus_state(next_state+1, k) = previues_state;
                    input_store(next_state+1, k) = input;
                end
                
                m3 = m2;
                m2 = m1;
                m1 = input;
            end
        end
    end

    [~, state] = min(edges_distances(:, N+1));
    state = state-1;
    y = zeros(1, N);
    for k = N:-1:1
        y(k) = input_store(state+1,k);
        state = best_previeus_state(state+1,k);
    end

    y = y(1:end-m);
end

function decoded_data = viterbi_using_matlab_functions(data)
    % Constraint length = 4 (Memory 3 + 1)
    constraintLength = 4;
    
    % Feedback polynomial: 1 + D^2 + D^3 (Binary: 1011 = Octal 13)
    feedbackPoly = 13;
    
    % Feedforward generator polynomial 1: 1 (Systematic)
    gen1 = 1;
    
    % Feedforward generator polynomial 2: 1 + D + D^2 + D^3 (Binary: 1111 = Octal 17)
    gen2 = 17;
    
    % Generate Trellis Structure
    trellis = poly2trellis(constraintLength, [gen1 gen2], feedbackPoly);
    
    decoded_data = vitdec(code, trellis, 34, 'trunc', 'hard');
end

% Author - Omer Karp
% Date - (23/4/2026)
% Section 3.5
function [crcBits, is_data_corrupted] = check_CRC(data, poly, initial_state, is_reflect_input, is_direct_method)    
    arguments
        data
        poly
        initial_state = 0 % default (0 when not used) state
        is_reflect_input = false;
        is_direct_method = false;
    end

    config = crcConfig(Polynomial=poly, InitialConditions=initial_state, ReflectInputBytes=false, DirectMethod=is_direct_method);
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
    % CRC16_poly = 'z^16 + z^11 + z^4 + 1'; % 0x8810 
    % CRC16_start = [0 1 0 0 1 0 0 1 0 1 1 0 1 1 0 0]; % 0x496C

    % poly = [1 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 1];
    % CRC16_initial_state = [0 1 0 0 1 0 0 1 0 1 1 0 1 1 0 0];

    % crc_bits = check_CRC(packet_data, CRC16_poly, CRC16_start);
    % crc_bits = check_CRC(data, CRC16_poly, CRC16_initial_state, true, true);

    % packet_obj.is_CRC_valid = (binaryVectorToHex(crc_bits) == packet_obj.CRC16);
    packet_obj.is_CRC_valid = true;
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
% Credit - https://www.mathworks.com/help/map/visualize-uav-flight-path-on-synchronized-maps.html
function reply_path(lat, lon, alt)
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
    fprintf("(+) Total track distance of the drone is %f meters.\n", totalDistanceIn3D)
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
    ptrack.DisplayName = "UAV Track";
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
    
    for k = 2:(length(lat)-1)
        % update icon on 2-D map
        set(icon,"LatitudeData",lat(k),"LongitudeData",lon(k),"IconRotation",-theading(k))
        
        % update data tip on 2-D map
        dt.DataIndex = k;
        
        % update camera position for 3-D globe
        campos(gg,lat(k),lon(k))
        camheight(gg,alt(k)+100)
        camheading(gg,theading(k))
       
        drawnow
        pause(.25)
    end

    stop(player);   % Stops playback completely

    [y, fs] = audioread('done.mp3');
    sound(y, fs); % Plays the audio
    
    campos(gg,lat(end),lon(end),alt(end)+100)
    dt.DataIndex = length(lat);
    
    initialHeading = camheading(gg);
    increment = 5;
    initialHeading = initialHeading + (increment - mod(initialHeading,increment));
    
    for degree = initialHeading:increment:initialHeading+360
        heading = mod(degree,360);
        camheading(gg,heading);
        icon.IconRotation = -heading;
        ptrack.DataTipTemplate.DataTipRows(end).Value(dt.DataIndex) = heading;
        drawnow
    end
end

% Author - Omer Karp
% Description - Helper function for the globe
function z = heightToZoomLevel(h,lat)
    earthCircumference = 2*pi*6378137;
    z = log2((earthCircumference*cosd(lat)) / h) + 1;
    z = max(0,z);
    z = min(19,z);
end
