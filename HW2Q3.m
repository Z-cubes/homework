clear all
close all
clc

% Load data from Excel file
Years = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'UIF, Net Evap Data', 'Range', 'B4:B73');
I1 = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'UIF, Net Evap Data', 'Range', 'C4:C73');

MAWS1 = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'DemandFunctions', 'Range', 'C6:C16');
MV1 = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'DemandFunctions', 'Range', 'D6:D16');
MAWS3 = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'DemandFunctions', 'Range', 'L4:L16');
MV3 = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'DemandFunctions', 'Range', 'M4:M16');

% Convert tables to arrays and convert units for easier calculations
Years = table2array(Years);
I1 = table2array(I1);

MAWS1 = table2array(MAWS1);
MV1 = table2array(MV1);
MV1 = MV1 / 43559.9 / 1e6 * 3.154e7; % milion $ / ft3
MAWS3 = table2array(MAWS3);
MV3 = table2array(MV3);
MV3 = MV3 / 43559.9 / 1e6 * 3.154e7; % milion $ / ft3

% Interpolate demamd data on to a much finer array
MAWS_interp = 0:0.01:1400;
MV1_interp = interp1(MAWS1, MV1, MAWS_interp);
MV3_interp = interp1(MAWS3, MV3, MAWS_interp);

% Initialize variables
Nt = length(Years); % Number of years

for U1original = 0:20:1220
    for U3original = 0:20:1220

        for t = 1:Nt
            U1 = U1original;
            U3 = U3original;
            while true
                % Node 1
                u1(t) = min(U1, I1(t)); % Outflow to city
                Q1(t) = I1(t) - u1(t); % Outflow from node 1

                % Node 2
                R2(t) = 0.45 * u1(t); % Return flow from the city
                Q2(t) = Q1(t) + R2(t); % Outflow from node 2

                % Node 3
                u3(t) = min(U3, Q2(t)); % Outflow to farms
                Q3(t) = Q2(t) - u3(t); % Outflow from node 3

                % Environmental flow deficit test
                if Q3(t) >= 0.5 * I1(t)

                    % Calculate EW1 and EW3 for each year
                    sieve1 = find(MAWS_interp <= u1(t));
                    EW1(t) = trapz(MAWS_interp(sieve1), MV1_interp(sieve1)) - MAWS_interp(sieve1(end)) * MV1_interp(sieve1(end)) - MV1_interp(sieve1(end))* max(U1original - u1(t), 0); % for one year

                    sieve3 = find(MAWS_interp <= u3(t));
                    EW3(t) = trapz(MAWS_interp(sieve3), MV3_interp(sieve3)) - MAWS_interp(sieve3(end)) * MV3_interp(sieve3(end)) - MV3_interp(sieve3(end)) * max(U3original - u3(t), 0); % for one year

                    break; % Advance simulation to the next year
                else
                    % Compute deficit
                    Dfct(t) = 0.5 * I1(t) - Q3(t);
                    % Compute reduction factors
                    F1 = U1 / (U1 + U3);
                    F3 = U3 / (U1 + U3);
                    % Update water share targets
                    U1 = max(U1 - F1 * Dfct(t), 0);
                    U3 = max(U3 - F3 * Dfct(t), 0);       
                end
            end
        end

    end
end


