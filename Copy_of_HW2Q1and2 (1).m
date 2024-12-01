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
MV1 = MV1 / 43559.9 / 1e6 * 3.154e7; % milion $ / ft3 for one year
MAWS3 = table2array(MAWS3);
MV3 = table2array(MV3);
MV3 = MV3 / 43559.9 / 1e6 * 3.154e7; % milion $ / ft3 for one year

% Interpolate demamd data on to a much finer array
MAWS_interp = 0:0.01:1400;
MV1_interp = interp1(MAWS1, MV1, MAWS_interp);
MV3_interp = interp1(MAWS3, MV3, MAWS_interp);

% Initialize variables
Nt = length(Years); % Number of years
U1_range = 0:20:1220; % Range for U1original
U3_range = 0:20:1220; % Range for U3original

% Initialize arrays to store results
Jmean_U1_matrix = zeros(length(U1_range), length(U3_range));
Vmean_U1_matrix = zeros(length(U1_range), length(U3_range));
Jmean_U3_matrix = zeros(length(U1_range), length(U3_range));
Vmean_U3_matrix = zeros(length(U1_range), length(U3_range));

% Simulation loop over ranges of U1original and U3original
for i = 1:length(U1_range)
    for j = 1:length(U3_range)
        U1original = U1_range(i);
        U3original = U3_range(j);
        
        % Reset arrays for each simulation
        u1 = zeros(1, Nt);
        Q1 = zeros(1, Nt);
        R2 = zeros(1, Nt);
        Q2 = zeros(1, Nt);
        u3 = zeros(1, Nt);
        Q3 = zeros(1, Nt);
        EW1 = zeros(1, Nt);
        EW3 = zeros(1, Nt);
        Dfct = zeros(1, Nt);
        
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
                    % Calculate EW1
                    sieve1 = find(MAWS_interp <= u1(t));
                    if ~isempty(sieve1) && length(sieve1) > 1
                        EW1(t) = (trapz(MAWS_interp(sieve1), MV1_interp(sieve1)) - ...
                                MAWS_interp(sieve1(end)) * MV1_interp(sieve1(end)) - ...
                                MV1_interp(sieve1(end)) * max(U1original - u1(t), 0));
                    else
                        EW1(t) = 0; % Default value when no valid sieve is found
                    end

                    % Calculate EW3
                    sieve3 = find(MAWS_interp <= u3(t));
                    if ~isempty(sieve3) && length(sieve3) > 1
                        EW3(t) = (trapz(MAWS_interp(sieve3), MV3_interp(sieve3)) - ...
                                MAWS_interp(sieve3(end)) * MV3_interp(sieve3(end)) - ...
                                MV3_interp(sieve3(end)) * max(U3original - u3(t), 0));
                    else
                        EW3(t) = 0; % Default value when no valid sieve is found
                    end

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
        
        % Compute Jmean and Vmean for U1original
        Jmean_U1 = sum(u1) / Nt;
        Vmean_U1 = sum(EW1) / Nt;

        % Compute Jmean and Vmean for U3original
        Jmean_U3 = sum(u3) / Nt;
        Vmean_U3 = sum(EW3) / Nt;

        % Store results in matrices
        Jmean_U1_matrix(i, j) = Jmean_U1;
        Vmean_U1_matrix(i, j) = Vmean_U1;
        Jmean_U3_matrix(i, j) = Jmean_U3;
        Vmean_U3_matrix(i, j) = Vmean_U3;
    end
end

% Generate a 2D grid for plotting
[U1_grid, U3_grid] = meshgrid(U1_range, U3_range);

% Plot Jmean for U1
figure;
surf(U1_grid, U3_grid, Jmean_U1_matrix');
xlabel('U1_{original}');
ylabel('U3_{original}');
zlabel('J_{mean}(U1)');
title('J_{mean}(U1_{original}, U3_{original})');
colorbar;

% Plot Vmean for U1
figure;
surf(U1_grid, U3_grid, Vmean_U1_matrix');
xlabel('U1_{original}');
ylabel('U3_{original}');
zlabel('V_{mean}(U1)');
title('V_{mean}(U1_{original}, U3_{original})');
colorbar;

% Plot Jmean for U3
figure;
surf(U1_grid, U3_grid, Jmean_U3_matrix');
xlabel('U1_{original}');
ylabel('U3_{original}');
zlabel('J_{mean}(U3)');
title('J_{mean}(U1_{original}, U3_{original})');
colorbar;

% Plot Vmean for U3
figure;
surf(U1_grid, U3_grid, Vmean_U3_matrix');
xlabel('U1_{original}');
ylabel('U3_{original}');
zlabel('V_{mean}(U3)');
title('V_{mean}(U1_{original}, U3_{original})');
colorbar;
