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
MV1 = MV1 / 43559.9 / 1e6; % milion $ / ft3
MAWS3 = table2array(MAWS3);
MV3 = table2array(MV3);
MV3 = MV3 / 43559.9 / 1e6; % milion $ / ft3

% Interpolate demamd data on to a much finer array
MAWS_interp = 0:0.01:1400;
MV1_interp = interp1(MAWS1, MV1, MAWS_interp);
MV3_interp = interp1(MAWS3, MV3, MAWS_interp);
figure(1)
plot(MAWS1, MV1)
figure(2)
plot(MAWS_interp, MV1_interp)


% Initialize variables
N = length(Years); % Number of years
U1original = 800; % Initial target water share for the city (adtust as needed)
U3original = 1100; % Initial target water share for the farms (adtust as needed)

% Simulation loop arbitrary watershares
for t = 1:N
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
            
            % look at interp data up to specified point 
            sieve1 = find(MAWS_interp <= u1(t));
            EW1(t) = trapz(MAWS_interp(sieve1), MV1_interp(sieve1));

            sieve3 = find(MAWS_interp <= u3(t));
            EW3(t) = trapz(MAWS_interp(sieve3), MV3_interp(sieve3));

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

% Calculate the mean values
mean_u1 = mean(u1)
mean_u3 = mean(u3)


% Plot Simulated flow for UIF, u1, u3, and Q3
figure(3) 
hold on 
plot(Years, I1, 'c')
plot(Years, u1, 'r')
plot(Years, u3, 'b')
plot(Years, Q3, 'g')

legend('UIF', 'u1(city)', 'u3(farm)', 'Q3')
ylabel('Flows, Water shares (cfs)')
xlabel('Years')
hold off

% Rank Flows, water shares
rank_I1 = sort(I1, 'descend');
rank_u1 = sort(u1, 'descend');
rank_u3 = sort(u3, 'descend');
rank_Q3 = sort(Q3, 'descend');

% Plot exceedance frequency for Simulated flow for UIF, u1, u3, and Q3
prob_exceedance = (1:N) / (N+1);
figure(4)
hold on 
plot(prob_exceedance(1,1:70), rank_I1, 'c')
plot(prob_exceedance(1,1:70), rank_u1, 'r')
plot(prob_exceedance(1,1:70), rank_u3, 'b')
plot(prob_exceedance(1,1:70), rank_Q3, 'g')
hold off

% Plot EW1 and EW3
figure(5)
hold on 
plot(Years, EW1)
plot(Years, EW3)
hold off

% Create Results table 
% ResultsTable = table(Years, I1, u1, Q1, R2, Q2, u3, Q3, ...
%     'VariableNames', {'TimeIndex', 'Inflow_I1', 'OutflowToCity_u1', ...
%     'OutflowFromNode1_Q1', 'ReturnFlow_R2', 'OutflowFromNode2_Q2', ...
%     'OutflowToFarms_u3', 'OutflowFromNode3_Q3'});
% 
% disp(ResultsTable);

% code Grave yard 

% % Use water supply to solve consumer surplus 
% for t = 1:N
%     U1 = U1original;
%     U3 = U3original;
%     for i = 1:length(MAWS_city)
%         if u1(t) >= MAWS_city(i)
%             Ci1 = i; 
%             Ci2 = i + 1;
%             if Ci2 > length(MAWS_farm)
%                 Ci1 = Ci1 - 1;
%                 Ci2 = Ci2 - 1;
%             end 
%         end 
%     end
%     % find where u1 is on the City Demand Function 
%     slope1 = (MV_city(Ci2) - MV_city(Ci1) / MAWS_city(Ci2) - MAWS_city(Ci1));
%     b1 = MV_city(Ci1) - (slope1 * MAWS_city(Ci1));
%     Marg_Val1(t) = u1(t) * slope1 + b1;
%     Con_Surplus1(t) = trapz(MAWS_city(1:Ci1), MV_city(1:Ci1)) - (MAWS_city(Ci1) * MV_city(Ci1));
%      if u1(t) >= U1
%         EW1(t) =  Con_Surplus1(t);
%     else
%         EW1(t) = Con_Surplus1(t) - (U1 - u1(t)) * Marg_Val1(t);
%     end
% 
% 
%     for i = 1:length(MAWS_farm)
%         if u3(t) >= MAWS_farm(i)
%             Fi1 = i; 
%             Fi2 = i + 1;
%             if Fi2 > length(MAWS_farm)
%                 Fi1 = Fi1 - 1;
%                 Fi2 = Fi2 - 1;
%             end 
%         end 
%     end
%     % find where u3 is on the farm Demand Function 
%     slope3 = (MV_farm(Fi2) - MV_farm(Fi1) / MAWS_farm(Fi2) - MAWS_farm(Fi1));
%     b3 = MV_farm(Fi1) - (slope3 * MAWS_farm(Fi1));
%     Marg_Val3(t) = u3(t) * slope3 + b3;
%     Con_Surplus3(t) = trapz(MAWS_farm(1:Fi1), MV_farm(1:Fi1)) - (MAWS_farm(Ci1) * MV_farm(Ci1));
%     if u3(t) >= U1
%         EW3(t) =  Con_Surplus3(t);
%     else 
%         EW3(t) = Con_Surplus3(t) - (U1 - u3(t)) * Marg_Val3(t);
%     end 
% end
% 
% % Calculate Market Price 
% EW1 = EW1 * 365.25 * 24 * 60 * 60 / (10^6) / (43559.9);
% EW3 = EW3 * 365.25 * 24 * 60 * 60 / (10^6) / (43559.9);
% 
% % Plot Economic welfere for water shares 
% figure(3)
% hold on 
% plot(Years, EW1, 'r')
% plot(Years, EW3, 'b')
% hold off