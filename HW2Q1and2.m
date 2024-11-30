clear all
close all
clc

% Load data from Excel file
Years = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'UIF, Net Evap Data', 'Range', 'B4:B73');
I1 = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'UIF, Net Evap Data', 'Range', 'C4:C73');

% Convert tables to arrays for easier calculations
Years = table2array(Years);
I1 = table2array(I1);

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
figure(1) 
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
figure(2)
hold on 
plot(prob_exceedance(1,1:70), rank_I1, 'c')
plot(prob_exceedance(1,1:70), rank_u1, 'r')
plot(prob_exceedance(1,1:70), rank_u3, 'b')
plot(prob_exceedance(1,1:70), rank_Q3, 'g')
hold off

%test
disp('Years')
disp(N)
disp('I1')
disp(length(I1))
disp('u1')
disp(length(u1))
disp('Q1')
disp(length(Q1))
disp('R2')
disp(length(R2))
disp('Q2')
disp(length(Q2))
disp('u3')
disp(length(u3))
disp('Q3')
disp(length(Q3))

% Create Results table 
ResultsTable = table(Years, I1, u1, Q1, R2, Q2, u3, Q3, ...
    'VariableNames', {'TimeIndex', 'Inflow_I1', 'OutflowToCity_u1', ...
    'OutflowFromNode1_Q1', 'ReturnFlow_R2', 'OutflowFromNode2_Q2', ...
    'OutflowToFarms_u3', 'OutflowFromNode3_Q3'});

disp(ResultsTable);

% code Grave yard 

