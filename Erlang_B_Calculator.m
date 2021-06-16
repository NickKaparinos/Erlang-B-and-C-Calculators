%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Erlang B Calculator         %
%                                         %
%              Nick Metalhead             %
%                                         %
% Uses:                                   %
% Calculation of P_B , P_n, P_>=n and     %
% the average number of users.            %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

%%% Variables %%%
% Erlangs %
rho = 3;

% Number of lines %
N = 4;

%%% Calculate blocking probablity %%%
P_B = calculate_P_B(rho, N);
disp("P_B =  " + num2str(P_B));

%%% Calculate the probability of n users in the system %%%
% Number of users %
n = 1;

P_n = calculate_P_n(rho, N, n);
disp("P_" + num2str(n) + " =  " + num2str(P_n));

%%% Calculate the probability that there are at least n users in the system %%%
% Number of users %
n = 2;

P_at_least_n = calculate_P_at_least_n(rho, N, n);
disp("P_n>=" + num2str(n) + " =  " + num2str(P_at_least_n));

%%% Average average number of users in the system %%%
avg = calculate_average(rho,N);
disp("Average number of users in the system = " + avg);



%%% Functions %%%
function P_n = calculate_P_n(rho, N, n)
    % Calculate P_n (probability of N users in the system)
    
    % Sum
    temp_sum = 0;
    for k = 0:N
        temp_sum = temp_sum + rho^k/factorial(k);
    end
    
    % Result
    P_n = rho^n/factorial(n) * 1/temp_sum;
end

function P_b = calculate_P_B(rho,N)
    % Calculate P_b (blocking probability)
    
    % Calculate probability of system being full, n=N
    P_b = calculate_P_n(rho,N,N);
end

function P_less_than_n = calculate_P_less_than_n(rho, N, n)
    % Calculate the probability than there are less than n users in the
    % system
    
    P_less_than_n = 0;
    for i = 0:n-1
        P_less_than_n = P_less_than_n + calculate_P_n(rho, N, i);
    end
end

function P_at_least_n = calculate_P_at_least_n(rho, N, n)
    % Calculate the probability than there are at least n users in the
    % system
    
    P_at_least_n = 1 - calculate_P_less_than_n(rho, N, n);
end


function avg = calculate_average(rho, N)
    % Calculate the average number of users in the system
    
    avg = 0;
    for i = 1:N
        avg = avg + i*calculate_P_n(rho, N, i);
    end
end
