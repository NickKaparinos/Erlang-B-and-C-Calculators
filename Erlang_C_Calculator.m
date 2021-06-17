%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Erlang C Calculator         %
%                                         %
%              Nick Kaparinos             %
%                                         %
% Uses:                                   %
% Calculation of P_0, P_T, E_T, E_N, E_M  %
% and E_W                                 %
%                                         %
% Variables:                              %
% Lambda, mu, N                           %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

%%% Variables %%%
lambda = 50;
mu = 64;
A = lambda/mu;

% Number of lines %
N = 1;

%%% Calculate P_0, probability of the system being empty %%%
P_0 = calculate_P_0(lambda, mu, N);
disp("P_0 =  " + num2str(P_0));

%%% Calculate P_T, probability of going in the queue %%%
P_T = calculate_P_T(lambda, mu, N);
disp("P_T =  " + num2str(P_T));

%%% Calculate E_W, average time of packets in queue %%%
E_W = calculate_E_W(lambda, mu, N);
disp("E_W =  " + num2str(E_W));

%%% Calculate E_T, average time of packets staying in the sytem %%%
E_T = calculate_E_T(lambda, mu, N);
disp("E_T =  " + num2str(E_T));

%%% Calculate E_N, average number of packets in the system %%%
E_N = calculate_E_N(lambda, mu, N);
disp("E_N =  " + num2str(E_N));

%%% Calculate E_M, average number of packets in queue %%%
E_m = calculate_E_m(lambda, mu, N);
disp("E_M =  " + num2str(E_m));


%%% Functions %%%
function P_0 = calculate_P_0(lambda, mu, N)
    % Calculate P_0 (probability of an empty system)
    
    % Calculate probability of the system being empty
    P_0 = 1/(1-lambda/(mu*N)) * (lambda/mu)^N/factorial(N);
    for i =0:N-1
        P_0 = P_0 + (lambda/mu)^i/factorial(i) ;
    end
    P_0 = 1/P_0;
end

function P_T = calculate_P_T(lambda, mu, N)
    % Calculate P_T (probability of entering the queue)
    
    % Calculate P_B for Erlang B
    P_B = calculate_Erlang_B_P_B(lambda/mu, N);
    
    % Calculate P_T
    P_T = N*P_B/(N-lambda/mu*(1-P_B));
end

function E_m = calculate_E_m(lambda, mu, N)
    % Calculate E_m (average number of packets in queue)
    A = lambda/mu;
    
    E_m = A*calculate_P_T(lambda, mu, N)/(N-A);

end

function E_W = calculate_E_W(lambda, mu, N)
    % Calculate E_w (average time in queue)
    
    E_W = calculate_E_m(lambda, mu, N)/lambda;
end

function E_T = calculate_E_T(lambda, mu, N)
    % Calculate E_T (average delay time)
    
    E_T = 1/mu + calculate_E_m(lambda, mu, N)/lambda;
end

function E_N = calculate_E_N(lambda, mu, N)
    % Calculate E_N (average number of packets in the system)
    rho = lambda/(N*mu);
    
    E_N = N*rho +  calculate_E_m(lambda, mu, N);
end


% Erlang B
function P_n = calculate_Erlang_B_P_n(rho, N, n)
    % Calculate P_n (probability of N users in the system)
    
    % Sum
    temp_sum = 0;
    for k = 0:N
        temp_sum = temp_sum + rho^k/factorial(k);
    end
    
    % Result
    P_n = rho^n/factorial(n) * 1/temp_sum;
end

function P_b = calculate_Erlang_B_P_B(rho,N)
    % Calculate P_b (blocking probability)
    
    % Calculate probability of system being full, n=N
    P_b = calculate_Erlang_B_P_n(rho,N,N);
end