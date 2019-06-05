%
% RK4.m
% function to implement fourth order Runge-Kuntta Integrator
%
% initial time, final time, initial y value and time interval as inputs
%
function sol = RK4( initial_t, final_t, initial_y, tau )

format long

f = @(t,y) (0.1*y); % given function

y(1) = initial_y; % initial value

T = final_t - initial_t; % total time
t = T/tau; % number of time steps

% Runge - Kutta - fourth order
for istep = 2:t+1
    
    k1 = tau * f((istep-1)*tau , y(istep-1)); 
    k2 = tau * f((istep-1)*tau + 2/9*tau , y(istep-1) + 2/9*k1);
    k3 = tau * f((istep-1)*tau + 1/3*tau , y(istep-1) + 1/12*k1 + 1/4*k2); 
    k4 = tau * f((istep-1)*tau + 3/4*tau , y(istep-1) + 69/128*k1 - 243/128*k2 + 135/64*k3);
    k5 = tau * f((istep-1)*tau + tau , y(istep-1) - 17/12*k1 + 27/4*k2 - 27/5*k3 + 16/15*k4);
    % solution for the individual time step
    y(istep) = y(istep-1) + 1/9*k1 + 9/20*k3 + 16/45*k4 + 1/12*k5;
    
end

%solution at final time
sol = y(t+1);
end

