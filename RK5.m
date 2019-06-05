%
% RK5.m
% function to implement fifth order Runge-Kuntta Integrator
%
% initial time, final time, initial y value and time interval as inputs
%
function sol = RK5( initial_t, final_t, initial_y, tau )

format long
f = @(t,y) (0.1*y); % given function

y(1) = initial_y; % initial value

T = final_t - initial_t; % total time
 t = T/tau; % number of time steps

% Runge - Kutta - fifth order
for istep = 2:t+1
    
    k1 = tau * f((istep-1)*tau , y(istep-1)); 
    k2 = tau * f((istep-1)*tau + 2/9*tau , y(istep-1) + 2/9*k1);
    k3 = tau * f((istep-1)*tau + 1/3*tau , y(istep-1) + 1/12*k1 + 1/4*k2); 
    k4 = tau * f((istep-1)*tau + 3/4*tau , y(istep-1) + 69/128*k1 - 243/128*k2 + 135/64*k3);
    k5 = tau * f((istep-1)*tau + tau , y(istep-1) - 17/12*k1 + 27/4*k2 - 27/5*k3 + 16/15*k4);
    k6 = tau * f((istep-1)*tau + 5/6*tau , y(istep-1) + 65/432*k1 - 5/16*k2 + 13/16*k3 + 4/27*k4 + 5/144*k5);
    % solution for the individual time step
    y(istep) = y(istep-1) + 47/450*k1 + 12/25*k3 + 32/225*k4 + 1/30*k5 + 6/25*k6; 
end

%solution at final time
sol = y(t+1);

end
