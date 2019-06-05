%
% RK4adaptive.m
% function to implement fouth order Runge-Kuntta Integrator
% with adaptive step size
%
% initial time-t, final time-T, initial y value and time interval and accuracy as inputs
%
function yrk4 = RK4adaptive( t, T, y, tau, accuracy)

format long

f = @(t,y) (0.1*y); % given function
% Runge - Kutta - fourth order with adaptive step size
while (t <= T)
    k1 = tau * f(t , y); 
    k2 = tau * f(t + 2/9*tau , y + 2/9*k1);
    k3 = tau * f(t + 1/3*tau , y + 1/12*k1 + 1/4*k2); 
    k4 = tau * f(t + 3/4*tau , y + 69/128*k1 - 243/128*k2 + 135/64*k3);
    k5 = tau * f(t + tau , y - 17/12*k1 + 27/4*k2 - 27/5*k3 + 16/15*k4);
    yrk4 = y + 1/9*k1 + 9/20*k3 + 16/45*k4 + 1/12*k5; 
    k6 = tau * f(t + 5/6*tau , y + 65/432*k1 - 5/16*k2 + 13/16*k3 + 4/27*k4 + 5/144*k5);
    yrk5 = y + 47/450*k1 + 12/25*k3 + 32/225*k4 + 1/30*k5 + 6/25*k6; 
    
    % Tuning the step-size to achieve desired accuracy
    delta = abs(yrk4 - yrk5);
    tau = tau * ( (accuracy/delta)^(1/5) ) * 0.9;
        
    if ( delta < accuracy)
        t = t + tau;
        y = yrk4;
    end
           
end

end