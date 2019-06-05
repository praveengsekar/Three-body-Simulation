% 
% exercise8_2.m
%
% program to calculate the motion of three bodies
% using fourth order Runge-Kutta with adaptive-step sizes
%

global m1 m2 m3

% given masses.
  m1 = 5.0;
  m2 = 3.0;
  m3 = 4.0;

% time rang
  t_start = 0.0;
  t_stop = 10.0;
  tau = 0.01;
  tau1(1) = tau;
  t(1) = t_start; 
  
 % accuracy
  accuracy = 1e-8;
   
 % initial conditions 
 % posx(1),posy(1),posz(1),velx(1),vely(1),velz(1),posx(2),posy(2),posz(2),velx(2),vely(2),velz(2),posx(3),posy(3),posz(3),velx(3),vely(3),velz(3) 
  y = [ 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 0.0, -2.0, -1.0, 0.0, 0.0, 0.0, 0.0];
  y1(1,:) = y(1:18);
  
 % counter variable to be used as index 
  c = 2;
  
 % Runge Kutta Integrator
  while(t_start < 10)
    k1 = tau * f( t_start , y(1:18)); 
    k2 = tau * f( t_start + 2/9*tau , y(1:18) + 2/9*k1(1:18));
    k3 = tau * f( t_start + 1/3*tau , y(1:18) + 1/12*k1(1:18) + 1/4*k2(1:18)); 
    k4 = tau * f( t_start + 3/4*tau , y(1:18) + 69/128*k1(1:18) - 243/128*k2(1:18) + 135/64*k3(1:18));
    k5 = tau * f( t_start + tau , y(1:18) - 17/12*k1(1:18) + 27/4*k2(1:18) - 27/5*k3(1:18) + 16/15*k4(1:18));
    yrk4 = y(1:18) + 1/9*k1(1:18) + 9/20*k3(1:18) + 16/45*k4(1:18) + 1/12*k5(1:18);
    
    k6 = tau * f(t_start + 5/6*tau , y(1:18) + 65/432*k1(1:18) - 5/16*k2(1:18) + 13/16*k3(1:18) + 4/27*k4(1:18) + 5/144*k5(1:18));
    yrk5 = y(1:18) + 47/450*k1(1:18) + 12/25*k3(1:18) + 32/225*k4(1:18) + 1/30*k5(1:18) + 6/25*k6(1:18); 
    
    % Tuning the step-size to achieve desired accuracy
    delta = max ( abs(yrk4-yrk5) );
    tau = tau * ((accuracy/delta)^(1/5)) * 0.9;
    tau1(c) = tau;
    
    if ( delta < accuracy)
        t_start = t_start + tau;
        t(c) = t_start;
        y = yrk4(1:18);
        y1(c, :) = y(1:18); % storing the solution vector in a matrix
        c = c+1;
    end
    
     
  end

% Plotting the trajectories
plot ( y1(:,1) , y1(:,2),'b.');
hold on
plot (y1(:,7) , y1(:,8),'r.');
hold on
plot (y1(:,13) , y1(:,14),'g.');
xlabel('x');
ylabel('y');
title('Trajectories of the three bodies');
legend('1','2','3');

%plot step-size as funtion of t
figure
plot(t, tau1);
xlabel('time');
ylabel('step size');
title('Step Size vs Time');


for i = 1 : c-1
    distance_12 = sqrt( (y1(i,1)- y1(i,7))^2 + (y1(i,2)-y1(i,8))^2 );
    distance_23 = sqrt( (y1(i,7)- y1(i,13))^2 + (y1(i,8)-y1(i,14))^2 ); 
    distance_31 = sqrt( (y1(i,13)- y1(i,1))^2 + (y1(i,14)-y1(i,2))^2 );
    min_distance(i) = min (distance_12,distance_23);
    min_distance(i) = min (min_distance(i), distance_31);
end

%plot minimum distance as funtion of t
figure
plot(t, min_distance);
xlabel('time');
ylabel('minimum distance');
title('Minimal distance between any of the three bodies vs time');

    