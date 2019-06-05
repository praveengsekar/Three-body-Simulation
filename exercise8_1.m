




























%
% exercise8_1.m
%
% Program to implement Runge Kutta Methods

format long

% time steps
tau = input ('Enter time interval : ');

% value at which the function has to be evaluated
T = input('Enter the value of t at which you want to evaluate y : ');

% desired accuracy
accuracy = input('Enter the required accuracy : ');

% initial time t and initial value y(0)
t=0;
y=1;

%function call
fprintf('\nThe solution using RK4 is %ld', RK4( t, T, y, tau ));
fprintf('\nThe solution using RK5 is %ld', RK5( t, T, y, tau ));
fprintf('\nThe solution using RK4 using adaptive steps is %ld\n', RK4adaptive( t, T, y, tau, accuracy));



    

