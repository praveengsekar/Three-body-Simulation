% 
% f.m
%
% program to calculate the differential of y
% 
%

function ydot = f ( t, y )

  global m1
  global m2
  global m3

  x1 = y(1:3);
  x2 = y(7:9);
  x3 = y(13:15);

  d1 = ( x3 - x2 ) / norm ( x3 - x2 )^3;
  d2 = ( x1 - x3 ) / norm ( x1 - x3 )^3;
  d3 = ( x2 - x1 ) / norm ( x2 - x1 )^3;

  ydot(1:3) = y(4:6); 
  ydot(4:6) = m2 * d3 - m3 * d2;
  ydot(7:9) = y(10:12);
  ydot(10:12) = m3 * d1 - m1 * d3;
  ydot(13:15) = y(16:18);
  ydot(16:18) = m1 * d2 - m2 * d1;

  return
end