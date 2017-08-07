function [xprime, yprime] = applyHomography(H, x, y)

%create row vectors
[num_rowsx,num_columnsx] = size(x);
xprime = zeros(1, num_columnsx);

[num_rowsy,num_columnsy] = size(y);
yprime = zeros(1, num_columnsy);

for p = 1:num_columnsx

   new_points = inv(H) * [x(p), y(p), 1.0]';
   
   %normalize
   xprime(p) = new_points(1) / new_points(3);
   yprime(p) = new_points(2) / new_points(3);

   
end