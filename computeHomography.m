function H = computeHomography(x, y, xprime, yprime )
%Uses the homogenous least squares method as described on slide 13 of
%lecture 7

A = zeros(9,9);

for p = 1:4
    A(2*p - 1, :) = [x(p), y(p), 1, 0, 0, 0, -xprime(p) * x(p), -xprime(p) * y(p), -xprime(p)];
    A(2*p, :) = [0, 0, 0, x(p), y(p), 1, -yprime(p) * x(p), -yprime(p) * y(p), -yprime(p)];
end

[U, S, V] = svd(A);
H = V(:,9);
H = reshape(H, 3, 3)'