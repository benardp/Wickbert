function f = ftest(x);
f(1) = x(1) ^ 2 - 4 * x(2);
f(2) = x(2) ^ 2 - 2 * x(1) + 4 * x(2);
f = f';