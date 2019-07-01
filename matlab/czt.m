function x_hat = czt(x)
n = max(size(x));
m = n;
x_hat = zeros(m, 1);
w = exp(-2j * pi / m);
a = 1;
chp = w .^ (((1-n):(max(m,n)-1)) .^ 2 / 2.0);
N2 = int32(2 ^ ceil(log2(m + n -1)));
xp = zeros(N2, 1);
for i=1:n
    xp(i) = x(i) * a ^ (-i) * chp(n-1+i);
end
ichp = zeros(N2, 1);
for i=1:(m+n-1)
    ichp(i) = 1 / chp(i);
end
r = ifft(fft(xp) .* fft(ichp));
for i=1:m
    x_hat(i) = r(n-1+i) * chp(n-1+i);
end