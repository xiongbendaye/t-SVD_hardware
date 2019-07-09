function A_hat = fft_tsvd(A)
[n1, n2, n3] = size(A);
A_hat = zeros(n1, n2, n3);

for i=1:n1
    for j=1:n2
        A_hat(i,j,:) = fft(squeeze(A(i,j,:)));
    end
end
end