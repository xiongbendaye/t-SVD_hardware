%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [U,S,V] = tsvd(A)
%
% Inputs:   A - 3 order tensor
%
%
% Output:   U,S,V - where U*S*Vt = A.
%
%
% Original author :  Misha Kilmer, Ning Hao
% Edited by       :  Miao Yin, Rutgers University, 2019/6/23
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,S,V] = tsvd(A, parOP)

% determine size of tensor
[n1,n2,n3]=size(A);

if ~exist('parOP','var')
    parOP = false;
end

%% Perform FFT along 3 to P axeses of Tensor

% conjugate symetric trick.
if n2 > n1
    transflag=1;
    A=tran(A);
    nn1=n1;
    n1=n2;
    n2=nn1;
end
U = zeros(n1,n1,n3);
S = zeros(n1,n2,n3);
V = zeros(n2,n2,n3);

fm = quantizer([16,6]);

A_fix = num2bin(fm, A);
A_hat = fft_tsvd(A);
A_hat_fix = num2bin(fm, A_hat);

% Do the conjugate symetric trick here.

endValue = int16(n3/2 + 1);        
[U,S,V] = takeSVDs(U,S,V,A_hat,endValue,parOP);

for j =n3:-1:endValue+1
    U(:,:,j) = conj(U(:,:,n3-j+2));
    V(:,:,j) = conj(V(:,:,n3-j+2));
    S(:,:,j) = S(:,:,n3-j+2);
end

%%
U_fix = num2bin(fm, U);
S_fix = num2bin(fm, S);
V_fix = num2bin(fm, V);

U = ifft(U,[],3);
S = ifft(S,[],3);
V = ifft(V,[],3);

U_fix = num2bin(fm, U);
S_fix = num2bin(fm, S);
V_fix = num2bin(fm, V);

if exist('transflag','var')
    Uold =U; U=V; S=tran(S); V=Uold;  
end


end
%% BEGIN SUBFUNCTIONS
%
%%
function [U, S, V] = takeSVDs(U,S,V,A,endI,runPar)

if ~exist('runPar','var')
    runPar = false;
end

if ~runPar || matlabpool('size') == 0

    for i=1:endI
        [U1,S1,V1]=svdj(A(:,:,i));
        U(:,:,i)=U1; S(:,:,i)=S1; V(:,:,i)=V1;
    end
else
    
    parfor i=1:endI
        [U1,S1,V1]=svdj(A(:,:,i));
        U(:,:,i)=U1; S(:,:,i)=S1; V(:,:,i)=V1;
    end
end

end
