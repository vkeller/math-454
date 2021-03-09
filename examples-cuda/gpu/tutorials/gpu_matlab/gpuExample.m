%%%%%%%%%%%%%%%%%%
%  gpuExample.m  %
%%%%%%%%%%%%%%%%%%
function gpuExample(A, B)
%function gpuExample(A, B)
% This function computes matrix product of A and B on Client and GPU
% Collects walltime and check if both results agree (to 5 decimal place)
% A - MATLAB Client array (N x N)
% B - MATLAB Client array (N x N)
% Usage example:
% >> N = 3000; A = rand(N); B = rand(N);
% >> gpuExample(A, B)

tic
C = A*B;    % matrix product on Client
tC = toc;
% copy A and B from Client to GPU
a = gpuArray(A); b = gpuArray(B);
tic
c = a*b;    % matrix product on GPU
tgpu = toc;
tic
CC = gather(c);   % copy data from GPU to Client
tg = toc;

disp(['Matrix multiply time on Client is ' num2str(tC)])
disp(['Matrix multiply time on GPU is ' num2str(tgpu)])
disp(['Time for gathering data from GPU back to Client is ' num2str(tg)])

% Verify that GPU and Client computations agree
tol = 1e-5;
if any(abs(CC-C) > tol)
    disp('Matrix product on Client and GPU disagree')
else
    disp('Matrix product on Client and GPU agree')
end

end   % function gpuExample
