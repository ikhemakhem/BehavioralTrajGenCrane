% Construct a block Hankel matrix from trajectory data.
%
% 2-D input (q x T matrix or T x q matrix, auto-transposed):
%   H = blkhank(w, i)      builds a (q*i x T-i+1) Hankel matrix
%   H = blkhank(w, i, j)   builds a (q*i x j)     Hankel matrix
%
% 3-D input (q x N x T array):
%   H = blkhank(w, i)      builds a (q*i x N*(T-i+1)) block Hankel matrix

function H = blkhank(w, i, j)
    if ndims(w) == 3
        [q, N, T] = size(w);
        if nargin < 3 || isempty(j), j = T - i + 1; end
        if j <= 0, error('Not enough data.'); end
        H = zeros(i * q, j * N);
        for ii = 1:i
            for jj = 1:j
                H(((ii-1)*q+1):(ii*q), ((jj-1)*N+1):(jj*N)) = w(:, :, ii+jj-1);
            end
        end
    else
        [q, T] = size(w);
        if T < q, w = w'; [q, T] = size(w); end
        if nargin < 3 || isempty(j), j = T - i + 1; end
        if j <= 0, error('Not enough data.'); end
        H = zeros(i * q, j);
        for ii = 1:i
            H(((ii-1)*q+1):(ii*q), :) = w(:, ii:(ii+j-1));
        end
    end
end
