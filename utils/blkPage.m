% Construct a block Page matrix from a stacked trajectory.
function P = blkPage(w, L)
    [T,q] = size(w);
    P = [];
    for i=1:L:size(w,1)
        if i+L > T
            break;
        end
        w_i = reshape(w(i:i+L-1,:)', q*L, 1);


        P = [P, w_i];
    end
end
