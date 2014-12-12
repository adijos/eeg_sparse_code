function [recon] = reconstruct_bruno(Phi, a)

    [N, M, w] = size(Phi);
    T = size(a, 3);
    
    recon = zeros(N, T);
    for t = 1:w
      tt = t:T - (w - t);
      recon(:, tt) = recon(:, tt) + Phi(:, :, t) * a;
    end

end