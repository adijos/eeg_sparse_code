function [a, Ihat, e, a0] = sparsify(Phi, gain, I, noise_var, beta, sigma, eta)

% Dimensions
[N, M, w] = size(Phi);
T = size(I, 2);

num_iterations = 20;

lambda_N = 1/noise_var; % ??
bos = beta/sigma;

% Padding Dimensions
marg = floor(w/2);
num_valid = T - 2 * marg;
sbuf = 1:marg;
ebuf = T - marg + 1:T;
I(:, sbuf) = 0;
I(:, ebuf) = 0;

% Initialize coefficients to ?
a = 0;
for t = 1:w
  tt = t:T - (w - t);
  a = a + Phi(:, :, t)' * I(:, tt);
end

% Normalize coefficients?
for i = 1:M
  a(i, :) = a(i, :) / (gain(i) * gain(i));
end


a0 = a; % starting coefficients

% Compute initial reconstruction as convolution of bases Phi and
% coefficients a
Ihat = zeros(N, T);
for t = 1:w
  tt = t:T - (w - t);
  Ihat(:, tt) = Ihat(:, tt) + Phi(:, :, t) * a;
end

% Get the reconstruction residual, set padded regions to 0
e = I - Ihat;
e(:, sbuf) = 0;
e(:, ebuf) = 0;

% Initalize objectiv function record
E = (0.5 * lambda_N * sum(e(:).^2) + beta * sum(S(a(:)/sigma)))/num_valid;
error_record = [];

for iter = 1:num_iterations
  error_record = [error_record, E];

  % Compute gradient of a
  grada = 0;
  for t = 1:w
    tt = t:T - (w - t);
    grada = grada + lambda_N * Phi(:, :, t)' * e(:, tt);
  end
  
  % Incorporate sparse penalty
  grada = grada - bos * Sp(a/sigma);
  da = eta * grada;
  
  % Update coefficients
  a = a + eta * da;

  % Compute reconstruction using updated coefficients
  for t = 1:w
    tt = t:T - (w - t);
    Ihat(:, tt) = Ihat(:, tt) + Phi(:, :, t) * da;
  end

  % Compute residual
  e = I - Ihat;
  e(:, sbuf) = 0;
  e(:, ebuf) = 0;
 
  % Re-evaluate objective function
  Eold = E;
  E = (0.5 * lambda_N * sum(e(:).^2) + beta * sum(S(a(:)/sigma)))/num_valid;
  fprintf('\rE=%10.4f, dEpct=%6.2f%%',E,100*(E-Eold)/Eold);
  
  % Return if objective function did not decrease
  if (E - Eold) > 0
    E
    a = [];
    return
  end
  

end

figure(2)
plot(error_record)
title('Coefficient Objective Function');
%{
figure(3)
for j = 1:16;
    subplot(4, 4, j)
    plot(a(j, :));
end
suptitle('Sparse Coefficients');


Ihat(:, sbuf) = 0;
Ihat(:, ebuf) = 0;
figure(4)
for j = 1:16;
    subplot(4, 4, j)
    hold on
    
    plot(I(j, :), 'k');
    plot(Ihat(j, :), 'b-');
    
end
suptitle('Reconstructions');
%}

% figure(302)
% for j = 1:16;
%     subplot(4, 4, j)
%     hold on
%     plot(I(j, marg+5:T - marg -5) - Ihat(j, marg+5:T - marg -5), 'b');
% end
% suptitle('Residual');


fprintf('\n')