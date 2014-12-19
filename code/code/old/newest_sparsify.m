% fita.m - function to fit coefficients to movie (gradient descent)
%
% [a, Ihat, e, a0] = fita(Phi,gain,I,noise_var,beta,sigma,eta)
%
% returns 'a' empty if solution diverges

function [a, Ihat, e, a0] = newest_sparsify(Phi,gain,I,noise_var,beta,sigma,eta)

[N, M, w] = size(Phi);

T = size(I, 2);

num_iterations = 5000;

lambda_N = 1/noise_var;
bos = beta/sigma;

marg = floor(w/2);
num_valid = T - 2 * marg;
sbuf = 1:marg;
ebuf = T - marg + 1:T;

I(:, sbuf) = 0;
I(:, ebuf) = 0;

a = 0;
for t = 1:w
  tt = t:T - (w - t);
  a = a + Phi(:, :, t)' * I(:, tt);
end

for i = 1:M
  a(i, :) = a(i, :) / (gain(i) * gain(i));
end

a0 = a;

Ihat = zeros(N, T);
for t = 1:w
  tt = t:T - (w - t);
  Ihat(:, tt) = Ihat(:, tt) + Phi(:, :, t) * a;
end

e = I - Ihat;
e(:, sbuf) = 0;
e(:, ebuf) = 0;

E = (0.5 * lambda_N * sum(e(:).^2) + beta * sum(S(a(:)/sigma)))/num_valid;
error_record = [];
for iter = 0:num_iterations - 1
  error_record = [error_record, E];

  grada = 0;
  for t = 1:w
    tt = t:T - (w - t);
    grada = grada + lambda_N * Phi(:, :, t)' * e(:, tt);
  end
  grada = grada - bos * Sp(a/sigma);
  da = eta * grada;
  
  a = a + da;

  for t = 1:w
    tt = t:T - (w - t);
    Ihat(:, tt) = Ihat(:, tt) + Phi(:, :, t) * da;
  end

  e = I - Ihat;
  e(:, sbuf) = 0;
  e(:, ebuf) = 0;
 
  Eold = E;
  E = (0.5 * lambda_N * sum(e(:).^2) + beta * sum(S(a(:)/sigma)))/num_valid;
  fprintf('\rE=%10.4f, dEpct=%6.2f%%',E,100*(E-Eold)/Eold);
  if (E - Eold) > 0
    E
    a = [];
    return
  end
  

end
figure(200)
plot(error_record)
suptitle('Coefficient Reconstruction Error');

figure(300)
for j = 1:16;
    subplot(4, 4, j)
    plot(a(j, :));
end
suptitle('Sparse Coefficients');

figure(301)
for j = 1:16;
    subplot(4, 4, j)
    hold on
    plot(I(j, :), 'b');
    plot(Ihat(j, :), 'g');
end
suptitle('Reconstructions');

figure(302)
for j = 1:16;
    subplot(4, 4, j)
    hold on
    plot(I(j, marg+5:T - marg -5) - Ihat(j, marg+5:T - marg -5), 'b');
end
suptitle('Residual');


fprintf('\n')
