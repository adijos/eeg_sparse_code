function derivative = S_prime(a,beta,sigma)
derivative = beta * 2 * a ./ (sigma^2 + a.^2);
end