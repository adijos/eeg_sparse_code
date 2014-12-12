beta = 10;
sigma = 2;
x = -5:0.01:5;
hold on
plot(x, exp(-(beta * log(1 + (x/sigma).^2))));
plot(x, beta * log(1 + (x/sigma).^2));
plot(x, beta * 2 * x ./ (sigma^2 + x.^2));
hold off