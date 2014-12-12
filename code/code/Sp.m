% Sp.m - derivative of sparse cost

function Sprime = Sp(a)

Sprime = 2*a./(1+a.*a);
