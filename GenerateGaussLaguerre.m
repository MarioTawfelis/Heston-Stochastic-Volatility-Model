function [x, w] = GenerateGaussLaguerre(n)
% Generate abscissas (x) and weights (w) for Gauss Laguerre integration

% The Laguerre polynomial
for k=0:n
	L(k+1) = (-1)^k/factorial(k)*nchoosek(n,k);
end

% Need to flip the vector to get the roots in the correct order
L = fliplr(L);

% Find the roots.  
x = flipud(roots(L));

% Find the weights
w = zeros(length(x),1);
for j=1:length(x)
	% The coefficients of the derivative of the Laguerre polynomial
	for k=1:n
		dL(k,j) = (-1)^k/factorial(k-1)*nchoosek(n,k)*x(j)^(k-1);
	end
	% The weight w(j)
	w(j) = 1/x(j)/sum(dL(:,j))^2;
	% The total weight w(j)exp(x(j))
	w(j) = w(j)*exp(x(j));
end

% Reference : The Heston Model and Its Extensions in Matlab and C# by Fabrice Douglas Rouah