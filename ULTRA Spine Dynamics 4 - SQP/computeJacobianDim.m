function A = computeJacobianDim(f, c, dim)
% numerical partial derivative of f with respect to dimension dim,
% centered at c{:}
% more or less, d f(c{:})/d c{dim}

x = c{dim};
nX = length(x);
x_eps = 0.01;

A = zeros(size(f(c{:}),1), nX);
for i=1:nX
	x_plus = x;
	x_plus(i) = x_plus(i) + x_eps;
    c_plus = c;
    c_plus{dim} = x_plus;
    x1_plus = f(c_plus{:});
	
    x_minus = x;
	x_minus(i) = x_minus(i) - x_eps;
	c_minus = c;
    c_minus{dim} = x_minus;
    x1_minus = f(c_minus{:});
    
	A(:,i) = (x1_plus - x1_minus) / (2*x_eps);
end
