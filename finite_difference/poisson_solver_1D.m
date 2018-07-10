function potential = poisson_solver_1D(f, x0, xn, phi_0, phi_n, n)

% Initializing x-grid and setting dirichlet boundary conditions
h = (xn - x0)/n
x_array   = x0:h:xn;
potential = zeros(n+1, 1);
potential(1) = phi_0;
potential(n + 1) = phi_n;

% Initializing coefficient matrix and functional matrix %

A = zeros(n - 1);
func_f = zeros(n - 1, 1);

% Top and bottom row values of coefficient matrix
A(1,1) = 2/power(h, 2);
A(1,2) = -1/power(h, 2);
A(n - 1,n - 1) = 2/power(h, 2);
A(n - 1,n - 2) = -1/power(h, 2);

% Top and Bottom Element of functional matrix, embodies the boundary conditions
func_f(1) = f(x_array(2)) + (phi_0/power(h, 2)); %  
func_f(n - 1) = f(x_array(n)) + (phi_n/power(h, 2)); %  

% Populating the other elements of the coefficient matrix
for i = 2:n - 2
  A(i, i) = 2/power(h, 2);
  A(i, i - 1) = -1/power(h, 2);
  A(i, i + 1) = -1/power(h, 2);
end

% Populating the functional matrix with the charge density values
for i = 2:n - 2
  func_f(i) = f(x_array(i + 1));
end

% Now we have to solve Au = f
% This will be done by Gaussian Elimination

% Creating Upper Triangular Matrix
for j = 1:n - 1
  for k = j + 1:n - 1
    lambda = A(k, j)/A(j, j);
    A(k, j:n-1) = A(k, j:n-1) - lambda*A(j,j:n-1);
    func_f(k) = func_f(k) - lambda*func_f(j);
  end
end

% Back Substitution
potential(n) = (func_f(n - 1)/A(n - 1, n - 1));
for j = n-2:-1:1
  potential(j + 1) = (1/A(j, j))*(func_f(j) - (A(j, j+1:n-1)*potential(j+2:n)));
end

% Plotting the potential obtained
plot(x_array, potential)
xlabel("x (m)")
ylabel("Potential (V)")
title("Electric Potential vs X (Dirichlet Boundary Conditions)")