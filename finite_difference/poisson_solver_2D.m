% f is the 2D charge density function
% top, left, bottom, right are the boundary condition function
function potential = poisson_solver_2D(f, top, left, bottom, right, xo, xn, yo, yn, h)
% Getting no of points from input values specified
nx = (xn - xo)/h
ny = (yn - yo)/h

if (floor(nx) != nx | floor(ny) != ny)
   disp("The choice of h is such that the no of grid points is not coming out to \
   be an integer. Please choose another value of h \
   such that the no of points is a valid int.")
   return
end

% Initializing Grid
grid = zeros(1, 2, ny + 1, nx + 1);
potential = zeros(ny + 1, nx + 1);

for i = 1: nx + 1
  for j = 1:ny + 1
    grid(:,:,j, i) = [xo + h*(i - 1), yn - h*(j - 1)];
  end
end

% Applying Boundary Conditions
% Top and Bottom Boundary
for i = 1:nx+1
  potential(1, i) = top(grid(1, 1, 1, i), grid(1, 2, 1, i));
  potential(ny+1, i) = bottom(grid(1, 1, ny+1, i), grid(1, 2, ny+1, i));
end

% Left and Right Boundary 
for j = 2:ny
  potential(j, 1) = left(grid(1, 1,j, 1), grid(1, 2, j, 1));
  potential(j, nx+1) = right(grid(1, 1, j, nx+1), grid(1, 2, j, nx+1));
end

% Initializing diagonal element of coefficient matrix
A = zeros(ny - 1);
A(1, 1) = 4; A(1, 2) = -1; A(ny - 1, ny - 1) = 4; A(ny - 1, ny - 2) = -1;

for i = 2:ny - 2
  A(i, i) = 4;
  A(i, i - 1) = -1;
  A(i, i + 1) = -1;
end

% Initializing coefficient matrix
coeff_matrix = zeros(ny - 1, ny - 1, nx - 1, nx - 1);

coeff_matrix(:,:,1,1) = A; coeff_matrix(:,:,1,2) = -eye(ny - 1);
coeff_matrix(:,:,nx-1,nx-1) = A; coeff_matrix(:,:,nx-1, nx-2) = -eye(ny - 1);

for j = 2:nx - 2
  coeff_matrix(:,:,j,j) = A;
  coeff_matrix(:,:,j,j-1) = -eye(ny - 1);
  coeff_matrix(:,:,j,j+1) = -eye(ny - 1);
end

% Initializing functional matrix
func_matrix = zeros(ny - 1, 1, nx - 1, 1);

% Populating the matrix with functional values and the top bottom boundary conditions
for i = 2:nx
  func = zeros(ny - 1, 1);
  top_bottom = zeros(ny - 1, 1);
  for j = 1:ny-1
    func(j) = power(h, 2)* f(grid(1, 1, j + 1, i), grid(1, 2, j + 1, i));
  end
  top_bottom(1) = potential(1, i);
  top_bottom(ny - 1) = potential(ny+1, i);
  func_matrix(:,:,i - 1) = func + top_bottom;
end

% Adding the left and right boudnary conditions to the functional matrix
left_boundary = zeros(ny - 1, 1);
right_boundary = zeros(ny - 1, 1);

for k = 1:ny - 1
  left_boundary(k) = potential(k + 1, 1);
  right_boundary(k) = potential(k + 1, nx + 1);
end

func_matrix(:,:,1) = func_matrix(:,:,1) + left_boundary;
func_matrix(:,:,nx - 1) = func_matrix(:,:,nx - 1) + right_boundary;
% Required matrices have been set up. Now we have to do Gaussian Elimination
% (Coeff matrix)U = (func_matrix) %

for i = 1:nx - 1
  for j = i + 1:nx - 1
    lambda = coeff_matrix(:,:,j,i)/coeff_matrix(:,:,i,i);
    for k = i:nx - 1
      coeff_matrix(:,:,j,k) = coeff_matrix(:,:,j,k) - lambda*coeff_matrix(:,:,i,k);
    end
    func_matrix(:,:,j) = func_matrix(:,:,j) - lambda*func_matrix(:,:,i);
  end
end

% Back Substitution
potential(2:ny, nx) = (eye(ny - 1)/coeff_matrix(:,:,nx-1, nx-1))*func_matrix(:,:,nx - 1);
for j = nx - 2:-1:1
  sum = zeros(ny - 1, 1);
  for k = j+1:nx-1
    sum = sum + coeff_matrix(:,:,j,k)*potential(2:ny, k + 1);
  end
  potential(2:ny,j+1) = (eye(ny - 1)/coeff_matrix(:,:,j,j))*(func_matrix(:,:,j) - sum);
end

% Making a surface plot of the potential obtained
x_array = grid(1,1,:,:);
y_array = grid(1,2,:,:);
x_array = reshape(x_array, [ny + 1, nx + 1]);
y_array = reshape(y_array, [ny + 1, nx + 1]);
surf(x_array, y_array, potential)
xlabel("x (m)")
ylabel("y (m)")
zlabel("Potential (V)")
title("Potential Values in the 2D grid")
end
