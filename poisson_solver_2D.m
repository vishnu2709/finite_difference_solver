% f is the 2D charge density function
% top, left, bottom, right are the boundary condition function
function potential = poisson_solver_2D(f, top, left, bottom, right, xo, xn, yo, yn, h)
% Getting no of points from input values specified
  try
   nx = typecast((xn - xo)/h, int16)
   ny = typecast((yn - yo)/h, int16)
  catch
   disp 'The choice of h is such that the no of grid points is not coming out be an integer. 
          Please choose another value of h such that the no of points is a valid12 int.'
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
  potential(nx+1, i) = top(grid(1, 1, nx+1, i), grid(1, 2, nx+1, i));
end

% Left and Right Boundary 
for j = 2:ny
  potential(j, 1) = left(grid(1, 1,j, 0), grid(1, 2, j, 0));
  potential(j, nx+1) = right(grid(1, 1, j, nx+1), grid(1, 2, j, nx+1));

% Initializing coefficient matrix
coeff_matrix = zeros(ny - 1, ny - 1, nx - 1, nx - 1);
% Initializing diagonal element of coefficient matrix
A = zeros(ny - 1);
A(1, 1) = 4; A(1, 2) = -1; A(ny - 1, ny - 1) = 4; A(ny - 1, ny - 2) = -1

for i = 2:ny - 2
  A(i, i) = 4
  A(i, i - 1) = -1
  A(i, i + 1) = -1
end

coeff_matrix(:,:,1,1) = A; coeff_matrix(:,:,1,2) = -eye(ny - 1)

