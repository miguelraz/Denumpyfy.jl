using pylib::FileWriteString
using std::fs::File
using std::fs::OpenOptions
# Code to construct puncture initial data for single black hole.

using numpy: zeros, size, sqrt, linspace
import scipy.linalg
struct EllipticSolver
n_grid::
delta::
rhs_1d::
A::
sol::
rad::
end

function __init__{T0, T1, T2}(self::EllipticSolver, x::T0, y::T1, z::T2)
println(join([" Setting up Poisson solver..."], " "));
self.n_grid = size(x)
self.delta = (x[1] - x[0])
nnn = pow(self.n_grid, 3)
self.rhs_1d = zeros(nnn)
self.A = zeros((nnn, nnn))
self.sol = zeros((self.n_grid, self.n_grid, self.n_grid))
self.rad = zeros((self.n_grid, self.n_grid, self.n_grid))
for i in (0:self.n_grid - 1)
for j in (0:self.n_grid - 1)
for k in (0:self.n_grid - 1)
rad2 = ((pow(x[i], 2) + pow(y[j], 2)) + pow(z[k], 2))
self.rad[(i, j, k)] = sqrt(rad2)
end
end
end
end

function setup_matrix{T0}(self::EllipticSolver, fct::T0)
n_grid = self.n_grid
i = 0
for j in (0:n_grid - 1)
for k in (0:n_grid - 1)
index = super_index(self, convert(, i), j, k)
self.A[(index, index)] = self.rad[(i, j, k)]
self.A[(index, (index + 1))] = -(self.rad[((i + 1), j, k)])
end
end
i = (n_grid - 1);
for j in (0:n_grid - 1)
for k in (0:n_grid - 1)
index = super_index(self, convert(, i), j, k)
self.A[(index, index)] = self.rad[(i, j, k)]
self.A[(index, (index - 1))] = -(self.rad[((i - 1), j, k)])
end
end
j = 0
for i in (1:(n_grid - 1) - 1)
for k in (0:n_grid - 1)
index = super_index(self, i, convert(, j), k)
self.A[(index, index)] = self.rad[(i, j, k)]
self.A[(index, (index + n_grid))] = -(self.rad[(i, (j + 1), k)])
end
end
j = (n_grid - 1);
for i in (1:(n_grid - 1) - 1)
for k in (0:n_grid - 1)
index = super_index(self, i, convert(, j), k)
self.A[(index, index)] = self.rad[(i, j, k)]
self.A[(index, (index - n_grid))] = -(self.rad[(i, (j - 1), k)])
end
end
k = 0
for i in (1:(n_grid - 1) - 1)
for j in (1:(n_grid - 1) - 1)
index = super_index(self, i, j, convert(, k))
self.A[(index, index)] = self.rad[(i, j, k)]
self.A[(index, (index + (n_grid*n_grid)))] = -(self.rad[(i, j, (k + 1))])
end
end
k = (n_grid - 1);
for i in (1:(n_grid - 1) - 1)
for j in (1:(n_grid - 1) - 1)
index = super_index(self, i, j, convert(, k))
self.A[(index, index)] = self.rad[(i, j, k)]
self.A[(index, (index - (n_grid*n_grid)))] = -(self.rad[(i, j, (k - 1))])
end
end
for i in (1:(n_grid - 1) - 1)
for j in (1:(n_grid - 1) - 1)
for k in (1:(n_grid - 1) - 1)
index = super_index(self, i, j, k)
self.A[(index, index)] = (-6.0 + (pow(self.delta, 2)*fct[(i, j, k)]))
self.A[(index, (index - 1))] = 1.0
self.A[(index, (index + 1))] = 1.0
self.A[(index, (index - n_grid))] = 1.0
self.A[(index, (index + n_grid))] = 1.0
self.A[(index, (index - (n_grid*n_grid)))] = 1.0
self.A[(index, (index + (n_grid*n_grid)))] = 1.0
end
end
end
end

function setup_rhs{T0}(self::EllipticSolver, rhs::T0)
n_grid = self.n_grid
for i in (1:(n_grid - 1) - 1)
for j in (1:(n_grid - 1) - 1)
for k in (1:(n_grid - 1) - 1)
index = super_index(self, i, j, k)
self.rhs_1d[index] = (pow(self.delta, 2)*rhs[(i, j, k)])
end
end
end
end

function solve{RT}(self::EllipticSolver)::RT
sol_1d = solve(la)
for i in (0:self.n_grid - 1)
for j in (0:self.n_grid - 1)
for k in (0:self.n_grid - 1)
index = super_index(self, i, j, k)
self.sol[(i, j, k)] = sol_1d[index]
end
end
end
return self.sol
end

function super_index{T0, T1, T2, RT}(self::EllipticSolver, i::T0, j::T1, k::T2)::RT
return (i + (self.n_grid*(j + (self.n_grid*k))))
end

struct Puncture
bh_loc::
lin_mom::
n_grid::
x_out::
delta::Float64
x::
y::
z::
solver::EllipticSolver
alpha::
beta::
u::
res::
end

function __init__{T0, T1, T2, T3}(self::Puncture, bh_loc::T0, lin_mom::T1, n_grid::T2, x_out::T3)
self.bh_loc = bh_loc
self.lin_mom = lin_mom
println(join([" Constructing class Puncture for single black hole"], " "));
println(join(["    at bh_loc = (", bh_loc[0], ",", bh_loc[1], ",", bh_loc[2], ")"], " "));
println(join(["    with momentum p = (", lin_mom[0], ",", lin_mom[1], ",", lin_mom[2], ")"], " "));
println(join([" Using", n_grid, "^3 gridpoints with outer boundary at", x_out], " "));
self.n_grid = n_grid
self.x_out = x_out
self.delta = ((2.0*x_out)/n_grid)
half_delta = (self.delta/2.0)
self.x = linspace((half_delta - x_out), (x_out - half_delta), n_grid)
self.y = linspace((half_delta - x_out), (x_out - half_delta), n_grid)
self.z = linspace((half_delta - x_out), (x_out - half_delta), n_grid)
self.solver = EllipticSolver(self.x, self.y, self.z)
self.alpha = zeros((n_grid, n_grid, n_grid))
self.beta = zeros((n_grid, n_grid, n_grid))
self.u = zeros((n_grid, n_grid, n_grid))
self.res = zeros((n_grid, n_grid, n_grid))
end

function construct_solution{T0, T1}(self::Puncture, tol::T0, it_max::T1)
setup_alpha_beta(self);
residual_norm = residual(self)
println(join([" Initial Residual = ", residual_norm], " "));
println(join([" Using up to", it_max, "iteration steps to reach tolerance of", tol], " "));
it_step = 0
while residual_norm > tol&&it_step < it_max
it_step += 1
update_u(self);
residual_norm = residual(self);
println(join([" Residual after", it_step, "iterations :", residual_norm], " "));
end
if residual_norm < tol
println(join([" Done!"], " "));
else

println(join([" Giving up..."], " "));
end
end

function update_u(self::Puncture)
n_grid = self.n_grid
fct = zeros((n_grid, n_grid, n_grid))
rhs = zeros((n_grid, n_grid, n_grid))
for i in (1:(n_grid - 1) - 1)
for j in (1:(n_grid - 1) - 1)
for k in (1:(n_grid - 1) - 1)
temp = ((self.alpha[(i, j, k)]*(1.0 + self.u[(i, j, k)])) + 1.0)
fct[(i, j, k)] = (((-7.0*self.beta[(i, j, k)])*self.alpha[(i, j, k)])/pow(temp, 8))
rhs[(i, j, k)] = -(self.res[(i, j, k)])
end
end
end
setup_matrix(self.solver, fct);
setup_rhs(self.solver, rhs);
delta_u = solve(self.solver)
self.u += delta_u
end

function residual(self::Puncture)::Float64
residual_norm = 0.0
for i in (1:(self.n_grid - 1) - 1)
for j in (1:(self.n_grid - 1) - 1)
for k in (1:(self.n_grid - 1) - 1)
ddx = ((self.u[((i + 1), j, k)] - (2.0*self.u[(i, j, k)])) + self.u[((i - 1), j, k)])
ddy = ((self.u[(i, (j + 1), k)] - (2.0*self.u[(i, j, k)])) + self.u[(i, (j - 1), k)])
ddz = ((self.u[(i, j, (k + 1))] - (2.0*self.u[(i, j, k)])) + self.u[(i, j, (k - 1))])
lhs = (((ddx + ddy) + ddz)/pow(self.delta, 2))
temp = ((self.alpha[(i, j, k)]*(1.0 + self.u[(i, j, k)])) + 1.0)
rhs = (-(self.beta[(i, j, k)])/pow(temp, 7))
self.res[(i, j, k)] = (lhs - rhs)
residual_norm += pow(self.res[(i, j, k)], 2)
end
end
end
residual_norm = (sqrt(residual_norm)*pow(self.delta, 3));
return residual_norm
end

function setup_alpha_beta(self::Puncture)
n_grid = self.n_grid
p_x = self.lin_mom[0]
p_y = self.lin_mom[1]
p_z = self.lin_mom[2]
for i in (0:n_grid - 1)
for j in (0:n_grid - 1)
for k in (0:n_grid - 1)
s_x = (self.x[i] - self.bh_loc[0])
s_y = (self.y[j] - self.bh_loc[1])
s_z = (self.z[k] - self.bh_loc[2])
s2 = ((pow(s_x, 2) + pow(s_y, 2)) + pow(s_z, 2))
s_bh = sqrt(s2)
l_x = (s_x/s_bh)
l_y = (s_y/s_bh)
l_z = (s_z/s_bh)
lP = (((l_x*p_x) + (l_y*p_y)) + (l_z*p_z))
fac = (3.0/(2.0*s2))
A_xx = (fac*(((2.0*p_x)*l_x) - ((1.0 - (l_x*l_x))*lP)))
A_yy = (fac*(((2.0*p_y)*l_y) - ((1.0 - (l_y*l_y))*lP)))
A_zz = (fac*(((2.0*p_z)*l_z) - ((1.0 - (l_z*l_z))*lP)))
A_xy = (fac*(((p_x*l_y) + (p_y*l_x)) + ((l_x*l_y)*lP)))
A_xz = (fac*(((p_x*l_z) + (p_z*l_x)) + ((l_x*l_z)*lP)))
A_yz = (fac*(((p_y*l_z) + (p_z*l_y)) + ((l_y*l_z)*lP)))
A2 = (((pow(A_xx, 2) + pow(A_yy, 2)) + pow(A_zz, 2)) + (2.0*((pow(A_xy, 2) + pow(A_xz, 2)) + pow(A_yz, 2))))
self.alpha[(i, j, k)] = (2.0*s_bh)
self.beta[(i, j, k)] = ((pow(self.alpha[(i, j, k)], 7)*A2)/8.0)
end
end
end
end

function write_to_file(self::Puncture)
filename = ((("Puncture_" + string(self.n_grid)) + "_") + string(self.x_out))
filename = (filename + ".data");
out = OpenOptions::new().write(true).open(filename)
if out
k = (self.n_grid / 2)
write(out, ("# Data for black hole at x = (%f,%f,%f)
" % (self.bh_loc[0], self.bh_loc[1], self.bh_loc[2])));
write(out, ("# with linear momentum P = (%f, %f, %f)
" % self.lin_mom));
write(out, ("# in plane for z = %e 
" % self.z[k]));
write(out, "# x            y              u              
");
write(out, "#============================================
");
for i in (0:self.n_grid - 1)
for j in (0:self.n_grid - 1)
write(out, ("%e  %e  %e
" % (self.x[i], self.y[j], self.u[(i, j, k)])));
end
write(out, "
");
end
close(out);
else

println(join([" Could not open file", filename, "in write_to_file()"], " "));
println(join([" Check permissions?"], " "));
end
end

function main_func()
println(join([" -------------------------------------------------------"], " "));
println(join([" --- puncture.py --- use flag -h for list of options ---"], " "));
println(join([" -------------------------------------------------------"], " "));
loc_x = 0.0
loc_y = 0.0
loc_z = 0.0
p_x = 1.0
p_y = 0.0
p_z = 0.0
n_grid = 16
x_out = 4.0
tol = 1e-12
it_max = 50
for i in (0:length(append!([PROGRAM_FILE], ARGS)) - 1)
if append!([PROGRAM_FILE], ARGS)[i] == "-h"
usage();
return
end
if append!([PROGRAM_FILE], ARGS)[i] == "-n_grid"
n_grid = Int64(append!([PROGRAM_FILE], ARGS)[(i + 1)]);
end
if append!([PROGRAM_FILE], ARGS)[i] == "-x_out"
x_out = float(append!([PROGRAM_FILE], ARGS)[(i + 1)]);
end
if append!([PROGRAM_FILE], ARGS)[i] == "-loc_x"
loc_x = float(append!([PROGRAM_FILE], ARGS)[(i + 1)]);
end
if append!([PROGRAM_FILE], ARGS)[i] == "-loc_y"
loc_y = float(append!([PROGRAM_FILE], ARGS)[(i + 1)]);
end
if append!([PROGRAM_FILE], ARGS)[i] == "-loc_z"
loc_z = float(append!([PROGRAM_FILE], ARGS)[(i + 1)]);
end
if append!([PROGRAM_FILE], ARGS)[i] == "-p_x"
p_x = float(append!([PROGRAM_FILE], ARGS)[(i + 1)]);
end
if append!([PROGRAM_FILE], ARGS)[i] == "-p_y"
p_y = float(append!([PROGRAM_FILE], ARGS)[(i + 1)]);
end
if append!([PROGRAM_FILE], ARGS)[i] == "-p_z"
p_z = float(append!([PROGRAM_FILE], ARGS)[(i + 1)]);
end
if append!([PROGRAM_FILE], ARGS)[i] == "-tol"
tol = float(append!([PROGRAM_FILE], ARGS)[(i + 1)]);
end
if append!([PROGRAM_FILE], ARGS)[i] == "-it_max"
it_max = Int64(append!([PROGRAM_FILE], ARGS)[(i + 1)]);
end
end
bh_loc = (loc_x, loc_y, loc_z)
lin_mom = (p_x, p_y, p_z)
black_hole = Puncture(bh_loc, lin_mom, n_grid, x_out)
construct_solution(black_hole, tol, it_max);
write_to_file(black_hole);
end

function usage()
println(join(["Constructs puncture initial data for single black hole."], " "));
println(join([""], " "));
println(join(["The following options can be used to over-write default parameters"], " "));
println(join(["	-n_grid: number of grid points [default: 16]"], " "));
println(join(["	-x_out: location of outer boundary [4.0]"], " "));
println(join(["	-loc_x, -loc_y, -loc_z: location of black hole [(0.0, 0.0, 0.0)]"], " "));
println(join(["	-p_x, -p_y, -p_z: lin. momentum of black hole [(1.0, 0.0, 0.0)]"], " "));
println(join(["	-tol: tolerance for elliptic solver [1.e-12]"], " "));
println(join(["	-it_max: maximum number of iterations [50]"], " "));
println(join(["For example, to construct data with x_out = 6.0, call"], " "));
println(join(["	python3 puncture.py -x_out 6.0"], " "));
end

function main()
main_func();
end

main()
