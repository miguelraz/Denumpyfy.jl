using Denumpyfy
using Test

@testset "Denumpyfy.jl" begin
    str = "x ** 2"
    @test replace_exponents(str)  == "x ^ 2"

    str = "range(0, self.n_grid):"
    @test replace_range(str) == "1:self.n_grid"
    str = "range(0, 1):"
    @test replace_range(str) == "1:2"
    str = "range(self.n_grid, 3):"
    @test replace_range(str) == "self.n_grid:4"
    str = "range(1, n_grid - 1):"
    @test replace_range(str) == "2:n_grid"
    str = "    for i in range(1, n_grid - 1):"
    @test replace_range(str) == "    for i in 2:n_grid"
    str = "range(1, 4):"
    @test replace_range(str) == "2:5"
    str = "range(1, 5):"
    @test replace_range(str) == "2:6"
    str = "range(1, n_grid):"
    @test replace_range(str) == "2:n_grid"
    str = "range(xs):"
    @test replace_range(str) == "1:xs"
    str = "    for i in range(1, 20):"
    @test replace_range(str) == "    for i in 2:21"
    str = "    for i in range(1, n_grid - 1):"
    @test replace_range(str) == "    for i in 2:n_grid"

    str = "len(xs)"
    @test replace_len(str) == "length(xs)"

    str = "linspace(half_delta - x_out, x_out - half_delta, n_grid)"
    @test replace_linspace(str) == "LinRange(half_delta - x_out, x_out - half_delta, n_grid)"

    str = "def construct(self, tol, it_max):"
    @test replace_def(str) == "function construct(self, tol, it_max)"

    str = "if x > 3:"
    @test replace_if(str) == "if x > 3"

    str = "while x < 4 && i == 0:"
    @test replace_while(str) == "while x < 4 && i == 0"

    str = "else:"
    @test replace_else(str) == "else"
    
    str = "zeros((n_grid, n_grid, n_grid))"
    @test replace_zeros(str) == "zeros(n_grid, n_grid, n_grid)"
    str = "zeros((n_grid, n_grid))"
    @test replace_zeros(str) == "zeros(n_grid, n_grid)"

    # TODO empty(()) ???

    str = "ones((n_grid, n_grid, n_grid))"
    @test replace_ones(str) == "ones(n_grid, n_grid, n_grid)"
    str = "ones((n_grid, n_grid))"
    @test replace_ones(str) == "ones(n_grid, n_grid)"

    str = "1 "
    @test handle_indexing_arg(str) == "2"
    str = "-1"
    @test handle_indexing_arg(str) == "end"
    str = " -2"
    @test handle_indexing_arg(str) == "end-1"
    strs = [" -2", "-1"]
    @test handle_indexing_arg.(strs) == ["end-1", "end"]

    strs = "alpha beta gamma sigma"
    @test replace_greek(strs) == "α β γ σ"

    str = "np.eye(5)"
    @test replace_eye(str) == "np.I(5)"

    str = "[0]"
    #@test replace_indexing(str) == "[1]"
    str = "[-1, n_grid, rho]"
    #@test replace_indexing(str) == "[end,n_grid,rho]"

    # Test multiple repeated indexes :D
    str = "ax[0] = ax[3,4] + x[3,4] + y[-1,y,5]"
    @test replace_indexing(str) == "ax[begin] = ax[4,5] + x[4,5] + y[end,y,6]"
end

