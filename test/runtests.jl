using Denumpyfy
using Test

@testset "Denumpyfy.jl" begin
    str = "x ** 2"
    @test replace_exponents(str)  == "x ^ 2"

    str = "range(0, self.n_grid):"
    @test replace_range(str) == "0:self.n_grid"

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
end
