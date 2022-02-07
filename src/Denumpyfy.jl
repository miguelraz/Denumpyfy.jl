module Denumpyfy

# Write your package code here.


replace_exponents(str) = replace(str, "**" => "^")

# range(0, self.n_grid):
replace_range(str) = replace(str, r"range\((\w+), (\w+)\):$" => s"\1:\2")

#replace(str, r"range\(\w")

# range(\d, \w - \d):
# range(1, n_grid - 1):
# range(\d, \w):
# range(1, n_grid):
#
# linspace(half_delta - x_out, x_out - half_delta, n_grid)
replace_linspace(str) = replace(str, "linspace(" => "LinRange(")

# def construct(self, tol, it_max):
replace_def(str) = replace(str, r"def (.*):$" => s"function \1")
    
# if condition:
replace_if(str) = replace(str, r"if (.*):$" => s"if \1")
# while cond:
replace_while(str) = replace(str, r"while (.*):$" => s"while \1")
# else:
replace_else(str) = replace(str, r"else:$" => "else")
#
# zeros((n_grid, n_grid, n_grid))
replace_zeros(str) = replace(str, r"zeros\(\((.*)\)\)" => s"zeros(\1)")
# zeros((n_grid, n_grid))
# empty(()) == what in Julia???
#replace_empty(str) = replace(str, r"empty\(\((.*)\)\)" => s"empty(\1)")
# ones(())
replace_ones(str) = replace(str, r"ones\(\((.*)\)\)" => s"ones(\1)")
# [\d] increase?
#
# greek letters?
# \_alpha => \_α    
# alpha beta rho gamma delta

replacements! = [
    replace_exponents,
    replace_linspace,
]
export replace_exponents
export replace_range
export replace_linspace
export replace_def
export replace_if
export replace_else
export replace_while
export replace_zeros
export replace_ones


end
