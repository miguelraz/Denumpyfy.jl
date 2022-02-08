module Denumpyfy

# Write your package code here.


replace_exponents(str) = replace(str, "**" => "^")

# range(0, self.n_grid):
function replace_range(str)
    isempty(str) && return str
    !occursin(r"range\(", str) && return str
    caps = match(r"range\((.*)\):$", str).captures[1]
    caps = filter(!isspace, caps)

    # Handle single argument case by checking for existence of a comma
    # TODO Assume there is a `len` call?
    !occursin(r",", caps) && return string("1:", caps)

    # Try and be smart about 1-based indexing
    first, second = split(caps, ',')

    # Handle either one being a positive digit and make it 1-based index
    if all(isdigit, first)
        first = parse(Int, first) + 1
    end
    if all(isdigit, second)
        second = parse(Int, second) + 1
    end
    if second == "-1"
        second = "end"
    end
    if second == "-2"
        second = "end-1"
    end
    if second == "-3"
        second = "end-2"
    end
    # If you need more may the lord have mercy on your code 
    return string(first,':', second)
end

replace_eye(str) = replace(str, "eye(" => "I(")

# INPUT:
# A possible indexing like
# 1, 2 => 1,2
# -1, ngrid - 1 => end,ngrid-1
# 0, 1 => 1,2
function handle_indexing_arg(str)
        str = filter(!isspace, str)
        # Special case "0"
        str == "0" && return "begin"
        # If is all digits, bump
        all(isdigit, str) && return string(parse(Int, str) + 1)
        # if uses "-1" or some other negative integer
        str == "-1" && return "end"
        str == "-2" && return "end-1"
        str == "-3" && return "end-2"
        str == "-4" && return "end-3"
        str == "-5" && return "end-4"
        str == "-6" && return "end-5"
        str == "-7" && return "end-6"
        str
end

replace_len(str) = replace(str, "len(" => "length(")

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

# zeros((n_grid, n_grid, n_grid))
replace_zeros(str) = replace(str, r"zeros\(\((.*)\)\)" => s"zeros(\1)")
# empty(()) == what in Julia???
# ones((n_grid, n_grid))
replace_ones(str) = replace(str, r"ones\(\((.*)\)\)" => s"ones(\1)")

# TODO handle slices

# [\d] increase?
function replace_indexing(str)
    isempty(str) && return str
    !occursin(r"\[|\]", str) && return str

    for m in eachmatch(r"\[(.+?)\]", str, overlap = false)
        cap = (m.match)[2:end-1] # get the stuff inside the [...]
        splits = split(cap, ',')
        idx = handle_indexing_arg.(splits)
        idx = '[' * join(idx, ',') * ']'
        str = replace(str, m.match => idx)
    end
    return str
end

# More/other letters possible but maybe are problematic with collisions :/
function replace_greek(str)
    for (greek, letter) in zip("αβργΔϵθκστω", ["alpha", "beta", "rho", "gamma", "delta", "epsilon", "theta", "kappa", "sigma", "tau", "omega"])
        str = replace(str, letter => greek)
    end
    str
end

replacements! = [
    replace_exponents,
    replace_range,
    replace_len,
    replace_linspace,
    replace_def,
    replace_if,
    replace_else,
    replace_while,
    replace_zeros,
    replace_ones,
    replace_greek,
    replace_eye,
    replace_indexing
]

function denumpyfy(file)
    jlfilename = replace(file, r"^(\w+).py$" => s"\1.jl")
    contents = readlines(file, keep = false) .|> String
    for f in replacements!
        contents .= f.(contents)
    end
    open(jlfilename, "w") do x
        for s in contents
            println(x, s)
        end
    end
end


export replace_exponents
export replace_range
export replace_len
export replace_linspace
export replace_def
export replace_if
export replace_else
export replace_while
export replace_zeros
export replace_ones
export replace_greek
export handle_indexing_arg
export replace_eye
export replace_indexing
export denumpyfy


end
