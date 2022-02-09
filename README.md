# Denumpyfy.jl

I hate translating numpy notebooks to Julia.

I then spent [about 2 afternoons packing 20 regexes in to a trenchcoat](https://xkcd.com/1205/) and calling it a Julia package so that others may join me in not fiddling with indexing of nasty, nasty "vectorized" code.

My code is quite naive - it tries to be a bit smart about line by line regex replacing, and doesn't really handle global context or `self` rewriting - and likely won't. It was written with **replacing tutorial numpy notebook tutorials in mind**. If you are looking for something that is much more robust/heavy duty, please check out [Py2Jl.jl -> py2many](https://github.com/adsharma/py2many). You can compare `Denumpyfy.jl` vs `py2many` with `puncture.jl` vs `tests/demopy2many.jl`, respectively. At time of writing, their AST parser is an unfinished copy/paste from their Rust transpiler, so please help out by filing PRs/issues over there.

Use entirely at your own risk, and see some of the `tests/runtests.jl` to see what the package can help you out with.

For example, this code in `test/puncture.py` 
```python
for i in range(0, self.n_grid-1):
    for j in range(0, self.n_grid-1):
        for k in range(0, self.n_grid-1):
            delta_x = x[i, -1] ** 2 + y[j, -2] ** 2 + z[k, -3] ** 2
```
will be turned into this with `denumpyfy("puncture.py")`
```julia
for i in 1:self.n_grid
    for j in 1:self.n_grid
        for k in 1:self.n_grid
            Δ_x = x[i, end] ^ 2 + y[j, end-1] ^ 2 + z[k, end-2] ^ 2
```
in a new file `puncture.jl`.

There's `replace` rules for array indexing (with 1-index bumping), some ascii -> Unicode replacements, `range(0,n-1):` to `1:n`, and the obvious arithmetic `**` to `^`.

TL;DR:

![You should be sponsoring me not reading this](friendshipendedwithnumpy.jpg "No more need for vectorized code")

-----

### Sponsorship

If you want to see more of this work, consider sponsoring me via [Github sponsors](https://github.com/sponsors/miguelraz/).

### Acknowledgments

This work was done via nerdsnipe and support of `Brenhin Keller`.
