# StandardBasisVectors

Re-orient a set of orthonormal vectors to create a proper right-handed coordinate systems. Based on the work in Damask, Jay. “A Consistently Oriented Basis for Eigenanalysis.” International Journal of Data Science and Analytics 10, no. 4 (October 2020): 301–19. [ttps://doi.org/10.1007/s41060-020-00227-z](https://doi.org/10.1007/s41060-020-00227-z).

## Usage

```julia
# left-handed system
Q = [-1.0 0.0;
      0.0 1.0]

 V = standardize_basis(Q)
```
