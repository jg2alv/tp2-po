# Pivot Point: a Julia simplex solver for optimal, infeasible, and unbounded programs

Pivot Point implements the simplex method directly in Julia, without matrix or linear-algebra libraries. It works with equality constraints and nonnegative variables, aiming to distinguish a finite optimum from infeasibility and unboundedness rather than returning only an objective value.

## 1. Academic context and scope

- **Institution:** Universidade Federal de Minas Gerais (UFMG)
- **Degree:** Ciência da Computação
- **Course:** Pesquisa Operacional
- **Assignment:** Trabalho Prático 2
- **Semester:** 2025-02

The program addresses linear programs of the form `maximize cᵀx`, subject to `Ax = b` and `x ≥ 0`. Its scope includes negative coefficients and right-hand sides. The single source file contains input parsing, auxiliary-problem construction, tableau pivoting, basis tracking, termination logic, and output formatting.

## 2. How to run

The source uses Julia and its standard `Printf` library. From the repository root, set `INPUT_FILE` to an existing instance path and `OUTPUT_FILE` to the destination path, then invoke the entry point:

```sh
julia simplex.jl "$INPUT_FILE" "$OUTPUT_FILE"
```

The parser expects whitespace-separated tokens: the number of constraints `n`, the number of variables `m`, then `m` objective coefficients. Each of the next `n` rows contains `m` coefficients followed by its right-hand side. Both paths are required positional arguments. The output file is created or overwritten.

The source writes `inviavel` for infeasibility, `ilimitada` for unboundedness, or `otima` followed by the objective with three decimal places. Those tokens remain Portuguese because they are part of the program's output contract. The repository does not pin a Julia version.

## 3. Implementation decisions

The solver first builds an auxiliary tableau with artificial variables, then reconstructs a tableau for the original objective. Explicit basic and nonbasic index vectors track which variables participate in the current basis. Rows with negative right-hand sides are multiplied by minus one during input preparation.

Entering-variable selection favors the smallest eligible variable index. The leaving-variable test compares positive pivot coefficients and breaks tied ratios using basis indices. Pivoting normalizes the selected row and eliminates that column from other rows through explicit loops, making the numerical mechanics visible.

## 4. Learning and limitations

Implementing the method exposed the connection between an abstract basis and its numerical tableau representation. Initialization, basis changes, and the distinction between infeasibility and an unbounded objective became concrete programming responsibilities.

The implementation uses `Float64`, an `EPS` tolerance of `1e-9` in pivot decisions, and a separate `1e-6` auxiliary-objective check. These choices mean numerical comparisons are approximate. The repository contains no instance fixtures or automated tests, so correctness across degeneracy, redundant constraints, and difficult numerical scales is not established by the available materials. The output reports status and objective rather than a full primal solution or certificate.
