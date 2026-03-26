# Harnessing Symbolic Regression and Nonlinear Programming for Computationally Expensive Multi-Objective Optimization

Symbolic regression provides the explicit analytical expression which allows us to use the efficient constraint handling approach available for explicit equations.

The proposed approach start by initial solutions using Latin Hypercube Sampling (LHS), followed by the true evaluations of objectives and constraints and archive them. Then, the symbolic regression using PySR toolbox has been built for all the objectives and constraints.
Then we perform a warm-up stage, nothing but getting the extremal solutions of the pareto front by solving for each objectives on a inexpensive symbolic regression and the solutions are truly evaluated and updated in the archive that is used to update the symbolic regression. 
Next, the scalarized subproblem on different reference points using Tchebycheff scalarization solved on the cheap symbolic regression to get the candidate solutions. Then, a single candidate solution is selected using Distance based Subset Selection (DSS) for true evaluation.
The evaluated solution is used to update the surrogate and the process continues until the evaluation budget exhausted.
