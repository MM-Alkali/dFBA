% Configure solver - use GLPK if Gurobi not available
try
    changeCobraSolver('gurobi', 'LP');
    fprintf('Using Gurobi solver\n');
catch
    changeCobraSolver('glpk', 'LP');
    fprintf('Using GLPK solver (Gurobi not available)\n');
end