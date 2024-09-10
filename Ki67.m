
    % Define symbolic variables
    syms a t ta AI PI
    % Define the symbolic equation
    equation = a * (1 - AI - PI) - ((1/t) * (PI + PI^2) - (1/ta) * AI * PI) == 0;
    % Constants
    a_val = 0.00273;
    t_val = 24/24;
    ta_val = 8/24;
    AI_val = 0.7/100;
    % Substitute specific values into the equation
    equation_substituted = subs(equation, [a, t, ta, AI], [a_val, t_val, ta_val, AI_val]);
    
    % Solve the equation for PI
    solutions = solve(equation_substituted, PI);
    
    % Filter out the negative solution
    positive_solution = max(double(solutions)); % Convert to double and take the max, assuming at least one positive solution
    
    % Check if the solution is positive
    if positive_solution < 0
        error('No positive solution found.');
    else
        % Assign the positive solution to the output variable
        ki67 = positive_solution*100;
    end
