clear;
% Define the function
function y = func1(x)
    y = exp(x) - cos(x) - 2;
end

% Bisection Method
function root = bisection_method(a, b, tol)
    if func1(a) * func1(b) >= 0
        disp('Bisection method fails.');
        root = NaN;
        return;
    end
    c = a;
    iter = 1;
    error_approx = (b - a) / 2.0;
    while  error_approx > tol
        c = (a + b) / 2.0;
        if func1(c) == 0
            break;
        elseif func1(a) * func1(c) < 0
            b = c;
        else
            a = c;
        end
        fprintf('Bisection Method: Iteration %d: Approximate Root = %.10f, Error = %.10e\n', iter, c,error_approx);
        iter=iter + 1;
    end
    root = c;
end

% Inverse Quadratic Interpolation Method with detailed output
function root = inverse_quadratic_interpolation(func2, x0, x1, x2, tol, max_iter)
    % Evaluate the function at initial guesses
    f0 = func2(x0);
    f1 = func2(x1);
    f2 = func2(x2);
    
    for iter = 1:max_iter
        % Check if any of the function values are sufficiently close to zero
        if abs(f0) < tol
            root = x0;
            fprintf('IQI Method: Iteration %d: Approximate Root = %.10f, Error = %.10e\n', iter, x0, abs(f0));
            return;
        elseif abs(f1) < tol
            root = x1;
            fprintf('IQI Method: Iteration %d: Approximate Root = %.10f, Error = %.10e\n', iter, x1, abs(f1));
            return;
        elseif abs(f2) < tol
            root = x2;
            fprintf('IQI Method: Iteration %d: Approximate Root = %.10f, Error = %.10e\n', iter, x2, abs(f2));
            return;
        end
        
        % Inverse quadratic interpolation formula
        denominator = (f0 - f1) * (f0 - f2) * (f1 - f2);
        if denominator == 0
            error('Denominator became zero, method failed.');
        end
        
        % Compute the next approximation for the root
        x_new = x0 * f1 * f2 / ((f0 - f1) * (f0 - f2)) ...
              + x1 * f0 * f2 / ((f1 - f0) * (f1 - f2)) ...
              + x2 * f0 * f1 / ((f2 - f0) * (f2 - f1));
        
        % Update values for next iteration
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
        x2 = x_new;
        f2 = func2(x_new);
        
        % Print status for each iteration
        fprintf('IQI Method: Iteration %d: Approximate Root = %.10f, Error = %.10e\n', iter, x_new, abs(f2));
        
        % Check convergence based on function value
        if abs(f2) < tol
            root = x_new;
            fprintf('IQI Method: Root found at Iteration %d: Approximate Root = %.10f, Error = %.10e\n', iter, x_new, abs(f2));
            return;
        end
    end
    
    % If maximum iterations are reached without convergence
    error('Maximum number of iterations reached without convergence.');
end

% Bisection method
a = 0;
b = 1;
tol = 1e-6; % Tolerance
root_bis = bisection_method(a, b, tol);
fprintf('Newton Bisection Method Final Root: %.10f\n', root_bis);

% Define the polynomial function for IQI
func2 = @(x) exp(x) - cos(x) - 2; % Define the polynomial function
x0 = 0;     % First initial guess
x1 = 0.93;   % Second initial guess
x2 = 1;     % Third initial guess
max_iter = 100; % Maximum number of iterations

% Call inverse quadratic interpolation function
root_iqi = inverse_quadratic_interpolation(func2, x0, x1, x2, tol, max_iter);
fprintf('Inverse Quadratic Interpolation Final Root: %.10f\n', root_iqi);

% Call MATLAB's built-in fzero function
% fzero requires an initial guess or interval. You can provide one of your initial guesses.
root_fzero = fzero(func2, [x0 x2]); % Providing the interval between x0 and x2
fprintf('MATLAB fzero Root: %.10f\n', root_fzero);

% Compare the results
fprintf('Difference between IQI and fzero: %.10e\n', abs(root_iqi - root_fzero));
