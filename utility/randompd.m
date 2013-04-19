% A MATLAB script to generate random data of strictly convex bqp
function [] = randompd(n,density,rc,type,outname)
    % Proprocess the input arguments
    if nargin < 4
        type = 1;
	end
    if nargin < 5
        outname = 'tmp';
    end

    % Generate Hessian
    H = sprandsym(n,density,rc,type);
    % Generate linear coefficients
    c = randn(n,1);
    % Generate upper bounds
    u = randn(n,1);
    % Save data
    save(outname,'H','c','u');
end
