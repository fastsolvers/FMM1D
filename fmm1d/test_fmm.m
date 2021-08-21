addpath(genpath('./'))

fun = 1;                               % see Ker below
local_shift = false;            % using local shifting or not                      
scaling = 1;                        % using diagonal scaling or not
T = 0;  

n = 1e4;
n = 2*n;
s = 100*abs(randn(n,1));
s = sort(s);
y = s(1:2:end);
x = s(2:2:end);

nx = length(x);
ny = length(y);
org = [1:nx]';
gap = x-y;

k = floor(nx/2);
J = randsample(ny-1,k);
org(J) = org(J)+1;
gap(J) = x(J) - y(J+1);
q = ones(size(y,1),1);


switch fun
case 1
    if T == 0
        Ker = 1 ./ (x  - y.');
    else
        Ker = 1 ./ (y  - x.');
    end
case 2
    if T == 0
        Ker = 1 ./ (x  - y.').^2;
    else
        Ker = 1 ./ (y  - x.').^2;
    end
case 3
    if T == 0
        Ker = log(abs((x  - y.')));
    else
        Ker = log(abs((y  - x.')));
    end
end  
Ker(isinf(Ker)) = 0;
zl0 = tril(Ker) * q;
zu0 = triu(Ker,1) * q;


for r = 20:10:60
    fprintf('r = %d\n', r )    
    
    if T == 0 
        
        if local_shift 
            [zl, zu] = trifmm1d_local_shift( r, x, y, q, gap, org, fun, scaling);
            
        else
            [zl, zu] = trifmm1d( r, x, y, q, fun, scaling ); 
        end

    elseif T == 1
        
        if local_shift 
            [zl, zu] = trifmm1d_local_shift_2( r, y, x, q, gap, org, fun, scaling);
            
        else
            [zl, zu] = trifmm1d( r, y, x, q, fun, scaling ); 
        end

    end
    
    error_lower = norm( zl0  - zl ) / norm( zl0 );
    error_upper = norm( zu0 - zu ) / norm( zu0 );
    error = norm( zl0 + zu0 - zl - zu ) / norm( zl0 + zu0 );
    
    fprintf('fmm error: |phi-z|/|phi| = %e\n', max( error, max (error_lower, error_upper) ) );
    fprintf('\n\n')
end
