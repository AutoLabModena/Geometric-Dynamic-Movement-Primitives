%
% Written by Giovanni Braglia and Davide Tebaldi, 2023
% University of Modena and Reggio Emilia
% website: https://www.automatica.unimore.it/
%
%---------------------------------------------------------------------
function ddphi = second_der( rbf, drbf, ddrbf, sum_rbf, dsum_rbf, ddsum_rbf )
%---------------------------------------------------------------------
% Calculates analitically the second derivative of normalized RBF


    a = @(s) ddrbf(s) .* sum_rbf(s) - rbf(s) .* ddsum_rbf(s);
    b = @(s) drbf(s) .* sum_rbf(s) - rbf(s) .* dsum_rbf(s);


    num = @(s) a(s) .* ( sum_rbf(s).^2 ) -  b(s) .* ( 2 * sum_rbf(s) ) .* dsum_rbf(s);
    den = @(s) ( sum_rbf(s).^4 );

    ddphi = @(s) num(s) ./ den(s);

end
