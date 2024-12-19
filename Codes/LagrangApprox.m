%
% Written by Giovanni Braglia and Davide Tebaldi, 2023
% University of Modena and Reggio Emilia
% website: https://www.automatica.unimore.it/
%
%---------------------------------------------------------------------
function P = LagrangeApprox( yd, PSI, Ayg, byg, N )
%---------------------------------------------------------------------
% Computes analitically the weighting coefficients for the
% approximation of 'yd' with RBF basis functions

    L  = length( yd );
    W  = eye( L ) .* 1e0;
    Bi = Ayg;
    Ri = byg;
    Ba = PSI;
    Ra = yd;

    Ba2 = Ba' * W * Ba + eye(N).*realmin;
    BRa = Ba' * W * Ra;

    [ k, ~ ] = size( Bi );
    Q = ( Bi / Ba2 ) * Bi'  + eye(k).*realmin ;

    l = Q \ ( ( Bi / Ba2 ) * BRa - Ri );
    P = ( Ba2 \ BRa ) - ( Ba2 \ Bi' ) * l ;

end
