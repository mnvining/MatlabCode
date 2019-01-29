function [h_sol_1,h_sol_2]=hsol(z,Al)
    % z is fine of pts
    %Ar1=(exp(1/sqrt(Al)*(-z-1)));
    %Ar2=exp(1/sqrt(Al)*(z-1));
    %val1=norm(Ar1);
    %val2=norm(Ar2);
    val=z(end);
    syms x
    h_sol_1=exp(1/sqrt(Al)*(-x-1));%/val1
    h_sol_2=exp(1/sqrt(Al)*(x-val));%/val2
end