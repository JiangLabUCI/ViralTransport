function pinf = fp(dose,P,mua)
%FP Vectorized fractional Poisson dose-response model.
%   C = FP(dose, P, mua) computes probability of infection assuming 
%   a fractional Poisson dose-response model. dose can be a vector. P and
%   mua must be scalars.
    pinf = P * (1-exp(-dose/mua));
end