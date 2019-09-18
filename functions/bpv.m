function pinf = bpv(dose, alpha, beta)
    %BPV Vectorized beta-Poisson dose-response model.
    %   C = BPV(dose, alpha, beta) computes probability of infection assuming
    %   a beta-Poisson dose-response model. dose can be a vector. alpha and
    %   beta must be scalars.
    pinf = 1 - (1 + dose / beta).^(-alpha);
end
