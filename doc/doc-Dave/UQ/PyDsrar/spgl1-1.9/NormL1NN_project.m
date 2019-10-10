function x = NormL1NN_project(x,weights,tau)
% Projection onto the non-negative part of the L1 ball

%idx = (x > 0);
%x(idx) = NormL1_project(x(idx),weights(idx),tau);
x(x < 0) = 0;
x = NormL1_project(x,weights,tau);
