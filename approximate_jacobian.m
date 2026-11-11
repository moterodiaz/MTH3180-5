% Copied unchanged (apart from this attribution) from ../MTH3199-Assignment-3/approximate_jacobian.m
% Source commit: 2c614d9f8d148b162ed373f79ef71a0706bce013.
% Original work: this user's Assignment 3 numerical-methods implementation.
%Implementation of finite difference approximation
%for Jacobian of multidimensional function
%INPUTS:
%fun: the mathematical function we want to differentiate
%X: the input value of fun that we want to compute the derivative at
%OUTPUTS:
%J: approximation of Jacobian of fun at x
% num_evals: the number of times a rate function has been called
%Implementation of finite difference approximation
%for Jacobian of multidimensional function
function [J,num_evals] = approximate_jacobian(fun,X)
    dx_scalar = 1e-6;

    f0 = fun(X);
    num_evals = 1; % starting a new count w the first look

    num_in = length(X);
    num_out = length(f0);

    dX = zeros(size(X));
    J = zeros(num_out,num_in);

    for n = 1:num_in
        dX(n) = dx_scalar;

        f_plus = fun(X + dX);
        f_minus = fun(X - dX);
        num_evals = num_evals + 2;

        dfdx = (f_plus - f_minus)/(2*dx_scalar);
        J(:,n) = dfdx;

        dX(n) = 0;
    end
end
% I removed the counter so counts do not carry between indiv calls
