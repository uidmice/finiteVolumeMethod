function [rhoR, etaR] = test (x, rho0, eta0, W11,W12,W21,W22, e, T)

clear model_function
N = length(x)-1;
rhoR = zeros(T+1, N);
etaR = zeros(T+1, N);

rhoN = rho0;
etaN = eta0;

rhoR(1,:) = rho0;
etaR(1,:) = eta0;

tic;
for i=1:T
    if mod(i,1000)==1
        clear model_function;
    end
    [t1, t2] = model_function (x, rhoN, etaN, W11,W12,W21,W22, e);
    rhoR(1+i,:) = t1;
    etaR(1+i,:) = t2;
    rhoN = t1;
    etaN = t2;
end

t = toc;
disp([ num2str(t) ' sec for ' num2str(T) ' runs']);
end

