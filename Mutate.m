function y=Mutate(x,mu,VarMin,VarMax,nVar)

    nmu=ceil(mu*nVar);  % Lam tron
    
    j=randsample(nVar,nmu); % Lay so ngau nhien tu nVar den nmu
    
    sigma=0.1*(VarMax(j,1)-VarMin(j,1)); % He so dot bien
    
    y=x;
    y(j)=x(j)+sigma*randn(size(j));
    
    y=max(y,VarMin(j,1));
    y=min(y,VarMax(j,1));

end