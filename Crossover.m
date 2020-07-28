function [y1, y2]=Crossover(x1,x2,gamma,VarMin,VarMax,nVar)

    alpha=unifrnd(-gamma,1+gamma,size(x1)); % he so lai ghep
    for j = 1:nVar
    y1(1,j)=alpha(1,j)*x1(1,j)+(1-alpha(1,j))*x2(1,j);
    y2(1,j)=alpha(1,j)*x2(1,j)+(1-alpha(1,j))*x1(1,j);
    end
    for j = 1:nVar
    y1(1,j)=max(y1(1,j),VarMin(j,1));
    y1(1,j)=min(y1(1,j),VarMax(j,1));
    
    y2(1,j)=max(y2(1,j),VarMin(j,1));
    y2(1,j)=min(y2(1,j),VarMax(j,1));
    end
end