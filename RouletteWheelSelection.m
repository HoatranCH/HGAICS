function i=RouletteWheelSelection(P)

    r=rand;      % Tao ngau nghien R 
    
    c=cumsum(P); % Cong tong theo cot cua P
    
    i=find(r<=c,1,'first'); % Tim 1 so sao cho r<=c, 

end