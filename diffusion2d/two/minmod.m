function [m]=minmod(a,b,c)
V = [a,b,c];
if any(diff(sign(V(V~=0))))
    m = 0;
else
    if sum(V)>0
        m = min(V);
    else
        m = max(V);
    end
end
    
end
