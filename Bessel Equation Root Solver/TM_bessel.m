%% Finding zero Xnp of bessel
syms x; 
n=0;                                        % nth TM mode
p=1;                                        % pth zero 
i=0;                                        % Variable
c=1;                                        % Counter
while i<35                                  % Set loop limit high for higher modes
    s=vpasolve(besselj(n,x)==0,x,i);        % Guessing solution
    if s>=10^-20                            
        values(c)=s;        
        if ~isempty(values(c))
            c=c+1;
       	end
    end
    i=i+1;
end
values=sort(values);                        % Sorting values found
c=1;
d=-1;
i=1;
while i<=length(values)                     % Elimination of repeated zeros
    if abs(d-values(i))>=10^-20
	d=values(i);
        sols(c)=d;
        c=c+1;
    end
    i=i+1;
end
Xnp=sols(p)                                 % pth zero of nth order of TMnp