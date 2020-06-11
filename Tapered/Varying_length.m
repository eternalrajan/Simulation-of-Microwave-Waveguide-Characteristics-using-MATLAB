%% Finding zero X'np of differential of bessel
syms x; 
n=0;                                        % nth TE mode
p=1;                                        % pth zero 
i=0;                                        % Variable
c=1;                                        % Counter
while i<50                                  % Set loop limit high for higher modes
    s=vpasolve(diff(besselj(n,x))==0,x,i);  % Guessing solution
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
Xdnp=sols(p);                               % pth zero of nth order of TEnp
% Effect of variation in length of waveguide on cutoff frequency
rW1=.00373;                                 % Initial radius at z=0
rW2=.00737;                                 % Initial radius at z=l
l=0;                                        % Waveguide length
while l<=.1
    z=0:.001:l;                             % Axial distance
    tanang=(rW2-rW1)/l;                     % Tangent of tapering angle
    rz=rW1+rW1*tanang*z./(rW2-l*tanang);    % Radius at z axial distance 
    f=3e8*Xdnp./(2*pi*rz);                  % Cutoff frequency
    hold on;
    plot(z,f);
    xlabel('Axial distance z');
    ylabel('Cutoff frequency f');
    l=l+.01;                                % Increasing waveguide length
end
legend('L=0','L=0.01','L=0.02','L=0.03','L=0.04','L=0.05','L=0.06','L=0.07','L=0.08','L=0.09','L=0.1');
