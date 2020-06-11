n=0;                                            % nth order of TEnp mode
while n<10
    % Finding zero X'np of differential of bessel
    syms x;
    p=1;                                        % pth zero
    i=0;                                        % Variable
    c=1;                                        % Counter
    while i<35                                  % Set loop limit high for higher modes
        s=vpasolve(diff(besselj(n,x))==0,x,i);  % Guessing solution
        if s>=10^-20                            
            values(c)=s;        
            if ~isempty(values(c))
                c=c+1;
            end
        end
        i=i+1;
    end
    values=sort(values);                        % Sorting zeros found
    c=1;
    d=-1;
    i=1;
    while i<=length(values)                     % Elimination of repeated values
        if abs(d-values(i))>=10^-20
            d=values(i);
            sols(c)=d;
            c=c+1;
        end
        i=i+1;
    end
    Xdnp=sols(p);                               % pth zero of nth order of TEnp
    % Dispersion Relation
    syms w b;
    c=3*10^8;                                   % Speed of light
    rW2=0.00737;                                % Final radius
    kc=Xdnp/rW2;                                % Cutoff wavenumber
    vt=1.32*10^8;                               % Transversal electron beam velocity
    vz=1.32*10^8;                               % Axial electron beam velocity
    s=1;                                        % Beam mode harmonics                
    wc=c*kc;                                    % Angular cutoff frequency
    g=(1-(vt^2+vz^2)/c^2)^-.5;                  % Relativistic mass factor
    e=1.6*10^-19;                               % Electronic charge
    uo=4*pi*10^-7;                              % Permeablity of free space
    No=10^4;                                    % Number of electrons per unit axial length
    m=0;                                        % Angular harmonic mode number
    rH=0.0001865;                               % Average hollow beam radius
    rL=.000373;                                 % Electron Larmor Radius
    meo=9.1*10^-31;                             % Mass of electron
    rW1=.00373;                                 % Inital radius
    Jds=@(x) (besselj(s+1,x)-besselj(s-1,x))/2; % J's(x)
    Jdm=@(x) (besselj(m+1,x)-besselj(m-1,x))/2; % J'm(x)
    num=uo*e^2*No*(vt/c)^2*(besselj(s-m,kc*rH))^2*Jds(kc*rL)^2;                  
    den=g*meo*pi*rW1^2*(1-(m/(kc*rW1))^2)*Jdm(kc*rW1);                                 
    f=@(B,w) ((w./c).^2-B.^2-kc^2).*(w-B.*vz-s*wc/g).^2+num*(w.^2-B.^2.*c^2)./den;
    hold on;
    fimplicit(f);                               % Plotting dispersion relation
    xlabel('Propagation constant B');           % X axis label
    ylabel('Angular frequency w');              % Y axis label
    xlim([-2000 2000]);                         % X limit
    ylim([0 8*10^11]);                          % Y limit
    n=n+1;                                      % Increasing n
end
legend('TE01','TE11','TE21','TE31','TE41','TE51','TE61','TE71','TE81','TE91');
