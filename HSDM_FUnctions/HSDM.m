function [Err,varargout]=...
    HSDM(N, x, t, Dose, rp, qo, isothermFunc, equiParams, c, nr, co, otherArgs)

arguments
    
    N (1,1) double {mustBeInteger, mustBeFinite, mustBePositive}
    x (1,:) double {mustBePositive, mustBeFinite}
    t (:,1) double {mustBeNonnegative, mustBeFinite}
    Dose (1,1) double {mustBePositive, mustBeFinite}
    rp (1,1) double {mustBePositive, mustBeFinite}
    qo (1,:)  double {mustBePositive, mustBeFinite}
    isothermFunc function_handle
    equiParams (:,:) double {mustBePositive, mustBeFinite}
    c (:,:) double {mustBeNonnegative, mustBeFinite} = zeros(size(t,1),N)
    nr (1,1) double {mustBePositive, mustBeFinite, mustBeInteger, mustBeGreaterThan(nr,51)}= 101
    co (1,:)  double {mustBeNonnegative, mustBeFinite} = c(1,:)
    otherArgs.scale (1,1)  double {mustBePositive, mustBeFinite} = 1e-15
    otherArgs.Display (1,1) string {mustBeMember(otherArgs.Display,{'off', 'iter', 'final'})} = 'off'
    
end

persistent Display

Ds = x( 1:N)*otherArgs.scale

DB = Dose*qo./co;
g = Ds./Ds(1);

Tb = Ds(1)*t*60/rp^2;
Yini = [zeros(1,N*nr), ones(1,N)];
r = linspace(0,1,nr);


try
    
    [To,Yo]=ode15s(@model,Tb,Yini);
    
    X=Yo(:,end-N+1:end);
    Yo=Yo(:,1:end-N);
    
    if N~=1
        
        Yo=permute(reshape(Yo,size(To,1),nr,N),[2,1,3]);
        
    else
        
        Yo=Yo';
        
    end
    
    qavg = reshape((3*trapz(r,((r'.^2).*Yo))),length(t),N).*qo;
    
    cp = (X.*co);
    Err=sum(sum((c-cp).^2));%+sum(sum((Xl-cp).^2));
    
catch ME
    
    Err = inf;
    cp = zeros(size(c));
    Yo = zeros([length(t),nr,N]);
    qavg = zeros([length(t),N]);
    rethrow(ME)
end

if ~strcmp(otherArgs.Display, 'off')
    
    if isempty(Display)
        Display = table(1, Ds, Err);
        Display.Properties.VariableNames = ["Iteration", strjoin("Ds"+string(1:N),' & '), "Error"];
        disp(Display)
    end
    clc
    Display(end+1,:) = {Display.Iteration(end)+1, Ds, Err};
    disp(Display)
    
end
[varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}]=deal(t,c,cp,Yo,qavg);

    function Yt=model(~,Y)
        
        X = Y(end-N+1:end);
        Y = Y(1:(end-N));
        Y = reshape(Y,nr,N)';
        Y(:,end) = (isothermFunc(equiParams,(X'.*co))./qo)';
        Yr = d1AG(0,1,nr,Y);
        Yr(:,1) = 0;
        Yrr = d2AG(0,1,nr,Y,Yr,2,2);
        
        op = 2:nr;
        Yt = g'.*[3*Yrr(:,1),(Yrr(:,op)+(2./r(:,op)).*Yr(:,op))];
        Xt = -(3*g.*DB.*((Yr(:,end)')));
        
        Yt = Yt';
        Yt = [Yt(:);Xt'];
        
    end
end
