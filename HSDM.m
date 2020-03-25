function [Err,varargout]=...
    HSDM(x,N,t,c,nr,dose,rp,qo,co,scale,isotherm,eqp,sim)
global Xti Tbx

Xti=[];Tbx=[];
Ds=x( 1:N).*scale(1);
[DB,g]=deal((dose*qo./co),Ds./Ds(1));
[Tb,Yini,r]=deal((Ds(1)*t*60/rp^2),([zeros(1,N*nr),(ones(1,N))]),linspace(0,1,nr));


try
[To,Yo]=ode15s(@model,Tb,Yini);
X=Yo(:,end-N+1:end);
Yo=Yo(:,1:end-N);
if N~=1
Yo=permute(reshape(Yo,size(To,1),nr,N),[2,1,3]);    
else
    Yo=Yo';
end
qavg=reshape((3*trapz(r,((r'.^2).*Yo))),length(t),N).*qo;
cp=(X.*co);
if sim==2

   Err=0;
    c=cp;
else
    Err=sum(sum((c-cp).^2));%+sum(sum((Xl-cp).^2));
end
catch
    Err=inf;
    cp=zeros(size(c));
    Yo=zeros([length(t),nr,N]);
    qavg=zeros([length(t),N]);
end



disp([x,Err])

[varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}]=deal(t,c,cp,Yo,qavg);

    function Yt=model(ts,Y)
        X=Y(end-N+1:end);
        Y=Y(1:(end-N));
        Y=reshape(Y,nr,N)';
        Y(:,end)=(Equilibrium(isotherm,'ce2qe',eqp,(X'.*co))./qo)';
        Yr=d1AG(0,1,nr,Y);
        Yr(:,1)=0;
        Yrr=d2AG(0,1,nr,Y,Yr,2,2);
        
        op=2:nr;
        Yt=g'.*[3*Yrr(:,1),(Yrr(:,op)+(2./r(:,op)).*Yr(:,op))];
        Xt=-(3*g.*DB.*((Yr(:,end)')));
        
        Xti=[Xti;Xt];Tbx=[Tbx;ts];
        
        Yt=Yt';
        Yt=[Yt(:);Xt'];
    end
end
