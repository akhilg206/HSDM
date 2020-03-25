clear all
clc

%% Data acquisition
% Open data files saved as .mat files in Data folder
% Data only of what was used in EHO paper
% Find Home directory
cur_dir=cd;                     % Current folder location
idcd=strfind(cur_dir,'\');      % Find the last '\' to navigate to home
home_dir=cur_dir(1:idcd(end)-1);% Set Home directory


addpath(home_dir+"\Data",home_dir+"\DSS_files",home_dir+"\EHO_fun") % Add path to required folders


% Acquire Run details as given in Table 3 of EHO paper
load("Run_det")

% Acquire Run data as per the details
load("Runs")

use_new_data='No';

if strcmp(use_new_data,'Yes')
    Run_Det=[]; % Same format as given in Table 3 of EHO paper
    Runs={[]};  % column vector of [time (min), concetration (mg/L)] in cell
    Co  =[]; %mg/L
    Dose=[]; %g/L
    RPM =[]; % Rate of agitation
    rhop=[]; % Density of Adsorbent, kg/m3
    rp  =[]; % Radius of adsorbent,  m
    kF  =[]; % Freundlich isotherm parameter
    nF  =[]; % Freundlich isotherm parameter
else
    
    
    %% Process conditions
    % Uncomment as needed
    % FOr iterating through all conditions
    Dose=[2,5,8];
    Co  =[250,500,750];
    RPM =[0,150,300];
    
    %     % For iteraing only at given conditions
    %     Dose=[2];
    %     Co  =[750];
    %     RPM =[150];
    %
    %% Adsorbent details
    rhop    = 1203.13;         % Density of Adsorbent, kg/m3
    rp      = 0.5*.71e-3;      % Radius of adsorbent,  m
    kF      = 30.86;            % Freundlich isotherm
    nF      = 0.28;
end

%% Finite difference mesh setup
nr= 21; % number of nodes
r= linspace(0,1,nr); % Discretized radial vector

%% Select optimization algorihtm
%% Optization algos as given in EHO papaer
%% NMA/SA/GA/PSO/EHO/ cycle through all
OA=["EHO"];

%% Execute HSDM Parameter estimation with set values
% if OA is all then run all optimization algos
if strcmp(OA,'all')
    OA=["NMA","SA","GA","PSO","EHO"];
end

% Number of function calls for compartitive study
global nfcall
for i=1:numel(OA)
    %     for mv=Dose
    %         for co=Co
    %             for rpm=RPM
    %                 try %catch error if any
%     Run=Run_det(Run_det(:,2)==mv & Run_det(:,3)==co & Run_det(:,4)==rpm);
    %                                         Run=4;
    for j=1:17
        mv=Run_det(j,2);
        co=Run_det(j,3);
        rpm=Run_det(j,4);
        Run=j
        for iteration=1:5
            iteration
            ct_exp_data=data{j};
            t=ct_exp_data(:,1);
            c=ct_exp_data(:,2);
            LB=[1
                1e-4];
            UB=[3e3
                2e3];
            nVar=2;
            nfcall=0;
            IG=unifrnd(LB,UB);  % Random Intitial guess
            Tol=1e-2;
            tic % For computation time
            OA(i)
            switch OA(i)
                
                case "NMA"
                    options=optimset('TolFun',Tol,'Display','iter');
                    x    =fminsearch( @(x) HSDM(x,t,c,kF,nF,rhop,rp,mv,co,nr,r),IG',options);
                case "SA"
                    options=optimoptions('simulannealbnd','FunctionTolerance',Tol,'Display','iter');
                    x=simulannealbnd( @(x) HSDM(x,t,c,kF,nF,rhop,rp,mv,co,nr,r),IG',LB,UB,options);
                case "GA"
                    options=optimoptions('ga','FunctionTolerance',Tol,'Display','iter');
                    x=            ga(@(x)  HSDM(x,t,c,kF,nF,rhop,rp,mv,co,nr,r),nVar,[],[],[],[],LB,UB,[],options);
                case "PSO"
                    options=optimoptions('particleswarm','FunctionTolerance',Tol,'Display','iter');
                    x= particleswarm(@(x)  HSDM(x,t,c,kF,nF,rhop,rp,mv,co,nr,r),nVar,LB,UB,options);
                case "EHO"
                    options=struct('numMotherE',12,'numCalves',5,'Tolerance',Tol,'OmegaStar',10,'Display','iter');
                    x=ElephantHerd(@(x)  HSDM(x,t,c,kF,nF,rhop,rp,mv,co,nr,r),nVar,LB,UB,options);
                    
            end
            t_comp=toc; % Computational time for Parameter estimation
            tic; % For simulation time
            [Er,cp]=HSDM(x,t,c,kF,nF,rhop,rp,mv,co,nr,r);
            t_sim=toc; % Time required per simulation
            fprintf('\n\n\t Computational Time %d \t\t Simulation time %d\n\n\t\t Estimated parameters: %d %d \n\t\t SSE: %d\n',t_comp,t_sim,x,Er)
            
            % Write file to excel sheet with name like
            % 30 May 2019 Run1.xlsx
            % saves time, c, cpred along with estimated
            % parameter values and error values
            fn=home_dir+"\Results\"+string(datetime('now','Format','d MMM uuuu'))+OA(i)+".xlsx";
            if ~exist(fn,'file')
                xlswrite(fn,[iteration,j,mv,co,rpm,x,Er])
            else
                t=xlsread(fn);
                t=[t;[iteration,j,mv,co,rpm,x,Er]];
                xlswrite(fn,t);
            end
            
            % Plot the results
            % open new figure
            %                         figure
            %                         title("Run "+string(j)) % Set tile
            %                         plot(t,c,'ok',t,cp,'-r','LineWidth',2); % Plot
            %                         xlabel('Time, min')
            %                         ylabel('Concetration, mg/L')
            %                         legend('Exp',OA(i)+" Fit")
            %                         legend('boxoff') % Remove box around legend
            %                         drawnow % Plot immediately
        end
    end
    %                 catch
    %                     fprintf('\n Erroneous conditions \n')
    %                     disp([mv,co,rpm])
end
%             end
%         end
%     end
% end



function [Err,varargout]=HSDM(x,t,c,kF,nF,rhop,rp,mv,co,nr,r)

global nfcall

% Append number of function calls
nfcall=nfcall+1;

% Distribute x to kc and Ds
kc = x(1)*1e-7;
Ds = x(2)*1e-14;
% disp([kc,Ds])
% Compute qo in equilibrium with co
qo=kF*co^nF;

% Biot number ensure units are consistant
Bi = kc*rp*co/(Ds*rhop*qo);

% Distribution coefficient
DB = mv*qo/co;

% Dimensionless time
TB = (Ds*t*60)/(rp^2);

% Intial value vector
Yo = zeros(1,nr);

% Try to solve model function
try
    % Call ode function
    [T,Yo] = ode15s(@pde2ode,TB,Yo);
    
    
    % Compute Y average
    % for i=1:length(T)
    %     temp=fit(r',(r.^2.*Yo(i,:))','cubicspline') ;
    %     Ybar(i)=3*integral(@(x) temp(x)',0,1);
    % end
    Ybar=3*trapz(r,Yo'.*(r'.^2));
    
    % Find X using material balance equation
    X = 1-DB.*Ybar;
    
    % Conver to dimensional form
    cp = X*co;
    
    % Compute SSE
    Err = sum((c-cp').^2);
catch % catch error and set defaults
    cp=0*c;
    Err=inf;
    
end
% Pass variable output as cp value
varargout{1}=cp';

% disp([x,Err])
%% Model Function description
    function Yt = pde2ode(~,Y)
        
        % concert Y to row vector
        Y=Y';
        %         Ybar=fit(r',(r.^2.*Y)','smoothingspline');
        %         Ybar=3*integral(@(x) Ybar(x)',0,1);
        Ybar=3*trapz(r',r.^2.*Y);
        X = 1-DB*Ybar;     % caluclate Xbar
        
        % Find concetration in equilibrium with Y surface
        Xs = (Y(nr))^(2/nF)^(1/2);
        
        % Compute first derivative
        Yr=dss004(0.0,1,nr,Y);   % First order spatial derivative
        
        % Set boundary condition at surface
        Yr(nr)=Bi*(X-Xs);        % last spatial derivative at r=1
        
        % Set BC at center
        Yr(1)=0.0;               % first spatial derivative at r=0 is 0
        Yrr=dss004(0.0,1,nr,Yr); % Second order Spatial Derivative
        
        % Set time derivative vector
        Yt= [3*Yrr(1), Yrr(2:end)+(2./r(2:end)).*Yr(2:end)];
        Yt=Yt'; % Transpose to a column vector as needed by ode15s
    end

end

