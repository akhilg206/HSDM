clear all;clc;close all
format shortG

fn = {'Parameter Estimation','Simulation'};
[index,tf] = listdlg('PromptString',{'Select a task to perform',''},...
    'SelectionMode','single','ListString',fn);

Task = fn{index}

N = 3;

Rd = readtable('BG_3C.xlsx','Sheet',"Run Details");


Est_params_final = array2table(zeros(1,N+2));
Est_params_final.Properties.VariableNames = ["Run","Ds"+string(1:N),"SSE"];

rp = 6e-4;
isoFunc = @ExtendedSipps;
load EquiParams_BG.mat

nr=100;

co = 100*ones(1,N);

%% For Parameter estimation use this section

Inputs = inputdlg({'Carbon (0 or 1): ',...
    'Temperature',...
    ' pH','Save results'},'Give Inputs',[1,100],{'0','27','3','1'});
save = str2num(Inputs{4});
for C = str2num(Inputs{1})
    for T = str2num(Inputs{2})
        for pH = str2num(Inputs{3})
            
            %%Find the corresponding Run Number
            Run = find(Rd.Carbon == C & Rd.Temperature == T & Rd.pH == pH);
            Dose = Rd.Dose(Run);
            
            %%Get the equilibrium isotherm parameters
            eqp = equiParams.Parameters{equiParams.Carbon == C &...
                equiParams.Temperature == T &...
                equiParams.pH == pH};
            
            
            
            switch Task
                
                case 'Parameter Estimation'
                    
                    
                    %%Read The data for the corresponding run
                    Data = readtable('BG_3C.xlsx', 'Sheet', string(Run));
                    t = Data{:,1};
                    c = Data{:,2:end};
                    
                    
                    %%Calculate reference variables for non dimensionalizing
                    co = c(1,:);
                    qo = isoFunc(eqp,co);
                    
                    %%Create variables for search domain
                    [LB,UB,scale]=deal(ones(1*N,1),1e3*ones(1*N,1),1e-14);
                    
                    
                    %%Call EHO to estimate parameters
                    
                    funcCall = @(x) HSDM(N, x, t, Dose, rp, qo, isoFunc, eqp, c, nr, co,'scale',scale,'Display','off');
                    [x,fval] = ElephantHerd(funcCall, N, LB, UB,'Display','iter','Hybrid','on','functionTolerance',1e-3,'Display','iter');
                    clear HSDM
                    
                    [Er,t,c,cp,Yo,qbp] =  funcCall(x);
                    Est_params_final(end+1,:) = array2table([Run, x(1,:)*scale, Er]);
                    
                    %%Plot Results
                    mark=["d","s","o"];
                    line=["-",":","-."];
                    cole=["r","b","g"];
                    
                    figure
                    for i=1:N
                        hold on
                        plot(t, c(:,i),mark(i)+cole(i),...
                            t,cp(:,i),line(i)+cole(i))
                    end
                    
                    
                    if save
                        
                        writetable(Est_params_final,fullfile(cd,"Results","Est_Results.xlsx"),'Sheet',string(datetime('now','Format','d MMM uuuu')),'WriteMode','append')
                        
                    end
                    
                    
                    
                    
                    
                    
                    
                    %% USe this section for SImulation
                    
                case 'Simulation'
                    
                    %%Set the stop time
                    Inputs = inputdlg({'Simulation time (>0), min: ',...
                        ' Size of Time Steps',...
                        ' Ds to Simulate with (in m^2/s give as space separated values for multicomponent i.e. Ds1 Ds2 Ds3 )',...
                        ' Multiplier Value'},'Give Inputs',[1,100],{'1000','100','','1e14'});
                    if ~isempty(Inputs{2}) && ~isempty(Inputs{1}) && ~isempty(Inputs{3})
                        t = 0:str2num(Inputs{2}):str2num(Inputs{1})';
                        x = str2num(Inputs{3});
                        scale = str2num(Inputs{4});
                    else
                        
                        error('Invalid Inputs');
                        
                    end
                    
                    
                    %%Calculate reference variables for non dimensionalizing
                    qo = isoFunc(eqp,co);
                    Ds = x;
                    c = zeros(length(t),N);
                    [Er,t,c,cp,Yo,qbp,q] =  HSDM(N, Ds, t, Dose,     rp, qo, isoFunc, eqp, c, nr, co,'scale',1/scale,'Display','iter');
                    
                    % q is a table which can be accessed as q.CA for CA
                    % q.CA{:,:} will give nr rows matrix with  length of t columns 
                    
                    %% PLot Results
                    mark=["d","s","o"];
                    line=["-",":","-."];
                    cole=["r","b","g"];
                    
                    figure
                    for i=1:N
                        hold on
                        plot(t,cp(:,i),line(i)+cole(i))
                    end
                    
                    
                    if save
                        
                        writetable(Est_params_final,fullfile(cd,"Results","Sim_Results.xlsx"),'Sheet',string(datetime('now','Format','d MMM uuuu')),'WriteMode','append')
                        writetable(q,fullfile(cd,"Results","q_Results.xlsx"),'Sheet',string(datetime('now','Format','d MMM uuuu'))+"C"+C+"T"+T+"pH"+pH,'WriteMode','append')
                    end
                    
            end
        end
    end
end





