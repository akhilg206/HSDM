function [x,varargout]=EHO(func,nvars,lb,ub,options)

arguments
func function_handle
nvars (1,1) double {mustBeInteger, mustBePositive, mustBeScalarOrEmpty, mustBeNonempty}
lb (:,1) double {mustBeNonempty, mustBeReal} = -1e50*ones(nvars,1)
ub (:,1) double {mustBeNonempty, mustBeReal} =  1e50*ones(nvars,1)
options.functionTolerance (1,1) double {mustBePositive, mustBeScalarOrEmpty, mustBeNonempty} = 1e-5
options.omegaStar (1,1) double {mustBeInteger, mustBePositive, mustBeScalarOrEmpty, mustBeNonempty} = 5*nvars
options.numMothers (1,1) double {mustBeInteger, mustBePositive, mustBeScalarOrEmpty, mustBeNonempty} = max(10,nvars^2)
options.numCalves (1,1) double {mustBeInteger, mustBePositive, mustBeScalarOrEmpty, mustBeNonempty} = max(2,2*nvars)
options.pValue (1,1) double {mustBeGreaterThan(options.pValue,0), mustBeLessThan(options.pValue, 1)} = 0.95
options.dbstate {ismember(options.dbstate,{'off','basic','advance','full'})} = 'off'
options.vectorized {ismember(options.vectorized,{'off','on'})} = 'off'
end


%% Search Domain variables
epsilon=abs(ub-lb);
sigma=log10(epsilon)/4;

%% create balnk elephants using createElephants function 

[ME,CE] = createElephants(options, nvars);
ME = initialize(ME, 'Mother', lb, ub, options.numMothers, nvars)

% 

% Mat=M;
% ME=repmat(M,1,nMother);
% ME=init(func,ME,nMother,lb,ub,epsilon);
% [gbest,id]=min([ME.curfval]);
% Mat=ME(id);
% 
% C=struct('realM',[],'curM',[],'curloc',zeros(nVars,1),'curfval',[],'bestloc',[],'bestfval',[]);
% CE=repmat(C,nMother,nCalf);
% CE=init(func,CE,nCalf,lb,ub,epsilon,nMothers);
% 
% 
% while Mat.Omega<Omegastar 
%    
%     
%     
% end
% 
% end
% function E=init(func,E,nE,lb,ub,epsilon,nAE)
% if nargin<4
%     init_loc=lb+epsilon.*rand(size([E.curloc]));
%     for i=1:nE
%         E(i).SN=i;
%         E(i).curloc=init_loc(:,i);
%         E(i).bestloc=E(i).curloc;
%         try
%             E(i).curfval=func(E(i).curloc');
%         catch userFcn_ME
%             msg = message('Objective function evaluation failed');
%             psw_ME = MException(msg.Identifier, getString(msg));
%             userFcn_ME = addCause(userFcn_ME, psw_ME);
%             rethrow(userFcn_ME)
%         end
%         if numel(E(i).curfval) ~= 1
%             error(message('Objective function should return scalar'));
%         end
%         E(i).bestfval=E(i).curfval;
%     end
% else
%     
%     for j=1:nAE
%         init_loc=lb+epsilon.*rand(size([E(j,:).curloc]));
%         for i=1:nE
%             E(i).realM=j;
%             E(i).curM=E(i).realM;
%             E(i).curloc=init_loc(:,i);
%             E(i).bestloc=E(i).curloc;
%             try
%                 E(i).curfval=func(E(i).curloc');
%             catch userFcn_ME
%                 msg = message('Objective function evaluation failed');
%                 psw_ME = MException(msg.Identifier, getString(msg));
%                 userFcn_ME = addCause(userFcn_ME, psw_ME);
%                 rethrow(userFcn_ME)
%             end
%             if numel(E(i).curfval) ~= 1
%                 error(message('Objective function should return scalar'));
%             end
%             E(i).bestfval=E(i).curfval;
%         end
%     end
% end
% end

%% Funtion to create blank elephants at the start
    function [M,C] = createElephants(options, nvars)
        
        M = cell2table(repmat({[],[zeros(nvars,1)], [], [zeros(nvars,1)], [],0},options.numMothers,1));
        M.Properties.VariableNames = {'SNo','curloc','curfval','bestloc','bestfval','Omega'};
        M.SNo(1:options.numMothers) = num2cell(1:options.numMothers);
        
        C = cell2table(repmat({[],[zeros(nvars,1)], [], [zeros(nvars,1)], []},options.numCalves,1));
        C.Properties.VariableNames = {'SNo','curloc','curfval','bestloc','bestfval'};
        C = repmat({C},1,options.numMothers);
    end

%% FUnction to Initialize Elephants
    function init_E = initialize(Elephants, Group, lb, ub, nE, nvar)

        switch Group
           
            case 'Mother'
                
                Elephants.curloc(:) = table2cell(table((lb+rand(nvar,nE).^1.*(ub-lb))'));
                Elephants.bestloc(:) = Elephants.curloc(:);
                
                for i = 1:nE
                   
                    
                    
                end
                
            case 'Calf'
            
            
        end
%         
%         if nargin<4
%     init_loc=lb+epsilon.*rand(size([E.curloc]));
%     for i=1:nE
%         E(i).SN=i;
%         E(i).curloc=init_loc(:,i);
%         E(i).bestloc=E(i).curloc;
%         try
%             E(i).curfval=func(E(i).curloc');
%         catch userFcn_ME
%             msg = message('Objective function evaluation failed');
%             psw_ME = MException(msg.Identifier, getString(msg));
%             userFcn_ME = addCause(userFcn_ME, psw_ME);
%             rethrow(userFcn_ME)
%         end
%         if numel(E(i).curfval) ~= 1
%             error(message('Objective function should return scalar'));
%         end
%         E(i).bestfval=E(i).curfval;
%     end
% else
%     
%     for j=1:nAE
%         init_loc=lb+epsilon.*rand(size([E(j,:).curloc]));
%         for i=1:nE
%             E(i).realM=j;
%             E(i).curM=E(i).realM;
%             E(i).curloc=init_loc(:,i);
%             E(i).bestloc=E(i).curloc;
%             try
%                 E(i).curfval=func(E(i).curloc');
%             catch userFcn_ME
%                 msg = message('Objective function evaluation failed');
%                 psw_ME = MException(msg.Identifier, getString(msg));
%                 userFcn_ME = addCause(userFcn_ME, psw_ME);
%                 rethrow(userFcn_ME)
%             end
%             if numel(E(i).curfval) ~= 1
%                 error(message('Objective function should return scalar'));
%             end
%             E(i).bestfval=E(i).curfval;
%         end
%     end
end
end
