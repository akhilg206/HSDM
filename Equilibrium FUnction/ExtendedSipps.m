function [out,varargout] = ExtendedSipps(parameters,ce)

[Qm,k,n] = deal(parameters(1,:),parameters(2,:),parameters(3,:));

out = Qm.*(ce).^n./(1+sum(k.*ce.^n,2));


end