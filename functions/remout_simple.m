function solset_new = remout_simple(solset,burnind,nouts)
%REMOUT_SIMPLE Naive outlier removal from a DE/DEMC solset object. 
%   solset_new = REMOUT_SIMPLE(solset,burnind,nouts) removes the outliers 
%   from a solset object returned by the optimization algorithm. 
%   burnind is the burn-in index, typically half the number of generations. 
%   nouts is the number of outliers and the chains with the nout highest 
%   mean objective function values (after burnin) are removed. 
    
    % Rank chains by the mean objective function value after burnin
    means = mean(solset.Flist(burnind:end,:),1);
    [~, inds] = sort(means);
    goodinds = sort(inds(1:end-nouts));
    
    % Slice out the 'good' chains
    solset_new.F = solset.F(goodinds);
    solset_new.Flist = solset.Flist(:,goodinds);
    solset_new.X = solset.X(goodinds,:);
    solset_new.Xlist = solset.Xlist(:,goodinds,:);

end