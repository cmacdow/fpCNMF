function H = SyntheticH(K,T,L,frac_active_time,no_overlap_flag)
        
H = zeros(K,T);
if no_overlap_flag %force non overlapping motifs
    timepoints = 1:L+ceil(L/2):T;
    num_splits = ceil(numel(timepoints)/K);
    for i = 1:size(H,1) %use as many motifs as you can and ~equally distribute occurances
        if numel(timepoints)>=num_splits
            idx = timepoints(randperm(numel(timepoints),num_splits));
            H(i,idx) = 1 * (rand(1,numel(idx))+0.5); %randomly weight
            timepoints(ismember(timepoints,idx))=[]; %remove the used indices
        else %populate with the remaining
            H(i,timepoints) = 1 * (rand(1,numel(timepoints))+0.5); %randomly weight
        end
    end
    H = H(randperm(size(H,1),size(H,1)),:);%randomly shuffle so last motif is not systematically removed
else %randomly distribute H weightings with a highly skewed distribution
    H = exprnd(1,K,T) .* binornd(10,frac_active_time, K,T);
end

%smooth each transient with gaussian
for i = 1:size(H,1)
    H(i,:) = smoothdata(H(i,:),'gaussian',round(L/3)); %smooth with a gaussian to make less binary. Keep <0.5L to maintain no_overlap condition    
end

end %function end