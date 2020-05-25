function FitKandL(X,K_range,L_range)
%we want to choose a K and L that maximize out fit to the data
L_range = [5:5:25];
K_range = [1:7];
X = data{1};
sweep = combvec(K_range,L_range);
sweep = K_range;
sweep =L_range;
%%
diss = [];
pev = [];
for i = 1:size(sweep,2)
    fprintf('\n%d of %d',i,size(sweep,2))
    for rep = 1:5
%         [W, H] = fpCNMF(X,'K',3,'L',sweep(1,i),'non_penalized_iter',10,'penalized_iter',0,'speedy',1);            
        [W, H] = fpCNMF(X,'K',sweep(1,i),'L',50,'non_penalized_iter',10,'penalized_iter',0,'speedy',1);  
        if rep>1
           diss(i,rep-1) = DissimilarityX(H,W,H_prev,W_prev); 
        end
        W_prev = W;
        H_prev = H;
        stats = CNMFStats(W,H,X,0);
        pev(i,rep) = stats.pev;
    end
end

%pev near the plateau and lower variance in the pev

%APPROXIMATE L By looking at the duration of bursts in the data

figure; hold on; plot(K_range,nanmean(diss,2)); plot(K_range,nanmean(pev,2));
figure; hold on; plot(L_range,nanmean(diss,2)); plot(L_range,nanmean(pev,2));
%%

figure; hold on;
plot(sweep(1,:),pev)
plot(sweep(2,:),pev)


%sweep; we don't care how sparsely we fit the data. 


%%TOMORROW WRITE THE MULTIPLICATIVE UPDATE ALGORITHM FOR YOUR VERSION

%plot the explained variance

 [W, H] = fpCNMF(X,'K',5,'L',20,'non_penalized_iter',10,'penalized_iter',0,'speedy',1);
 VisualizeData(tensor_convolve(W,H),W,H)
 