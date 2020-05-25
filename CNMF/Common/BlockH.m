function block_h = BlockH(h,L)
    [K,T] = size(h);
    block_h = zeros(L*K,T); 
    for lag = 0:L-1
       temp = circshift([zeros(K,L),h,zeros(K,L)],lag,2); %pad and shift
       block_h(1+K*lag:K*(lag+1),:)= temp(:,L+1:size(temp,2)-L);%remove pad and save
    end
end %function