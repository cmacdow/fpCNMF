function W = SyntheticW(K,N,L,no_overlap_flag)
    W = zeros(N,K,L);    

    for k = 1:K
        if no_overlap_flag %force non overlapping motifs        
            M = floor(N/K); %number of pixels active per motif. equal across motifs
        else
            M = randi([ceil(N/2),N]); 
        end
        
        a = zeros(M,L);
        a(randi([1,M]),randi([1,L]))=1;    
        for m = 1:randi([ceil(M/2),M])
            a = a + circshift(a,[1,1]);
            a = a /max(a(:));
        end                
        W((k-1)*M+1:M*k,k,:) = a;  
    end   
end
                