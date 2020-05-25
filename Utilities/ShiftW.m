function [W,shiftvec] = ShiftW(W,shift_in)
    %adapted from from seqNMF toolbox
    % shift factors by center of mass
    if nargin <2 %if the desired shift is not specified then shift by center of mass
        shift_in = NaN(1,size(W,2));
    end
        
    
    % get size of W and H
    [N,K,L] = size(W);
    
    if L>1 % if L=1, no room to shift
    
    center = max(floor(L/2),1); % middle bin
    
    % pad with zeros, for use with circshift. data is zeropadded within seqNMF, so don't need to pad H
    Wpad = cat(3,zeros(N,K,L),W,zeros(N,K,L));
    shiftvec = zeros(K);
    for k = 1:K
        % compute center of mass
        temp = nansum(squeeze(W(:,k,:)),1);
        cmass = max(floor(nansum(temp.*(1:length(temp)))/nansum(temp)),1); 
        if isnan(shift_in(k)) %Use center of mass
            Wpad(:,k,:) = circshift(squeeze(Wpad(:,k,:)),[0,center-cmass]);
        else %Use the input shift value
            Wpad(:,k,:) = circshift(squeeze(Wpad(:,k,:)),[0,shift_in(k)]);
        end
        shiftvec(k) = center-cmass; %outputs the shifts that you used
    end
 

    % undo zero pad
    W = Wpad(:,:,(L+1):(end-L));

    end