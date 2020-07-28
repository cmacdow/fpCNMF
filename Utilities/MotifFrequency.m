function [pwr,freq] = MotifFrequency(X,fps)

%convert NaN (periods of no motif activity) to zero
X(isnan(X))=0;

pwr = cell(1,size(X,1));
for i = 1:size(X,1)
    temp = X(i,:);
    N = length(temp);
    xdft = pmtm(temp); %pmtm used to be fft
    xdft = xdft(1:N/2+1);
    psdx = (1/(fps*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:fps/length(temp):fps/2;
    pwr{i} = psdx';
%     figure; hold on; plot(freq,pwr{i});
%     grid on
%     xlabel('Frequency (Hz)')
%     ylabel('Power (dB)')
%     xlim([0 2])
%     axis square;
%     pause()
   
end


end