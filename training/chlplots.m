
for i=1:3,
    per_rms=100*sqrt((est(:,i)-nom(:,i)).^2)./(ones(length(est),1)*(max(nom(:,i))-min(nom(:,i))));
    mean_per_rms=mean(per_rms)    
end;

figure;for i=1:3;subplot(3,1,i);plot([est(:,i) nom(:,i)]);end

figure;for i=1:3; subplot(3,1,i);plot(W(:,((21+2)*100+2)*i-1)); end
