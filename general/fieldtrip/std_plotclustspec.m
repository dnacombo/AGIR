function std_plotclustspec(STUDY,ALLEEG,clusts)

STUDY = std_readspec(STUDY,ALLEEG,'clusters',clusts);

aploter = [];
for i = 1:numel(clusts)
    aploter(i,:) = mean(STUDY.cluster(clusts(i)).specdata{1},2)';
end
set(gca,'colororder',varycolor(numel(clusts)))
hold on;
plot(STUDY.cluster(clusts(1)).specfreqs,aploter');
hold off;
xlabel('Frequency (Hz)')
ylabel('Power (10 x log_1_0({\mu}V^2)','interpreter','tex')
xlim([min(STUDY.cluster(clusts(1)).specfreqs) max(STUDY.cluster(clusts(1)).specfreqs)]);

% legstr = {};
% for i = 1:numel(clusts)
%     legstr = [legstr {['Cluster ' num2str(clusts(i)) ]}];
% end
legstr = num2str(clusts(:)-1);
legend(legstr);

