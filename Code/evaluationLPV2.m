close all
xxyy = 'xy';
figure(fig)
fig = fig + 1;
leg = [];
plotesti = {'-',':','--'};
for ii = 1:3
    subplot(18,1,1:6)
    plot(t,networknormeta(ii,:),['k' plotesti{ii}])
    hold on
    leg = [leg; Estimator(ii,:)];
    subplot(18,1,4*ii+4:4*ii+6)
    stem(t,trigcountplot(ii,:),'k'); %,'k.','LineStyle','none')
    hold on
%     title(Estimator(ii,:))
    ylabel(Estimator(ii,:))
    axis([0 t(end) 0 3])
    grid on
end

subplot(18,1,1:6)
grid on
legend(leg)
ylabel('Formation error norm')

% subplot(16,1,10:12)
% ylabel('Trigger events')

subplot(18,1,16:18)
xlabel('t')

if savefigs
    figtitle = 'EtaTrigger';
    cleanfigure('targetResolution', 300);
    saveas(gcf, [figtitle '.png'])
    matlab2tikz('showInfo', false, [figtitle '.tikz'], 'height', '5cm', 'width', '7cm')
end
%%
figure(fig)
fig = fig + 1;
plotesti = {'-',':','--'};
% colresti = {'b','g','
for ii = 1:3
    for xy = 1:2
        subplot(11,1,6*xy-5:6*xy-1)
        plot(t,posesti(xy,:,ii),plotesti{ii},'Color',[0 0.4470 0.7410])
        hold on
        grid on
        axis([0 t(end) -1 1])
        %     title([xxyy(xy) ' - Direction'])
        plot(t,posesti(xy+2,:,ii),plotesti{ii},'Color',[0.8500 0.3250 0.0980])
        plot(t,posesti(xy+4,:,ii),plotesti{ii},'Color',[0.9290 0.6940 0.1250])
        
        ylabel(xxyy(xy))
    end
end
    xlabel('t')
    subplot(11,1,1:5)
    legend('Agent 1','Agent 2','Agent 3')

if savefigs
    figtitle = 'XYPosition';
    cleanfigure('targetResolution', 300);
    saveas(gcf, [figtitle '.png'])
    matlab2tikz('showInfo', false, [figtitle '.tikz'], 'height', '3.5cm', 'width', '7cm')
end