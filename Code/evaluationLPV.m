close all
xxyy = 'xy';
figure(fig)
fig = fig + 1;

tplot = [];
leg = [];
normeta = zeros(N,length(t));
for k = 1:length(t)
    for ii = 1:N
        normeta(ii, k) = sqrt( eta(n.y*(ii-1)+1:n.y*ii,k)'*eta(n.y*(ii-1)+1:n.y*ii,k));
        iinormeta(k) = sqrt(eta(:,k)'*eta(:,k));
        networknormeta(esti+1,k) = sqrt(eta(:,k)'*eta(:,k));
        grid on;
    end
end
subplot(7,1,[1 2 3 4])
for ii = 1:N
    plot(t,normeta(ii,:))
    hold on
    leg = [leg; 'Agent ' num2str(ii)];
end
grid on
legend(leg)
% title(['Trigger events: ' num2str( sum(sum(trigcount))) ' Trigger rate: ' num2str(sum(sum(trigcount))/N/length(t))]) 
ylabel('Formation error norm')

trigcountplot(esti+1,:) = sum(trigcount);
subplot(7,1,[6 7])
for k = 1:length(t)
    if sum(trigcount(:,k)) ~= 0
        tplot = [tplot, [t(k) ;sum(trigcount(:,k))]];
    end
end

stem(t,sum(trigcount)); %,'k.','LineStyle','none')
axis([0 t(end) 0 3])
grid on
xlabel('t')
ylabel('Trigger events')

if savefigs
    figtitle = [Estimator(esti+1,:) '_EtaTrigger'];
    cleanfigure('targetResolution', 300);
    saveas(gcf, [figtitle '.png'])
    matlab2tikz('showInfo', false, [figtitle '.tikz'], 'height', '3cm', 'width', '7cm')
end

%%
figure(fig)
fig = fig + 1;
for xy = 1:2
    subplot(7,1,[4*xy-3 4*xy-2 4*xy-1])
    plot(t,pos(xy,:))
    hold on
    grid on
    axis([0 t(end) -1 1])
%     title([xxyy(xy) ' - Direction'])
    plot(t,pos(xy+2,:))
    plot(t,pos(xy+4,:))
    
    ylabel(xxyy(xy))
end
    xlabel('t')

if savefigs
    figtitle = [Estimator(esti+1,:) '_XYPosition'];
    cleanfigure('targetResolution', 300);
    saveas(gcf, [figtitle '.png'])
    matlab2tikz('showInfo', false, [figtitle '.tikz'], 'height', '3cm', 'width', '7cm')
end
%%
figure(fig)
fig = fig + 1;
for ii = 1:N
    subplot(8,3,[ii ii+3])
    plot(pos(2*ii-1,:),pos(2*ii,:))
    hold on
    grid on
    title(['Agent ' num2str(ii)])
    xlabel('x')
    ylabel('y')
    axis(max(max(abs(pos)))*1.1*[-1 1 -1 1]);
end
for xy = 1:2
    subplot(8,2,[xy+6 xy+8])
    plot(t,pos(xy,:))
    hold on
    grid on
    title([xxyy(xy) ' - Direction'])
    plot(t,pos(xy+2,:))
    plot(t,pos(xy+4,:))
    xlabel('t')
    ylabel(xxyy(xy))
    legend('Agent 1','Agent 2','Agent 3')
end
subplot(8,1,[7 8])
stem(t,sum(trigcount))%,'k.','LineStyle','none')
axis([0 t(end) 0 3])
grid on
xlabel('t')
ylabel('Trigger events')
title('Simultaneous trigger events')
if savefigs
    figtitle = [Estimator(esti+1,:) '_Performance'];
    cleanfigure('targetResolution', 300);
    saveas(gcf, [figtitle '.png'])
end
%%
Results(esti+2,:) = {Estimator(esti+1,:), sum(iinormeta)/length(t), sum(sum(trigcount)), sum(sum(trigcount))/N/length(t)};