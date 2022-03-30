cleanfigure('targetResolution', 300);
saveas(gcf, [pwd '\figures\' figtitle '.png'])
matlab2tikz('showInfo', false, [pwd '\figures\' figtitle '.tikz'], 'height', figheight, 'width', '12cm')
% try
%     thesisfigs = 'C:\Users\geral\Desktop\Project thesis\LaTeX\StudentMaster-masternew\StudentMaster-master\thesis\figures\';
%     matlab2tikz('showInfo', false, [thesisfigs figtitle '.tikz'], 'height', figheight, 'width', '12cm')
% end