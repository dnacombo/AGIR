function mylimo_show_design_matrix(LIMO)


figure('Name','design matrix','Color','w','NumberTitle','off')
Xdisplay = LIMO.design.X;
if   LIMO.design.nb_continuous
    REGdisplay = LIMO.design.X(:, LIMO.design.nb_conditions+1:size(LIMO.design.X,2)-LIMO.design.intercept);
    REGdisplay = REGdisplay + max(abs(min(REGdisplay)));
    Xdisplay(:, LIMO.design.nb_conditions+1:size(LIMO.design.X,2)-LIMO.design.intercept) = REGdisplay ./ max(max(REGdisplay));
end
imagesc(Xdisplay); colormap('gray'); drawnow;
title('Design matrix','fontsize',18); xlabel('regressors','fontsize',18);ylabel('trials','fontsize',18);
set(gca,'XTick',1:size( LIMO.design.X,2),'fontsize',18)
drawnow


return