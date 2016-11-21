function limo_show(Y,LIMO)


% dimensions of Y are supposed to be : channels, time, trials

nb_conditions = 2;
nb_continuous = 0;
X = zeros(size(Y,3),nb_conditions+nb_continuous+1);

X(1:50,1) = 1;
X(51:end,2) = 1;

X(:,end) = 1;

LIMO.design.nb_conditions = nb_conditions;
LIMO.design.nb_continuous = nb_continuous;
LIMO.design.X = X;
LIMO.design.nb_items = sum(X(:,1:nb_conditions));
LIMO.design.conditions = 1:nb_conditions;
LIMO.design.name = 'Test';
LIMO.electrode = 27;

[LIMO ] = refreshmodel(Y,LIMO);
refreshfig(Y,LIMO);

function refreshfig(Y,LIMO)
electrode = LIMO.electrode;
figure(4536);

set(gcf,'userdata',struct('LIMO',LIMO,'Y',Y));

colormap('gray')
subplot(1,4,1)
title('Y')
h = imagesc(squeeze(Y(electrode,:,:))');
set(h,'tag','sY')
view(2);shading interp
subplot(1,4,2)
title('B')
h = imagesc(reshape(LIMO.Betas(electrode,:,:),size(Y,2),size(LIMO.design.X,2)));
set(h,'tag','sB');
view(2);shading interp
subplot(1,4,3)
title('X')
h = imagesc(LIMO.design.X);
set(h,'buttondownfcn',@update_design)
set(h,'tag','sX');
view(2);shading interp
subplot(1,4,4)
title('Res')
h = imagesc(squeeze(LIMO.Res(electrode,:,:))');
set(h,'tag','sR');
view(2);shading interp


return

function update_design(hObject,evt)

ud = get(gcf,'userdata');
LIMO = ud.LIMO;
Y = ud.Y;
sY = findobj('tag','sY');
sB = findobj('tag','sB');
sR = findobj('tag','sR');
sX = findobj('tag','sX');

X = get(sX,'CData');
pos = get(gca,'currentpoint');
x = pos(1,1);
y = pos(1,2);

Xx = round(x);
Yy = round(y);
state = 1-X(Yy,Xx);
newstate = 1-state*ones(1,LIMO.design.nb_conditions);
newstate(:,Xx) = state;
X(Yy,1:LIMO.design.nb_conditions) = newstate;
set(sX,'CData',X);
refreshdata;

LIMO.design.X = X;

LIMO = refreshmodel(Y,LIMO);
refreshfig(Y,LIMO)

function LIMO = refreshmodel(Y,LIMO)
electrode = LIMO.electrode;

for ielec = 1:numel(electrode)
    warning off; fprintf('analyzing electrode %g',electrode(ielec)); disp(' ');
    model = limo_glm(squeeze(Y(electrode(ielec),:,:))',LIMO); warning on;
    
    % update the LIMO.mat (do it only once)
    if ielec == 1
        LIMO.model.model_df = model.df;
        if LIMO.design.nb_conditions ~=0
            LIMO.model.conditions_df  = model.univariate.conditions.df;
        end
        if LIMO.design.nb_continuous ~=0
            LIMO.model.continuous_df  = model.univariate.continuous.df;
        end
    end
    
    % update the files to be stored on the disk
    fitted_data = LIMO.design.X*model.betas;
    LIMO.Yhat(electrode(ielec),:,:) = fitted_data';
    LIMO.Res(electrode(ielec),:,:)  = squeeze(Y(electrode(ielec),:,:)) - fitted_data'; clear fitted_data
    LIMO.R2(electrode(ielec),:,1) = model.R2_univariate; R2(electrode(ielec),:,2) = model.F; R2(electrode(ielec),:,3) = model.p;
    LIMO.Betas(electrode(ielec),:,:,1) = model.betas';
    
    if LIMO.design.nb_conditions ~=0
        LIMO.Condition_effect(electrode(ielec),:,1) = model.univariate.conditions.F;
        LIMO.Condition_effect(electrode(ielec),:,2) = model.univariate.conditions.p;
    end
    
    if LIMO.design.nb_continuous ~=0
        for i=1:LIMO.design.nb_continuous
            LIMO.Continuous(electrode(ielec),:,i,1) = model.univariate.continuous.F(i,:);
            LIMO.Continuous(electrode(ielec),:,i,2) = model.univariate.continuous.p(i,:);
        end
    end
end




