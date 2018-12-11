clear all
close all
%% create path
delta_t = 1e-1;

ay_cell_model();
[~,Ak,Bk,Qk]=ay_path_simulator();
load('Path_X');
load('Cell_Field_Model');
[MS,PS] = ay_cell_simulator(X,CellModel,delta_t);

% create cell
figure(1)
load('Path_X');
load('Cell_Field_Model');
load('MS_PS');
%% plot the result
figure(1)
subplot(1,3,1)
for m=1:length(CellModel)
    for k=1:length(CellModel{m})
        xa=CellModel{m}(k).M(1);
        ya=CellModel{m}(k).M(2);
        plot(xa,ya,'o','LineWidth',4,'MarkerSize',18,'Color',[0.2 0.1+0.04*m    0.2+0.02*m ]);hold on;
        alpha 0.1
        hold on
    end
    text(xa-3,ya-3,num2str(m),'Color',[1  0 0],'FontSize',14);
end
plot(X(:,1),X(:,2),'y-','LineWidth',3);
hold off
xlabel('X')
ylabel('Y')
title('Yellow line shows the path, and circles show center of cells receptive field')
%% simulate path
subplot(1,3,2)
delta_t = 1e-1;
imagesc(MS(:,2:end-2))
xlabel('Cell')
ylabel('Time Index')

%% decoder
%% GMM merge parameters
merge_c = 0.1;
drop_c  = 0.99;
%% Main loop
post_gmm(1).w= 1;
post_gmm(1).m= [MS(1,end-1); MS(1,end)];
post_gmm(1).s= Qk;
DecoderMean  = zeros(size(MS,1),2);
% Set up figure properties:
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
set(gcf, 'Toolbar', 'none', 'Menu', 'none');
subplot(1,3,3)
for i=1:size(MS,1)
    %% Sk
    Sk = MS(i,2:end-2);
    tic
    %% prior
    pre_gmm    = post_gmm;
    %% one-step 
    one_step   = ay_gmm_one_step(pre_gmm,Ak,0*Bk,1.4*Qk);
    %% filter
    post_gmm   = ay_gmm_approx_fast(one_step,Sk,Mk,CellModel,FieldModel,delta_t,drop_c,2);
    %% merge up to 1
    [i length(post_gmm)]
    if length(post_gmm) > 1
        post_gmm = ay_gmm_merge_alpha_optimized(post_gmm,post_gmm,merge_c);
    end
    length(post_gmm)
    time_proc=toc;
    %% calculate mean of estimate
    DecoderMean(i,:) = ay_gmm_mean(post_gmm);
    
    plot(X(:,1),X(:,2),'Color',[0.9 0.9 0.9],'LineWidth',2);hold on
    plot(X(1:i,1),X(1:i,2),'b-');
    plot(DecoderMean(1:i,1),DecoderMean(1:i,2),'r','LineWidth',2);
    if sum(Sk)
        plot(X(i,1),X(i,2),'go','LineWidth',3,'MarkerSize',12);
    end
    legend('Whole Path','Path','Decode','Location','Best')
    xlabel('X')
    ylabel('Y')
    title(['Time Index:' num2str(i) ', Time: ' num2str(time_proc*1000)  ' msec'   ])
    
    hold off
    pause(0.01)
end
save('DecoderMean','DecoderMean');