clear all
close all
%% generate model
load('params_epc0_bond_mdec0906-8000-8_2-1.mat')
temp_ind = 0;
Xm=[];
Sm=[];
Wm=[]; 
temp_gmm =[];
for i  = 1:1
    ind = 0;
    temp = [];
    j_ind =0;
    for j=1:M
            ind = ind +1;
         %   sx  = sp_prmPstMd(ind) * pdf('normal',sp_prmPstMd(ind+1),sp_prmPstMd(ind+1),sqrt(sp_prmPstMd(ind+2)));
            if (sp_prmPstMd(ind) > 1e-9) && (sp_prmPstMd(ind+2)<1e3) && (sp_prmPstMd(ind+2)>1e-5)
                j_ind = j_ind +1;
                temp(j_ind).W = sp_prmPstMd(ind);
                % mark domain
                temp(j_ind).Mm = mk_prmPstMd_u(j,:);
                temp(j_ind).Cm = squeeze(mk_prmPstMd_Sg(j,:,:));


                % point process
                ind = ind + 1;
                temp(j_ind).Mp = sp_prmPstMd(ind);
                
                ind = ind + 1;
                temp(j_ind).Cp = sp_prmPstMd(ind);
                
                Xm = [Xm;temp(j_ind).Mp];
                Sm = [Sm;temp(j_ind).Cp];
                Wm = [Wm;temp(j_ind).W];

                temp_ind = temp_ind + 1;
                
                temp_gmm(j_ind).w = temp(j_ind).W;
                temp_gmm(j_ind).s = temp(j_ind).Cp;
                temp_gmm(j_ind).m = temp(j_ind).Mp;
            else
                ind = ind + 2;
            end
 
    end
    CellModel{i}= temp;
end
for i=1:length(temp_gmm)
    temp_gmm(i).w=temp_gmm(i).w/sum(Wm);
end

temp_gmm   = ay_gmm_merge_alpha_optimized(temp_gmm,temp_gmm,0.001);
temp_all   = [];
for i=1:length(temp_gmm)
    temp_all(i).Cp = temp_gmm(i).s;
    temp_all(i).Mp = temp_gmm(i).m;
    temp_all(i).W  = temp_gmm(i).w*sum(Wm);
end
FieldModel = temp_all;
delta_t    = 1e-3;

%% define Ak,Bk, Qk
Ak = 1;
Bk = 0;
Qk = 1e-3;

% create cell
load('MS_PS_1');
%% decoder
start_ind =  192760; % 220297151454;%96380;
end_ind   =  220297; %165222;%123917;
%% GMM merge parameters
merge_c = 0.01;
drop_c  = 0.99;
%% Main loop
post_gmm(1).w = 1;
post_gmm(1).m = MS(start_ind,1);
post_gmm(1).s = 10*Qk;
DecoderMean   = zeros(end_ind-start_ind+1,1);
Xs = [-10:0.01:10];
DecoderPdf   = zeros(end_ind-start_ind+1,length(Xs));
DecoderLL    = zeros(end_ind-start_ind+1,length(Xs));
MixtureNo    = zeros(end_ind-start_ind+1,1);

% Set up figure properties:
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.

%% e stands for exact
e_post       = ay_gmm_posterior_1d(post_gmm,Xs);
eDecoderPdf  = DecoderPdf;
eDecoderMean = DecoderMean;

sum_num  = 0;
sum_time = 0;
sum_mix  = 0;
sum_spk_mix = 0;
sum_spk     = 0;
ML_NOSPIKE = ay_gmm_ll(0,[0 0 0 0],CellModel,FieldModel,delta_t,Xs);
fig_ind  = 1;
sample_ind = 0;
for i=start_ind:end_ind            %size(MS,1)
    [i length(post_gmm)]
    sample_ind = sample_ind + 1; 
    %% Sk
    Sk = MS(i,2);
    Mk = MS(i,3:end);
    tic
    %% prior
    pre_gmm    = post_gmm;
    %% one-step 
    one_step   = ay_gmm_one_step(pre_gmm,Ak,0*Bk,Qk);
    %% filter
    % here, we don't check posteior covariance matrix to be PSD (this is correct,as we have changed the ground truth using merging)
    post_gmm   = ay_gmm_approx_fast(one_step,Sk,Mk,CellModel,FieldModel,delta_t,drop_c,2);
    %% we can use this if we want to check each covariance matrix to be PSD
    % post_gmm   = ay_gmm_approx(one_step,Sk,Mk,CellModel,FieldModel,delta_t,drop_c);
    
    px_gmm     = post_gmm;
    %% merge up to 1
    if length(post_gmm) > 1
        [post_gmm,xxx_gmm] = ay_gmm_merge_alpha_optimized(post_gmm,post_gmm,merge_c);
    end
    time_proc=toc;
    sum_num  = sum_num+1;
    sum_time = sum_time + time_proc;
    sum_mix  = sum_mix  + length(post_gmm);
    
    
    %% Exact Method
%     e_one_step = ay_exact_one_step(e_post,Ak,0*Bk,Qk,Xs);
%     e_ll       = ay_gmm_ll(Sk,Mk,CellModel,FieldModel,delta_t,Xs);
%     e_post     = e_ll.*e_one_step;
%     e_post     = e_post/sum(e_post);
%     
%     eDecoderPdf(sample_ind,:) = e_post;
%     eDecoderMean(sample_ind,:) = e_post*Xs';

    %% calculate mean of estimate
    [DecoderPdf(sample_ind,:),DecoderMean(sample_ind,:)]=ay_gmm_posterior_1d(post_gmm,Xs);
    MixtureNo(sample_ind) = length(post_gmm);
  %  if Sk
   if i==end_ind
        sum_spk_mix = sum_spk_mix + length(post_gmm);
        sum_spk     = sum_spk +1;
    
        subplot(4,2,5)
        plot(delta_t *(1:sample_ind),DecoderMean(1:sample_ind,:),'r','LineWidth',2);hold on
        plot(delta_t *(1:sample_ind),eDecoderMean(1:sample_ind,:),'b','LineWidth',2);
        plot(delta_t *(1:sample_ind),MS(start_ind:i,1),'k','LineWidth',2);
        if sum(Sk)
            plot(delta_t *sample_ind,MS(i,1),'go','LineWidth',3,'MarkerSize',12);
        end
        hold off
        legend('GMM Mean','Exact Mean','Rat Position','Location','northwest')
        xlabel('Time')
        ylabel('X')
        title(['Time Index:' num2str(i) ', Avg Proc. Time: ' num2str(sum_time*1000/sum_num)  ' msec, Avg No. of Mix.:' num2str(sum_mix/sum_num) '(' num2str(sum_spk_mix/sum_spk) ')' ]);
        ylim([-8 8])

        subplot(4,1,1)
      %  imagesc(delta_t *(1:i),Xs,sqrt(DecoderPdf(1:i,:)'),[eps*eps  0.025]);
        imagesc(delta_t *(1:sample_ind),Xs,sqrt(DecoderPdf(1:sample_ind,:)'),[1e-3 1e-1]);
        hold on
        plot(delta_t *(1:sample_ind),DecoderMean(1:sample_ind,:),'r','LineWidth',2);
        plot(delta_t *(1:sample_ind),MS(start_ind:i,1),'k','LineWidth',2);hold off
        ylim([-8 8])
        colormap parula
        ylabel('X')
        xlabel('Time')
        title('GMM Post. Distribution')
        
        subplot(4,1,2)
        imagesc(delta_t *(1:sample_ind),Xs,sqrt(eDecoderPdf(1:sample_ind,:)'),[1e-3 1e-1]);
        hold on
        plot(delta_t *(1:sample_ind),eDecoderMean(1:sample_ind,:),'r','LineWidth',2);
        plot(delta_t *(1:sample_ind),MS(start_ind:i,1),'k','LineWidth',2);hold off
        ylim([-8 8])
        colormap parula
        ylabel('X')
        xlabel('Time')
        title('Exact Post. Distribution')
        
        subplot(4,3,10)
        plot(Xs,e_ll/max(e_ll),'k','LineWidth',2);hold on;
        temp = ay_gmm_posterior_1d(one_step,Xs);
        plot(Xs,temp/max(temp),'b','linewidth',2);
        plot(Xs,DecoderPdf(sample_ind,:)/max(DecoderPdf(sample_ind,:)),'r','linewidth',2);hold on;
        hold off
        title([ 'Spike Time, Pre:' num2str(length(px_gmm))  ' - Post:' num2str(length(post_gmm))] )
        legend('Likelihood','one-step','Decoder','Location','south','Orientation','horizontal')
        xlabel('X')
        xlim([-8 8])
        
        subplot(4,3,11)
        plot(Xs,e_ll/max(e_ll),'k','LineWidth',2);hold on;
        plot(Xs,e_one_step/max(e_one_step),'b','linewidth',2)
        plot(Xs,eDecoderPdf(sample_ind,:)/max(eDecoderPdf(sample_ind,:)),'r','linewidth',2);
        hold off
        title(['Spike Time Exact Post'] )
        legend('Likelihood','one-step','Decoder','Location','south','Orientation','horizontal')
        xlabel('X')
        xlim([-8 8])
        
        subplot(4,3,12)
        plot(Xs,DecoderPdf(sample_ind,:)/max(DecoderPdf(sample_ind,:)),'r','linewidth',2);
        hold on;
        plot(Xs,eDecoderPdf(sample_ind,:)/max(eDecoderPdf(sample_ind,:)),'b','linewidth',2);
        hold off;
        legend('GMM','Exact','Location','south','Orientation','horizontal')
        title('GMM & Exact')
        xlabel('X')
        xlim([-8 8])
        
        
        subplot(4,2,6)
        plot(Xs,ML_NOSPIKE,'linewidth',2)
        xlabel('X')
        title('Likelihood Non-Spike')
        xlim([-8 8])
        
        hold off
        pause(0.1)
   
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 32 24])
        print('-dpng',['images\result_' num2str(fig_ind) '.png']);
        fig_ind = fig_ind + 1;
    end
    pause(0.01)
end
ind_a=1;
ind_b=fig_ind-1;
%%------------------
outputVideo = VideoWriter(fullfile('images','mark_decoder_march25_2018.avi'));
outputVideo.FrameRate = 5;
open(outputVideo);
for ii = ind_a:ind_b
   if exist(['images\result_' num2str(ii) '.png'], 'file')
    img = imread(['images\result_' num2str(ii) '.png']);
    writeVideo(outputVideo,img)
   end
end
close(outputVideo);

save('DecoderMean','DecoderMean');