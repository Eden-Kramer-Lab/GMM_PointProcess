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
d_step = 0.01;
Xs = [-12:d_step:12];
% create cell
load('MS_PS_1');

%% decoder
start_ind =  96380;
end_ind   =  123917;

%% Run the exact solution
post_gmm(1).w = 1;
post_gmm(1).m = MS(start_ind,1);
post_gmm(1).s = 10*Qk;

e_post        = ay_gmm_posterior_1d(post_gmm,Xs);

hpd      = 0;
sum_time = 0;

eDecoderPdf   = zeros(end_ind-start_ind+1,length(Xs));
eDecoderMean  = zeros(end_ind-start_ind+1,1);
for i=start_ind:end_ind  
    i
    Sk  = MS(i,2);
    Mk  = MS(i,3:end);
    tic
    e_one_step = ay_exact_one_step(e_post,Ak,0*Bk,Qk,Xs);
    e_ll       = ay_gmm_ll(Sk,Mk,CellModel,FieldModel,delta_t,Xs);
    e_post     = e_ll.*e_one_step;
    e_post     = e_post/sum(e_post);
    time_proc=toc;
    sum_time = sum_time + time_proc;
    
    sample_ind = i-start_ind+1;
    eDecoderPdf(sample_ind,:)  = e_post;
    eDecoderMean(sample_ind) = e_post*Xs';
    
    [se_port,se_ind] = sort(e_post,'descend');
    cse_port         = cumsum(se_port);
    [~,ind] = min(abs(cse_port-0.95));
    if min(abs(Xs(se_ind(1:ind))-MS(i,1)))<d_step
        hpd = hpd +1;
    end
end
err = sqrt(mean((MS(start_ind:end_ind,1)-eDecoderMean).^2));
hpd = hpd/(end_ind-start_ind+1);
avg_time = sum_time/(end_ind-start_ind+1);
RESULT(1,1)=err;
RESULT(1,2)=hpd*100;
RESULT(1,3)=avg_time*1000;


%% GMM merge parameters
post_gmm(1).w = 1;
post_gmm(1).m = MS(start_ind,1);
post_gmm(1).s = 10*Qk;

merge_c = 0.001;
drop_c  = 0.9999;
aDecoderPdf   = zeros(end_ind-start_ind+1,length(Xs));
aDecoderMean  = zeros(end_ind-start_ind+1,1);
sum_mix       = zeros(end_ind-start_ind+1,3);
hpd      = 0;
sum_time = 0;

for i=start_ind:end_ind  
    i
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
    time_proc= toc;
    sum_time = sum_time + time_proc;
    
    
    sample_ind = i-start_ind+1;
    [aDecoderPdf(sample_ind,:),aDecoderMean(sample_ind,:)]=ay_gmm_posterior_1d(post_gmm,Xs);
    e_post = aDecoderPdf(sample_ind,:);
    sum_mix(sample_ind,1)  =  length(post_gmm);
    sum_mix(sample_ind,2)  =  Sk;
    sum_mix(sample_ind,3)  =  length(px_gmm);
    
    [se_port,se_ind] = sort(e_post,'descend');
    cse_port         = cumsum(se_port);
    [~,ind] = min(abs(cse_port-0.95));
    if min(abs(Xs(se_ind(1:ind))-MS(i,1)))<d_step
        hpd = hpd +1;
    end
    
end
err = sqrt(mean((MS(start_ind:end_ind,1)-aDecoderMean).^2));
hpd = hpd/(end_ind-start_ind+1);
avg_time = sum_time/(end_ind-start_ind+1);
RESULT(2,1)=err;
RESULT(2,2)=hpd*100;
RESULT(2,3)=avg_time*1000;
ind = find(sum_mix(:,2)==1);
RESULT(2,4)=mean(sum_mix(:,1));
RESULT(2,6)=mean(sum_mix(ind,3));
RESULT(2,5)=mean(sum_mix(ind,1));
RESULT(2,7)=max(sum_mix(ind,1));
RESULT(2,8)=max(sum_mix(ind,3));

%% GMM merge parameters
post_gmm(1).w = 1;
post_gmm(1).m = MS(start_ind,1);
post_gmm(1).s = 10*Qk;

merge_c = 0.05;
drop_c  = 0.9999;
bDecoderPdf   = zeros(end_ind-start_ind+1,length(Xs));
bDecoderMean  = zeros(end_ind-start_ind+1,1);
sum_mix       = zeros(end_ind-start_ind+1,3);
hpd      = 0;
sum_time = 0;

for i=start_ind:end_ind  
    i
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
    time_proc= toc;
    sum_time = sum_time + time_proc;
    
    
    sample_ind = i-start_ind+1;
    [bDecoderPdf(sample_ind,:),bDecoderMean(sample_ind,:)]=ay_gmm_posterior_1d(post_gmm,Xs);
    e_post = bDecoderPdf(sample_ind,:);
    sum_mix(sample_ind,1)  =  length(post_gmm);
    sum_mix(sample_ind,2)  =  Sk;
    sum_mix(sample_ind,3)  =  length(px_gmm);
    
    [se_port,se_ind] = sort(e_post,'descend');
    cse_port         = cumsum(se_port);
    [~,ind] = min(abs(cse_port-0.95));
    if min(abs(Xs(se_ind(1:ind))-MS(i,1)))<d_step
        hpd = hpd +1;
    end
    
end
err = sqrt(mean((MS(start_ind:end_ind,1)-bDecoderMean).^2));
hpd = hpd/(end_ind-start_ind+1);
avg_time = sum_time/(end_ind-start_ind+1);
RESULT(3,1)=err;
RESULT(3,2)=hpd*100;
RESULT(3,3)=avg_time*1000;
ind = find(sum_mix(:,2)==1);
RESULT(3,4)=mean(sum_mix(:,1));
RESULT(3,6)=mean(sum_mix(ind,3));
RESULT(3,5)=mean(sum_mix(ind,1));
RESULT(3,7)=max(sum_mix(ind,1));
RESULT(3,8)=max(sum_mix(ind,3));

%% GMM merge parameters
post_gmm(1).w = 1;
post_gmm(1).m = MS(start_ind,1);
post_gmm(1).s = 10*Qk;

merge_c = 0.01;
drop_c  = 0.999;
cDecoderPdf   = zeros(end_ind-start_ind+1,length(Xs));
cDecoderMean  = zeros(end_ind-start_ind+1,1);
sum_mix       = zeros(end_ind-start_ind+1,3);
hpd      = 0;
sum_time = 0;

for i=start_ind:end_ind  
    i
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
    time_proc= toc;
    sum_time = sum_time + time_proc;
    
    
    sample_ind = i-start_ind+1;
    [cDecoderPdf(sample_ind,:),cDecoderMean(sample_ind,:)]=ay_gmm_posterior_1d(post_gmm,Xs);
    e_post = cDecoderPdf(sample_ind,:);
    sum_mix(sample_ind,1)  =  length(post_gmm);
    sum_mix(sample_ind,2)  =  Sk;
    sum_mix(sample_ind,3)  =  length(px_gmm);
    
    [se_port,se_ind] = sort(e_post,'descend');
    cse_port         = cumsum(se_port);
    [~,ind] = min(abs(cse_port-0.95));
    if min(abs(Xs(se_ind(1:ind))-MS(i,1)))<d_step
        hpd = hpd +1;
    end
    
end
err = sqrt(mean((MS(start_ind:end_ind,1)-cDecoderMean).^2));
hpd = hpd/(end_ind-start_ind+1);
avg_time = sum_time/(end_ind-start_ind+1);
RESULT(4,1)=err;
RESULT(4,2)=hpd*100;
RESULT(4,3)=avg_time*1000;
ind = find(sum_mix(:,2)==1);
RESULT(4,4)=mean(sum_mix(:,1));
RESULT(4,6)=mean(sum_mix(ind,3));
RESULT(4,5)=mean(sum_mix(ind,1));
RESULT(4,7)=max(sum_mix(ind,1));
RESULT(4,8)=max(sum_mix(ind,3));
%% GMM merge parameters
post_gmm(1).w = 1;
post_gmm(1).m = MS(start_ind,1);
post_gmm(1).s = 10*Qk;

merge_c = 0.05;
drop_c  = 0.995;
dDecoderPdf   = zeros(end_ind-start_ind+1,length(Xs));
dDecoderMean  = zeros(end_ind-start_ind+1,1);
sum_mix       = zeros(end_ind-start_ind+1,3);
hpd      = 0;
sum_time = 0;

for i=start_ind:end_ind  
    i
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
    time_proc= toc;
    sum_time = sum_time + time_proc;
    
    
    sample_ind = i-start_ind+1;
    [dDecoderPdf(sample_ind,:),dDecoderMean(sample_ind,:)]=ay_gmm_posterior_1d(post_gmm,Xs);
    e_post = dDecoderPdf(sample_ind,:);
    sum_mix(sample_ind,1)  =  length(post_gmm);
    sum_mix(sample_ind,2)  =  Sk;
    sum_mix(sample_ind,3)  =  length(px_gmm);
    
    [se_port,se_ind] = sort(e_post,'descend');
    cse_port         = cumsum(se_port);
    [~,ind] = min(abs(cse_port-0.95));
    if min(abs(Xs(se_ind(1:ind))-MS(i,1)))<d_step
        hpd = hpd +1;
    end
    
end
err = sqrt(mean((MS(start_ind:end_ind,1)-dDecoderMean).^2));
hpd = hpd/(end_ind-start_ind+1);
avg_time = sum_time/(end_ind-start_ind+1);
RESULT(5,1)=err;
RESULT(5,2)=hpd*100;
RESULT(5,3)=avg_time*1000;
ind = find(sum_mix(:,2)==1);
RESULT(5,4)=mean(sum_mix(:,1));
RESULT(5,6)=mean(sum_mix(ind,3));
RESULT(5,5)=mean(sum_mix(ind,1));
RESULT(5,7)=max(sum_mix(ind,1));
RESULT(5,8)=max(sum_mix(ind,3));


%% GMM merge parameters
post_gmm(1).w = 1;
post_gmm(1).m = MS(start_ind,1);
post_gmm(1).s = 10*Qk;

merge_c = 0.5;
drop_c  = 0.9999;
sDecoderPdf   = zeros(end_ind-start_ind+1,length(Xs));
sDecoderMean  = zeros(end_ind-start_ind+1,1);
sum_mix       = zeros(end_ind-start_ind+1,3);
hpd      = 0;
sum_time = 0;

for i=start_ind:end_ind  
    i
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
    if Sk ==1
        [post_gmm,xxx_gmm] = ay_gmm_merge_alpha_optimized(post_gmm,post_gmm,merge_c);
    end
    time_proc= toc;
    sum_time = sum_time + time_proc;
    
    
    sample_ind = i-start_ind+1;
    [sDecoderPdf(sample_ind,:),sDecoderMean(sample_ind,:)]=ay_gmm_posterior_1d(post_gmm,Xs);
    e_post = sDecoderPdf(sample_ind,:);
    sum_mix(sample_ind,1)  =  length(post_gmm);
    sum_mix(sample_ind,2)  =  Sk;
    sum_mix(sample_ind,3)  =  length(px_gmm);
    
    [se_port,se_ind] = sort(e_post,'descend');
    cse_port         = cumsum(se_port);
    [~,ind] = min(abs(cse_port-0.95));
    if min(abs(Xs(se_ind(1:ind))-MS(i,1)))<d_step
        hpd = hpd +1;
    end
    
end
err = sqrt(mean((MS(start_ind:end_ind,1)-sDecoderMean).^2));
hpd = hpd/(end_ind-start_ind+1);
avg_time = sum_time/(end_ind-start_ind+1);
RESULT(6,1)=err;
RESULT(6,2)=hpd*100;
RESULT(6,3)=avg_time*1000;
ind = find(sum_mix(:,2)==1);
RESULT(6,4)=mean(sum_mix(:,1));
RESULT(6,6)=mean(sum_mix(ind,3));
RESULT(6,5)=mean(sum_mix(ind,1));
RESULT(6,7)=max(sum_mix(ind,1));
RESULT(6,8)=max(sum_mix(ind,3));

% Set up figure properties:
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.


%% figure 1
sample_ind = end_ind-start_ind+1;
% trajectory
imagesc(delta_t *(1:sample_ind),Xs,sqrt(eDecoderPdf(1:sample_ind,:)'),[1e-3 1e-1]);
hold on
plot(delta_t *(1:sample_ind),eDecoderMean(1:sample_ind,:),'r','LineWidth',2);
plot(delta_t *(1:sample_ind),MS(start_ind:end_ind,1),'w','LineWidth',2);hold off
ylim([-7 7])
colormap parula
ylabel('X')
xlabel('Time')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 4])
print('-dpng',['images\result_exact_map'  '.png']);

%% figure 2
% trajectory
imagesc(delta_t *(1:sample_ind),Xs,sqrt(aDecoderPdf(1:sample_ind,:)'),[1e-3 1e-1]);
hold on
plot(delta_t *(1:sample_ind),aDecoderMean(1:sample_ind,:),'r','LineWidth',2);
plot(delta_t *(1:sample_ind),MS(start_ind:end_ind,1),'w','LineWidth',2);hold off
ylim([-7 7])
colormap parula
ylabel('X')
xlabel('Time')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 4])
print('-dpng',['images\result_99_1_map'  '.png']);


%% figure 3
% trajectory
imagesc(delta_t *(1:sample_ind),Xs,sqrt(bDecoderPdf(1:sample_ind,:)'),[1e-3 1e-1]);
hold on
plot(delta_t *(1:sample_ind),bDecoderMean(1:sample_ind,:),'r','LineWidth',2);
plot(delta_t *(1:sample_ind),MS(start_ind:end_ind,1),'w','LineWidth',2);hold off
ylim([-7 7])
colormap parula
ylabel('X')
xlabel('Time')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 4])
print('-dpng',['images\result_99_5_map'  '.png']);

%% figure 3
% trajectory
imagesc(delta_t *(1:sample_ind),Xs,sqrt(cDecoderPdf(1:sample_ind,:)'),[1e-3 1e-1]);
hold on
plot(delta_t *(1:sample_ind),cDecoderMean(1:sample_ind,:),'r','LineWidth',2);
plot(delta_t *(1:sample_ind),MS(start_ind:end_ind,1),'w','LineWidth',2);hold off
ylim([-7 7])
colormap parula
ylabel('X')
xlabel('Time')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 4])
print('-dpng',['images\result_95_1_map'  '.png']);

%% figure 4
% trajectory
imagesc(delta_t *(1:sample_ind),Xs,sqrt(dDecoderPdf(1:sample_ind,:)'),[1e-3 1e-1]);
hold on
plot(delta_t *(1:sample_ind),dDecoderMean(1:sample_ind,:),'r','LineWidth',2);
plot(delta_t *(1:sample_ind),MS(start_ind:end_ind,1),'w','LineWidth',2);hold off
ylim([-7 7])
colormap parula
ylabel('X')
xlabel('Time')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 4])
print('-dpng',['images\result_95_5_map'  '.png']);

%% figure 5
plot(delta_t *(1:sample_ind),MS(start_ind:end_ind,1),'k','LineWidth',2);
set(gca,'Ydir','reverse')
axis tight
ylim([-7 7])
colormap parula
ylabel('X')
xlabel('Time')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 4])
print('-dpng',['images\result_trajectory_map'  '.png']);

save('DecoderResult','RESULT');