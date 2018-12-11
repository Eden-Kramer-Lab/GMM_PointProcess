clear all
close all
%% generate model
%load('bn15-2_1-1.mat')
% temp_ind = 0;
% Xm=[];
% Sm=[];
% Wm=[]; 
% temp_gmm =[];
% ind_tot  = 0;
% for j=1:9
%     j
%     if j==1
%         load('NewTetroid/bond_md2d1004-8_2-4_tet51.mat');
%         MS_1 = decodeme; 
%     end
%     if j==2
%         load('NewTetroid/bond_md2d1004-8_2-4_tet52.mat');
%         MS_2 = decodeme; 
%     end
%     if j==3
%         load('NewTetroid/bond_md2d1004-8_2-4_tet53.mat');
%         MS_3 = decodeme; 
%     end
%     if j==4
%         load('NewTetroid/bond_md2d1004-8_2-4_tet54.mat');
%         MS_4 = decodeme; 
%     end
%     if j==5
%         load('NewTetroid/bond_md2d1004-8_2-4_tet58.mat');
%         MS_5 = decodeme; 
%     end
%     if j==6
%         load('NewTetroid/bond_md2d1004-8_2-4_tet61.mat');
%         MS_6 = decodeme; 
%     end
%     if j==7
%         load('NewTetroid/bond_md2d1004-8_2-4_tet62.mat');
%         MS_7 = decodeme; 
%     end
%     if j==8
%         load('NewTetroid/bond_md2d1004-8_2-4_tet64.mat');
%         MS_8 = decodeme; 
%     end
%     if j==9
%         load('NewTetroid/bond_md2d1004-8_2-4_tet68.mat');
%         MS_9 = decodeme; 
%     end
%     
%     for i  = 1:length(l0)
%                     temp(i).W = l0(i);
%                     % mark domain
%                     temp(i).Mm = us(i,:);
%                     temp(i).Cm = squeeze(Covs(i,:,:));
% 
% 
%                     % point process
%                     temp(i).Mp = f(i,:);
% 
%                     temp(i).Cp = squeeze(sp_cov(i,:,:));
% 
%                     Xm = [Xm;temp(i).Mp];
%                     Sm = [Sm;temp(i).Cp];
%                     Wm = [Wm;temp(i).W];
% 
%                     ind_tot = ind_tot + 1;
%                     temp_gmm(ind_tot).w = temp(i).W;
%                     temp_gmm(ind_tot).s = temp(i).Cp;
%                     temp_gmm(ind_tot).m = temp(i).Mp;
% 
%     end
%     CellModel{j}= temp;
% end
% for i=1:length(temp_gmm)
%     temp_gmm(i).w=temp_gmm(i).w/sum(Wm);
%     temp_gmm(i).m =temp_gmm(i).m';
% end
% 
% temp_gmm   = ay_gmm_merge_alpha_optimized(temp_gmm,temp_gmm,0.00001);
% temp_all   = [];
% for i=1:length(temp_gmm)
%     temp_all(i).Cp = temp_gmm(i).s;
%     temp_all(i).Mp = temp_gmm(i).m';
%     temp_all(i).W  = temp_gmm(i).w*sum(Wm);
% end
% FieldModel = temp_all;
% delta_t    = 1e-3;
 load('Encoded-9-tetroeds0.5.mat')
%% define Ak,Bk, Qk
Ak = [1 0;0 1];
Bk = [0.00;0.00];
Qk = [1e-3 0;
      0    1e-3];
d_step = 0.025;
ds = -0.1:d_step:1.1;
for i=1:length(ds)
    for j=1:length(ds)
        Xs(j+(i-1)*length(ds),:) = [ds(i)  ds(j)];
    end
end
% create cell
%load('MS_PS_1');


%% decoder
start_ind =  1;
end_ind   =  26878;
Ind=0;

% %% Run the exact solution
% post_gmm(1).w = 1;
% post_gmm(1).m = MS_1(start_ind,[1 2])';
% post_gmm(1).s = Qk;
% e_post        = ay_gmm_posterior_2d(post_gmm,Xs);
% 
% hpd      = 0;
% sum_time = 0;
% 
% eDecoderPdf   = zeros(end_ind-start_ind+1,size(Xs,1));
% eDecoderMean  = zeros(end_ind-start_ind+1,2);
% for i=start_ind:end_ind 
%     i
%     Sk  = [MS_1(i,3)  MS_2(i,3) MS_3(i,3) MS_4(i,3)  MS_5(i,3) MS_6(i,3) MS_7(i,3)  MS_8(i,3) MS_9(i,3)];
%     Mk  = [MS_1(i,4:end);
%            MS_2(i,4:end);
%            MS_3(i,4:end);
%            MS_4(i,4:end);
%            MS_5(i,4:end);
%            MS_6(i,4:end);
%            MS_7(i,4:end);
%            MS_8(i,4:end);
%            MS_9(i,4:end)];
%     tic
%     e_one_step = ay_exact_one_step_2d(e_post,Ak,Bk,Qk,Xs);
%     e_ll       = ay_gmm_ll_2d(Sk,Mk,CellModel,FieldModel,delta_t,Xs);
%     
%     e_post     = e_ll.*e_one_step;
%     e_post     = e_post/sum(e_post);
%     time_proc=toc;
%     sum_time = sum_time + time_proc;
%     
%     sample_ind = i-start_ind+1;
%     eDecoderPdf(sample_ind,:)  = e_post;
%     eDecoderMean(sample_ind,:)   = (Xs'*e_post)';
%     
%     [se_port,se_ind] = sort(e_post,'descend');
%     cse_port         = cumsum(se_port);
%     [~,ind] = min(abs(cse_port-0.95));
%     dd = sqrt((Xs(se_ind(1:ind),1)-MS_1(i,1)).^2+(Xs(se_ind(1:ind),2)-MS_1(i,2)).^2);
%     if min(dd)<d_step
%         hpd = hpd +1;
%     end
%     
% %     
% %     if sum(Sk)>0
% %     subplot(1,2,1);
% %     temp = reshape(e_ll,[length(ds)  length(ds)]);
% %     imagesc(ds,ds,temp)
% %     hold on
% %     c=[.8,.8,.8];
% %      plot(MS_1(:,1),MS_1(:,2),'Color',[.5,.5,.5],'Linewidth',6);
% %     plot(MS_1(start_ind:end_ind,1),MS_1(start_ind:end_ind,2),'Color',c,'Linewidth',8);
% %     
% %     c=[0,0,0.8];
% %     plot(eDecoderMean(1:sample_ind-1,1),eDecoderMean(1:sample_ind-1,2),'g','Linewidth',2);
% %     
% %     h3=scatter(eDecoderMean(sample_ind,1),eDecoderMean(sample_ind,2),'d');
% %     h3.LineWidth = 6;
% %     h3.MarkerEdgeColor = 'r';
% %     h3.MarkerFaceColor = [0.5 0.5 0.0];
% %     
% %     h4=scatter(MS_1(i,1),MS_1(i,2),'d');
% %     h4.LineWidth = 6;
% %     h4.MarkerEdgeColor = 'y';
% %     hold off
% %     set(gca,'YDir','reverse');
% %     set(gca,'XDir','reverse');
% %     title('GMM-2D');
% %     
% %     subplot(1,2,2);
% %     temp = reshape(e_post,[length(ds)  length(ds)]);
% %     imagesc(ds,ds,temp)
% %     hold on
% %     c=[.8,.8,.8];
% %      plot(MS_1(:,1),MS_1(:,2),'Color',[.5,.5,.5],'Linewidth',6);
% %     plot(MS_1(start_ind:end_ind,1),MS_1(start_ind:end_ind,2),'Color',c,'Linewidth',8);
% %     
% %     c=[0,0,0.8];
% %     plot(eDecoderMean(1:sample_ind-1,1),eDecoderMean(1:sample_ind-1,2),'g','Linewidth',2);
% %     
% %     h3=scatter(eDecoderMean(sample_ind,1),eDecoderMean(sample_ind,2),'d');
% %     h3.LineWidth = 6;
% %     h3.MarkerEdgeColor = 'r';
% %     h3.MarkerFaceColor = [0.5 0.5 0.0];
% %     
% %     h4=scatter(MS_1(i,1),MS_1(i,2),'d');
% %     h4.LineWidth = 6;
% %     h4.MarkerEdgeColor = 'y';
% %     hold off
% %     set(gca,'YDir','reverse');
% %     set(gca,'XDir','reverse');
% %     title('Filter ')
% %     saveas(gcf,['Result/Imgs/Img',num2str(Ind),'.png'])
% %     Ind=Ind+1;
% %     end
% end
% err = sqrt(mean((MS_1(start_ind:end_ind,1)-eDecoderMean(:,1)).^2 + (MS_1(start_ind:end_ind,2)-eDecoderMean(:,2)).^2));
% hpd = hpd/(end_ind-start_ind+1);
% avg_time = sum_time/(end_ind-start_ind+1);
% RESULT(1,1)=err;
% RESULT(1,2)=hpd*100;
% RESULT(1,3)=avg_time*1000;
% 
% 
% save('Result/ExactResult.mat','eDecoderPdf','eDecoderMean','RESULT');
% 
% 

%% GMM 1
post_gmm = [];
post_gmm(1).w = 1;
post_gmm(1).m = MS_1(start_ind,[1 2])';
post_gmm(1).s = Qk;

merge_c  = 0.05;
drop_c   = 1-0.1;
aDecoderPdf   = zeros(end_ind-start_ind+1,length(Xs));
aDecoderMean  = zeros(end_ind-start_ind+1,2);
sum_mix       = zeros(end_ind-start_ind+1,3);
hpd       = 0;
sum_time  = 0;
Ind       = 0;
for i=start_ind:end_ind  
    i
    if i==1264
        stop=1;
    end
   %% Sk
    Sk  = [MS_1(i,3) MS_2(i,3) MS_3(i,3) MS_4(i,3) MS_5(i,3) MS_6(i,3) MS_7(i,3) MS_8(i,3) MS_9(i,3)];
    Mk  = [MS_1(i,4:end);
           MS_2(i,4:end);
           MS_3(i,4:end);
           MS_4(i,4:end);
           MS_5(i,4:end);
           MS_6(i,4:end);
           MS_7(i,4:end);
           MS_8(i,4:end);
           MS_9(i,4:end)];
    tic
    %% prior
    pre_gmm    = post_gmm;
    %% one-step 
    one_step   = ay_gmm_one_step(pre_gmm,Ak,Bk,Qk);
    %% filter
    if sum(Sk)>1
        ind        = find(Sk==1);
        t_one_step = one_step;
        for ss=1:sum(Sk)
            tSk    = 0*Sk;
            tSk(ind(ss))=1;
            t_one_step   = ay_gmm_approx_fast(t_one_step,tSk,Mk,CellModel,FieldModel,delta_t/sum(Sk),drop_c,2);
            if length(t_one_step) > 1
                    [t_one_step,xxx_gmm] = ay_gmm_merge_alpha_optimized(t_one_step,t_one_step,merge_c);
            end
        end
        post_gmm = t_one_step;
    else
        % here, we don't check posteior covariance matrix to be PSD (this is correct,as we have changed the ground truth using merging)
        post_gmm   = ay_gmm_approx_fast(one_step,Sk,Mk,CellModel,FieldModel,delta_t,drop_c,2);
        %% merge up to 1
        if length(post_gmm) > 1
            [post_gmm,xxx_gmm] = ay_gmm_merge_alpha_optimized(post_gmm,post_gmm,merge_c);
        end
    end
    
    %% we can use this if we want to check each covariance matrix to be PSD
    % post_gmm   = ay_gmm_approx(one_step,Sk,Mk,CellModel,FieldModel,delta_t,drop_c);
    px_gmm     = post_gmm;
    %%% MY
    [e_post0,hjh]=ay_gmm_posterior_2d(post_gmm,Xs);
 
    %% merge up to 1
    if length(post_gmm) > 1
        [post_gmm,xxx_gmm] = ay_gmm_merge_alpha_optimized(post_gmm,post_gmm,merge_c);
    end
    time_proc= toc;
    sum_time = sum_time + time_proc;
    
    
    sample_ind = i-start_ind+1;
    [aDecoderPdf(sample_ind,:),aDecoderMean(sample_ind,:)]=ay_gmm_posterior_2d(post_gmm,Xs);
    e_post = aDecoderPdf(sample_ind,:);
    sum_mix(sample_ind,1)  =  length(post_gmm);
    sum_mix(sample_ind,2)  =  sum(Sk);
    sum_mix(sample_ind,3)  =  length(px_gmm);
    
    [se_port,se_ind] = sort(e_post,'descend');
    cse_port         = cumsum(se_port);
    [~,ind] = min(abs(cse_port-0.95));
    dd = sqrt((Xs(se_ind(1:ind),1)-MS_1(i,1)).^2+(Xs(se_ind(1:ind),2)-MS_1(i,2)).^2);
    if min(dd)<d_step
        hpd = hpd +1;
    end
    
     
    if sum(Sk)>100
    subplot(1,2,1);
    temp = reshape(e_post0,[length(ds)  length(ds)]);
    imagesc(ds,ds,temp)
    hold on
    c=[.8,.8,.8];
     plot(MS_1(:,1),MS_1(:,2),'Color',[.5,.5,.5],'Linewidth',6);
    plot(MS_1(start_ind:end_ind,1),MS_1(start_ind:end_ind,2),'Color',c,'Linewidth',8);
    
    c=[0,0,0.8];
    plot(aDecoderMean(1:sample_ind-1,1),aDecoderMean(1:sample_ind-1,2),'g','Linewidth',2);
    
    h3=scatter(aDecoderMean(sample_ind,1),aDecoderMean(sample_ind,2),'d');
    h3.LineWidth = 6;
    h3.MarkerEdgeColor = 'r';
    h3.MarkerFaceColor = [0.5 0.5 0.0];
    
    h4=scatter(MS_1(i,1),MS_1(i,2),'d');
    h4.LineWidth = 6;
    h4.MarkerEdgeColor = 'y';
    hold off
    set(gca,'YDir','reverse');
    set(gca,'XDir','reverse');
    title('After drop');
    
    subplot(1,2,2);
    temp = reshape(e_post,[length(ds)  length(ds)]);
    imagesc(ds,ds,temp)
    hold on
    c=[.8,.8,.8];
     plot(MS_1(:,1),MS_1(:,2),'Color',[.5,.5,.5],'Linewidth',6);
    plot(MS_1(start_ind:end_ind,1),MS_1(start_ind:end_ind,2),'Color',c,'Linewidth',8);
    
    c=[0,0,0.8];
    plot(aDecoderMean(1:sample_ind-1,1),aDecoderMean(1:sample_ind-1,2),'g','Linewidth',2);
    
    h3=scatter(aDecoderMean(sample_ind,1),aDecoderMean(sample_ind,2),'d');
    h3.LineWidth = 6;
    h3.MarkerEdgeColor = 'r';
    h3.MarkerFaceColor = [0.5 0.5 0.0];
    
    h4=scatter(MS_1(i,1),MS_1(i,2),'d');
    h4.LineWidth = 6;
    h4.MarkerEdgeColor = 'y';
    hold off
    set(gca,'YDir','reverse');
    set(gca,'XDir','reverse');
    title('After Merge ')
    saveas(gcf,['Result/Img/Img',num2str(Ind),'.png'])
    Ind=Ind+1;
    end

        
end

err = sqrt(mean((MS_1(start_ind:end_ind,1)-aDecoderMean(:,1)).^2 + (MS_1(start_ind:end_ind,2)-aDecoderMean(:,2)).^2));
hpd = hpd/(end_ind-start_ind+1);
avg_time = sum_time/(end_ind-start_ind+1);
RESULTg(2,1)=err;
RESULTg(2,2)=hpd*100;
RESULTg(2,3)=avg_time*1000;
ind = find(sum_mix(:,2)==1);
RESULTg(2,4)=mean(sum_mix(:,1));
RESULTg(2,6)=mean(sum_mix(ind,3));
RESULTg(2,5)=mean(sum_mix(ind,1));
RESULTg(2,7)=max(sum_mix(ind,1));
RESULTg(2,8)=max(sum_mix(ind,3));

% save('Result/GMMResult0.1.mat','aDecoderPdf','aDecoderMean','RESULTg');

%% GMM 2
post_gmm = [];
post_gmm(1).w = 1;
post_gmm(1).m = MS_1(start_ind,[1 2])';
post_gmm(1).s = Qk;

merge_c  = 0.15;
drop_c   = 1-0.15;
bDecoderPdf   = zeros(end_ind-start_ind+1,length(Xs));
bDecoderMean  = zeros(end_ind-start_ind+1,2);
sum_mix       = zeros(end_ind-start_ind+1,3);
hpd       = 0;
sum_time  = 0;
Ind       = 0;
for i=start_ind:end_ind  
    i
    if i==1264
        stop=1;
    end
   %% Sk
    Sk  = [MS_1(i,3) MS_2(i,3) MS_3(i,3) MS_4(i,3) MS_5(i,3) MS_6(i,3) MS_7(i,3) MS_8(i,3) MS_9(i,3)];
    Mk  = [MS_1(i,4:end);
           MS_2(i,4:end);
           MS_3(i,4:end);
           MS_4(i,4:end);
           MS_5(i,4:end);
           MS_6(i,4:end);
           MS_7(i,4:end);
           MS_8(i,4:end);
           MS_9(i,4:end)];
    tic
    %% prior
    pre_gmm    = post_gmm;
    %% one-step 
    one_step   = ay_gmm_one_step(pre_gmm,Ak,Bk,Qk);
    %% filter
    if sum(Sk)>1
        ind        = find(Sk==1);
        t_one_step = one_step;
        for ss=1:sum(Sk)
            tSk    = 0*Sk;
            tSk(ind(ss))=1;
            t_one_step   = ay_gmm_approx_fast(t_one_step,tSk,Mk,CellModel,FieldModel,delta_t/sum(Sk),drop_c,2);
            if length(t_one_step) > 1
                    [t_one_step,xxx_gmm] = ay_gmm_merge_alpha_optimized(t_one_step,t_one_step,merge_c);
            end
        end
        post_gmm = t_one_step;
    else
        % here, we don't check posteior covariance matrix to be PSD (this is correct,as we have changed the ground truth using merging)
        post_gmm   = ay_gmm_approx_fast(one_step,Sk,Mk,CellModel,FieldModel,delta_t,drop_c,2);
        %% merge up to 1
        if length(post_gmm) > 1
            [post_gmm,xxx_gmm] = ay_gmm_merge_alpha_optimized(post_gmm,post_gmm,merge_c);
        end
    end
    
    %% we can use this if we want to check each covariance matrix to be PSD
    % post_gmm   = ay_gmm_approx(one_step,Sk,Mk,CellModel,FieldModel,delta_t,drop_c);
    px_gmm     = post_gmm;
    %%% MY
    [e_post0,hjh]=ay_gmm_posterior_2d(post_gmm,Xs);
 
    %% merge up to 1
    if length(post_gmm) > 1
        [post_gmm,xxx_gmm] = ay_gmm_merge_alpha_optimized(post_gmm,post_gmm,merge_c);
    end
    time_proc= toc;
    sum_time = sum_time + time_proc;
    
    
    sample_ind = i-start_ind+1;
    [bDecoderPdf(sample_ind,:),bDecoderMean(sample_ind,:)]=ay_gmm_posterior_2d(post_gmm,Xs);
    e_post = bDecoderPdf(sample_ind,:);
    sum_mix(sample_ind,1)  =  length(post_gmm);
    sum_mix(sample_ind,2)  =  sum(Sk);
    sum_mix(sample_ind,3)  =  length(px_gmm);
    
    [se_port,se_ind] = sort(e_post,'descend');
    cse_port         = cumsum(se_port);
    [~,ind] = min(abs(cse_port-0.95));
    dd = sqrt((Xs(se_ind(1:ind),1)-MS_1(i,1)).^2+(Xs(se_ind(1:ind),2)-MS_1(i,2)).^2);
    if min(dd)<d_step
        hpd = hpd +1;
    end
    
     
    if sum(Sk)>100
    subplot(1,2,1);
    temp = reshape(e_post0,[length(ds)  length(ds)]);
    imagesc(ds,ds,temp)
    hold on
    c=[.8,.8,.8];
     plot(MS_1(:,1),MS_1(:,2),'Color',[.5,.5,.5],'Linewidth',6);
    plot(MS_1(start_ind:end_ind,1),MS_1(start_ind:end_ind,2),'Color',c,'Linewidth',8);
    
    c=[0,0,0.8];
    plot(bDecoderMean(1:sample_ind-1,1),bDecoderMean(1:sample_ind-1,2),'g','Linewidth',2);
    
    h3=scatter(bDecoderMean(sample_ind,1),bDecoderMean(sample_ind,2),'d');
    h3.LineWidth = 6;
    h3.MarkerEdgeColor = 'r';
    h3.MarkerFaceColor = [0.5 0.5 0.0];
    
    h4=scatter(MS_1(i,1),MS_1(i,2),'d');
    h4.LineWidth = 6;
    h4.MarkerEdgeColor = 'y';
    hold off
    set(gca,'YDir','reverse');
    set(gca,'XDir','reverse');
    title('After drop');
    
    subplot(1,2,2);
    temp = reshape(e_post,[length(ds)  length(ds)]);
    imagesc(ds,ds,temp)
    hold on
    c=[.8,.8,.8];
     plot(MS_1(:,1),MS_1(:,2),'Color',[.5,.5,.5],'Linewidth',6);
    plot(MS_1(start_ind:end_ind,1),MS_1(start_ind:end_ind,2),'Color',c,'Linewidth',8);
    
    c=[0,0,0.8];
    plot(bDecoderMean(1:sample_ind-1,1),bDecoderMean(1:sample_ind-1,2),'g','Linewidth',2);
    
    h3=scatter(bDecoderMean(sample_ind,1),bDecoderMean(sample_ind,2),'d');
    h3.LineWidth = 6;
    h3.MarkerEdgeColor = 'r';
    h3.MarkerFaceColor = [0.5 0.5 0.0];
    
    h4=scatter(MS_1(i,1),MS_1(i,2),'d');
    h4.LineWidth = 6;
    h4.MarkerEdgeColor = 'y';
    hold off
    set(gca,'YDir','reverse');
    set(gca,'XDir','reverse');
    title('After Merge ')
    saveas(gcf,['Result/Img/Img',num2str(Ind),'.png'])
    Ind=Ind+1;
    end

        
end

err = sqrt(mean((MS_1(start_ind:end_ind,1)-bDecoderMean(:,1)).^2 + (MS_1(start_ind:end_ind,2)-bDecoderMean(:,2)).^2));
hpd = hpd/(end_ind-start_ind+1);
avg_time = sum_time/(end_ind-start_ind+1);
RESULTg(3,1)=err;
RESULTg(3,2)=hpd*100;
RESULTg(3,3)=avg_time*1000;
ind = find(sum_mix(:,2)==1);
RESULTg(3,4)=mean(sum_mix(:,1));
RESULTg(3,6)=mean(sum_mix(ind,3));
RESULTg(3,5)=mean(sum_mix(ind,1));
RESULTg(3,7)=max(sum_mix(ind,1));
RESULTg(3,8)=max(sum_mix(ind,3));


%% Plot Result model trajectory result
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

sample_ind = 24000;
% trajectory
% for i=start_ind:end_ind
subplot(1,3,1)
 temp = reshape(eDecoderPdf(sample_ind,:)',[length(ds)  length(ds)]);
 imagesc(ds,ds,temp)
 hold on
% end
plot(MS_1(:,1),MS_1(:,2),'color',[.4 .4 .4],'LineWidth',6);
scatter(eDecoderMean(sample_ind,1),eDecoderMean(sample_ind,2),'db','LineWidth',4);
scatter(MS_1(sample_ind,1),MS_1(sample_ind,2),'or','LineWidth',4);

hold off
xlim([.1 .9])
ylim([-.1 1.1])
colormap parula
ylabel('X')
xlabel('Time')
title('Exact Solution')

subplot(1,3,2)
 temp = reshape(aDecoderPdf(sample_ind,:)',[length(ds)  length(ds)]);
 imagesc(ds,ds,temp)
 hold on
% end
plot(MS_1(:,1),MS_1(:,2),'color',[.4 .4 .4],'LineWidth',6);
scatter(aDecoderMean(sample_ind,1),aDecoderMean(sample_ind,2),'db','LineWidth',4);
scatter(MS_1(sample_ind,1),MS_1(sample_ind,2),'or','LineWidth',4);

hold off
xlim([.1 .9])
ylim([-.1 1.1])
colormap parula
ylabel('X')
xlabel('Time')
title('GMM_1')

subplot(1,3,3)
 temp = reshape(bDecoderPdf(sample_ind,:)',[length(ds)  length(ds)]);
 imagesc(ds,ds,temp)
 hold on
% end
plot(MS_1(:,1),MS_1(:,2),'color',[.4 .4 .4],'LineWidth',6);
scatter(bDecoderMean(sample_ind,1),bDecoderMean(sample_ind,2),'db','LineWidth',4);
scatter(MS_1(sample_ind,1),MS_1(sample_ind,2),'or','LineWidth',4);

hold off
xlim([.1 .9])
ylim([-.1 1.1])
colormap parula
ylabel('X')
xlabel('Time')
title('GMM_2')


% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 4])
% % print('-dpng',['images\result_exact_map'  '.png']);
