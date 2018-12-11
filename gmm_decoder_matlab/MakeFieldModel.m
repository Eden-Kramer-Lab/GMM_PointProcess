clear all
close all
%% generate model
%load('bn15-2_1-1.mat')
temp_ind = 0;
Xm=[];
Sm=[];
Wm=[]; 
temp_gmm =[];
ind_tot  = 0;
for j=1:9
    j
    if j==1
        load('NewTetroid/bond_md2d1004-8_2-4_tet51.mat');
        MS_1 = decodeme; 
    end
    if j==2
        load('NewTetroid/bond_md2d1004-8_2-4_tet52.mat');
        MS_2 = decodeme; 
    end
    if j==3
        load('NewTetroid/bond_md2d1004-8_2-4_tet53.mat');
        MS_3 = decodeme; 
    end
    if j==4
        load('NewTetroid/bond_md2d1004-8_2-4_tet54.mat');
        MS_4 = decodeme; 
    end
    if j==5
        load('NewTetroid/bond_md2d1004-8_2-4_tet58.mat');
        MS_5 = decodeme; 
    end
    if j==6
        load('NewTetroid/bond_md2d1004-8_2-4_tet61.mat');
        MS_6 = decodeme; 
    end
    if j==7
        load('NewTetroid/bond_md2d1004-8_2-4_tet62.mat');
        MS_7 = decodeme; 
    end
    if j==8
        load('NewTetroid/bond_md2d1004-8_2-4_tet64.mat');
        MS_8 = decodeme; 
    end
    if j==9
        load('NewTetroid/bond_md2d1004-8_2-4_tet68.mat');
        MS_9 = decodeme; 
    end
    
    for i  = 1:length(l0)
                    temp(i).W = l0(i);
                    % mark domain
                    temp(i).Mm = us(i,:);
                    temp(i).Cm = squeeze(Covs(i,:,:));


                    % point process
                    temp(i).Mp = f(i,:);

                    temp(i).Cp = squeeze(sp_cov(i,:,:));

                    Xm = [Xm;temp(i).Mp];
                    Sm = [Sm;temp(i).Cp];
                    Wm = [Wm;temp(i).W];

                    ind_tot = ind_tot + 1;
                    temp_gmm(ind_tot).w = temp(i).W;
                    temp_gmm(ind_tot).s = temp(i).Cp;
                    temp_gmm(ind_tot).m = temp(i).Mp;

    end
    CellModel{j}= temp;
end
for i=1:length(temp_gmm)
    temp_gmm(i).w=temp_gmm(i).w/sum(Wm);
    temp_gmm(i).m =temp_gmm(i).m';
end

temp_gmm   = ay_gmm_merge_alpha_optimized(temp_gmm,temp_gmm,.5);
temp_all   = [];
for i=1:length(temp_gmm)
    temp_all(i).Cp = temp_gmm(i).s;
    temp_all(i).Mp = temp_gmm(i).m';
    temp_all(i).W  = temp_gmm(i).w*sum(Wm);
end
FieldModel = temp_all;
delta_t    = 1e-3;
save('Encoded-9-tetroeds0.5.mat')