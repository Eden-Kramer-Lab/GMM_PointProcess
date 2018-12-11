function ay_cell_model(X0,m_radius,Nc,max_mix)
if nargin==0
    X0  = [0;0];
    m_radius = 100;
    Nc       = 22;
    max_mix  = 4;
end
%% mixture setting
max_scatter = [30  0;
               0  30];
max_weight  = 1000;          
cell_sxy    = [70 0;
               0 70];
       
%% generate circles
angle   = 2 * pi * rand(Nc,1);
radi    = sqrt( rand(Nc,1));
temp_all = [];
temp_ind = 0;
for i  = 1:Nc
    %% center
    xc = X0(1)+ (m_radius * radi(i))*cos(angle(i));
    yc = X0(2)+ (m_radius * radi(i))*sin(angle(i));
    %% generate randi
    nc = randi(max_mix);
    %% now generate nc center
    cs = mvnrnd([xc;yc],max_scatter,nc);
    ws = max_weight * rand(nc,1);
    temp = [];
    for j=1:nc
        temp(j).W = ws(j);
        temp(j).M = cs(j,:);
        temp(j).C = cell_sxy;
        
        temp_ind = temp_ind + 1;
        temp_all(temp_ind).W = ws(j);
        temp_all(temp_ind).M = cs(j,:);
        temp_all(temp_ind).C = cell_sxy;
    end
    CellModel{i}= temp;
end
% Ws = temp_all(1).W;
% Ms = temp_all(1).M;
% for i=2:temp_ind
%     Ws = Ws + temp_all(i).W;
%     Ms = Ms + temp_all(i).M * temp_all(i).W;
% end
% temp_all = [];
% temp_all(1).W = eps;
% temp_all(1).M = [0 0];
% temp_all(1).C = [5000 0;0 5000];

% extra to make sure everything is covered in this circle

temp_ind = temp_ind +1;
temp_all(temp_ind).W= 1;
temp_all(temp_ind).M= X0';
temp_all(temp_ind).C= m_radius^2*eye(length(X0),length(X0));

FieldModel = temp_all;

save('Cell_Field_Model.mat','CellModel','FieldModel');
