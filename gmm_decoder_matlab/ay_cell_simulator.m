function [Ms,Ps]=ay_cell_simulator(X,CellModel,dt)

Ns  = size(X,1);
M   = length(CellModel);
Ms  = zeros(Ns,M+3);
Ps  = [];
for i=1:Ns
    lambda_cell=zeros(M,1);
    for j=1:M
        lambda_cell(j)=lambda_cal(X(i,:),CellModel{j},dt);
    end
    Ps =[Ps;lambda_cell'];
    Ms(i,2:(M+1))= lambda_cell>rand(M,1);
    
    Ms(i,1)    = i*dt;
    Ms(i,end-1)= X(i,1);
    Ms(i,end)  = X(i,2);
end
MS = Ms;
PS = Ps;
save('MS_PS','MS','PS');
    function lambda=lambda_cal(xm,cell,dt)
        lambda = 0;
        for ii=1:length(cell)
            lambda = lambda + cell(ii).W * mvnpdf(xm,cell(ii).M,cell(ii).C); 
        end
        lambda = lambda*dt;
    end

end