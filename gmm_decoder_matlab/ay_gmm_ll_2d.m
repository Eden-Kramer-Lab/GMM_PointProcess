function LL=ay_gmm_ll_2d(sk,mk,cell_model,field_model,dt,Xs)
%% spike time (we assume there is only 1 neuron which spikes in this time window)
LL = ones(size(Xs,1),1);
%% build new mixture models, and drop those with tiny chances
ind = find(sk==1);
if ~isempty(ind)   
    % http://compbio.fmph.uniba.sk/vyuka/ml/old/2008/handouts/matrix-cookbook.pdf
    for h=1:length(ind)
        sk_cell_model = cell_model(ind(h));
        m_2    = length(sk_cell_model{1});
        x_gmm  = struct([]);
        for j= 1:m_2
            x_gmm(j).w = dt * sk_cell_model{1}(j).W  * mvnpdf(mk(ind(h),:),sk_cell_model{1}(j).Mm,sk_cell_model{1}(j).Cm);
            x_gmm(j).s = sk_cell_model{1}(j).Cp;
            x_gmm(j).m = sk_cell_model{1}(j).Mp';
        end
        LL = LL .* ay_gmm_posterior_2d(x_gmm,Xs);
    end
end
t_gmm  = struct([]);
for j= 1:length(field_model)
        t_gmm(j).w = dt * field_model(j).W ;
        t_gmm(j).s = field_model(j).Cp;
        t_gmm(j).m = field_model(j).Mp';
end
[temp_l,~,temp_s]=ay_gmm_posterior_2d(t_gmm,Xs);
LL = LL.*exp(-temp_s*temp_l);


