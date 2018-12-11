function post_gmm=ay_gmm_approx(pre_gmm,sk,mk,cell_model,field_model,dt,drop_thr)
%% spike time (we assume there is only 1 neuron which spikes in this time window)

%% build new mixture models, and drop those with tiny chances
ind = find(sk==1);
if ~isempty(ind)   
    % http://compbio.fmph.uniba.sk/vyuka/ml/old/2008/handouts/matrix-cookbook.pdf
    org_gmm = pre_gmm;
    for h=1:length(ind)
        sk_cell_model = cell_model(ind(h));
        m_1    = length(org_gmm);
        m_2    = length(sk_cell_model{1});
        gmm_no = 1;
        ws     = zeros(m_1*m_2,1);
        x_gmm  = struct([]);
        
        for i= 1:m_1
            for j= 1:m_2
                m1 = org_gmm(i).m;
                m2 = sk_cell_model{1}(j).Mp';
                
                s1 = org_gmm(i).s;
                s2 = sk_cell_model{1}(j).Cp;
                
                c0 = mvnpdf(m1,m2,s1+s2);
                sc = pinv(pinv(s1)+pinv(s2));
                ms = sc*(pinv(s1)*m1 + pinv(s2)*m2 );
                
                
                x_gmm(gmm_no).w = org_gmm(i).w * sk_cell_model{1}(j).W * c0 * mvnpdf(mk(ind(h),:),sk_cell_model{1}(j).Mm,sk_cell_model{1}(j).Cm);
                x_gmm(gmm_no).s = sc;
                x_gmm(gmm_no).m = ms;
                
                ws(gmm_no) = x_gmm(gmm_no).w;           
                
                gmm_no = gmm_no +1;
             end
        end
        org_gmm = x_gmm;
    end
%%    tmp_gmm =ay_gmm_drop(x_gmm,x_gmm,drop_thr);
    
    tmp_gmm    = struct([]);
    [ws,ws_ind]= sort(ws,'descend');
    ws         = cumsum(ws)/max(realmin,sum(ws));
    tmp_ind    = find(ws>=drop_thr);
    if ~isempty(tmp_ind)
        for i=1:tmp_ind(1)
            tmp_gmm(i).w= x_gmm(ws_ind(i)).w *(dt)^length(ind);
            tmp_gmm(i).s= x_gmm(ws_ind(i)).s ;
            tmp_gmm(i).m= x_gmm(ws_ind(i)).m;
        end
    else
       tmp_gmm = pre_gmm;
    end
else
    tmp_gmm = pre_gmm;
end

%% build post gmm
rs    = 1:-0.1:0;
w_sum = 0;
for i=1:length(tmp_gmm)
    %% linearize around this point 
    s0 = tmp_gmm(i).s;
    m0 = tmp_gmm(i).m;
    %% calculate hessian and gradient
    Hs  = 0*s0;
    Grd = 0*m0;    
    for j=1:length(field_model)
        p_s0   = field_model(j).W * mvnpdf(m0,field_model(j).Mp',field_model(j).Cp);
        inv_s  = pinv(field_model(j).Cp) ;
        Grd    = Grd - p_s0 * inv_s* (m0-field_model(j).Mp');
        Hs     = Hs  + p_s0 * (inv_s*(m0-field_model(j).Mp')*(m0-field_model(j).Mp')'*inv_s-inv_s);
    end
    Grd = Grd * dt;
    Hs  = Hs  * dt;
    %% update variance, then mean
    do_loop = 1;
    ind     = 0;
    while do_loop
        ind = ind +1;
        tmp_gmm(i).s = pinv(pinv(s0)+rs(ind)*Hs);
        do_loop= ~all(eig(tmp_gmm(i).s)>eps);
    end
    tmp_gmm(i).s = 0.5* ( tmp_gmm(i).s+tmp_gmm(i).s');
    % make sure it is positive definite
    xa = eye(size(s0))+0.5*(1-rs(ind))*tmp_gmm(i).s*Hs;
    xb = m0-tmp_gmm(i).s*(Grd-0.5*(1-rs(ind))*Hs*m0);
    tmp_gmm(i).m = pinv(xa)*xb;
    tmp_gmm(i).w = tmp_gmm(i).w * mvnpdf(m0,m0,s0)/max(realmin,mvnpdf(m0,tmp_gmm(i).m,tmp_gmm(i).s));
    w_sum = tmp_gmm(i).w +w_sum;
end
for i=1:length(tmp_gmm)
    tmp_gmm(i).w=tmp_gmm(i).w/w_sum;
end
post_gmm = tmp_gmm;
