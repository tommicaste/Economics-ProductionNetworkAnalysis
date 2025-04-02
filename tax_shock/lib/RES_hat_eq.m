%%%%%% Basile Grassi, copyright 2020


function [RES,PSigmaY_hat,Equi_hat] = RES_hat_eq_v3(lvec,zi_hat,li_hat,ki_hat,betai_hat,theta,sigma,epsilon,N, Omega,eta,gamma,phi,Delta,psi,w)

%% Load the unknown
    %take the log out
    vec= exp(lvec);
    %vec=real(vec);
    
    


    %Load price
    pi_hat = vec(1:(N));
    
    %Load quantity
    yi_hat = vec((N+1):(2*N));
    
    %GDP
    PSigmaY_hat = vec(2*N+1);
    
    
    
    
%% Parameters
% Omega_ij = share of j in total intermediate inputs cost by i
% eta_i = value added share (of revenue) in sector i
% phi_i = share of final demand i in total output of i
% Delta_ij =  expenditure on j by i as a share of total production of j
% psi_i = share of final demand i in to total final demand
% gamma_i = labor income share in value added in sector i
    
    
%% Some usefull variables
    %power of prices
    pi_hat_1MoinsEpsi = pi_hat.^(1-epsilon);
    pi_hat_1MoinsThetai  = pi_hat.^(1-theta);
    
    
    %Intermediate Bundle Price
    if epsilon == 1
        Pmi_hat = exp( Omega * log(pi_hat)  );
    else
        Pmi_hat = ( Omega * pi_hat_1MoinsEpsi ).^(1/(1-epsilon));
    end
    
    %Value added bundle quantity 
    hi_hat = (li_hat).^gamma .* (ki_hat).^(1-gamma);
    
    %Value added bundle price
    vi_hat = zi_hat.^( (theta-1)/(theta) ) .* (  yi_hat ./ hi_hat  ).^(1/theta).* pi_hat;
    
    %wage
    wi_hat = vi_hat.*hi_hat ./ li_hat;
    %rental rate
    ri_hat = vi_hat.*hi_hat ./ ki_hat;
    
    %final demand
    fi_hat =  betai_hat.*pi_hat.^(-sigma) .* PSigmaY_hat;
    
    %intermediate demand
    xij_hat = ( zi_hat.^(theta-1) .* Pmi_hat.^(epsilon - theta) .* pi_hat.^(theta).* yi_hat   )*( pi_hat.^(-epsilon) )' ;
    %xij_hat = diag( ( zi_hat.^(theta-1) .* Pmi_hat.^(epsilon - theta) .* pi_hat.^(theta).* yi_hat  ) ) * ones(N,N) * diag( ( pi_hat.^(-epsilon) ) ) ;
    
%     xij_hat=zeros(N,N);
%     for i=1:N
%         for j=1:N
%             xij_hat(i,j) = zi_hat(i).^(theta-1) * Pmi_hat(i).^(epsilon - theta) * pi_hat(i).^(theta) * yi_hat(i) * pi_hat(j).^(-epsilon) ;
%         end
%     end
    
    %
    tau_i=zeros(N,1);
    tau_i(55,1)=0.20;


    
%% The residuals equation

%sector price = marginal cost
    if theta == 1
        RES1 = log(pi_hat) - ( -log(zi_hat) + eta.*log(vi_hat) + (1-eta).*log( Pmi_hat ) +log(1+tau_i) );
    else
        RES1 = pi_hat_1MoinsThetai -  (zi_hat).^(theta-1).*( eta.* vi_hat.^(1-theta) + (1-eta).* Pmi_hat.^(1-theta)  ).*(1+tau_i);
    end
    
%sector quantity = market clearing
    %RES2 = ( yi_hat' -( (phi.* fi_hat)' + ones(N,1)'*(Delta.* xij_hat) ) )';
    RES2 =  yi_hat -( (phi.* fi_hat) + (Delta .* xij_hat)'*ones(N,1) ) ;
    
% price index
    
    
    if sigma == 1
        RES3 = 0 - sum( betai_hat.*psi .* log(pi_hat) );
    else
        RES3 = sum( betai_hat.*psi.* pi_hat.^(1-sigma) ) - 1;
    end
    
%The residual
    RES = [ RES1;  RES2; RES3];

%% GDP Ita
%It is the sum of each sector total production minus the intermediate
%bundle used in production, we can do the same with changes in GDP in the
%framework of this model. We need to make a weighted average of the loss,
%where the weight is the relative importanc eof each secotr in Italian
%Output. 
    v_ita= vi_hat(1:54);
    h_ita= hi_hat(1:54);
    vhi_ita= v_ita.*h_ita;
    GDPIta= vhi_ita.*w;
    GDPIta= sum(GDPIta);

    
%% All the variables in a structure
    %Shocks
    Equi_hat.zi=zi_hat;
    
    Equi_hat.li=li_hat;
    
    Equi_hat.ki=ki_hat;


    %Other variables
    Equi_hat.GDPIta = GDPIta;

    Equi_hat.GDP = PSigmaY_hat;
    
    Equi_hat.domar = pi_hat.*yi_hat ./ PSigmaY_hat;
    
    Equi_hat.pyi = pi_hat.*yi_hat;
    
    Equi_hat.pfi =  pi_hat.*fi_hat;
    
    Equi_hat.wli = wi_hat.*li_hat;
    
    Equi_hat.rki = ri_hat.*ki_hat;
    
    Equi_hat.wi=wi_hat;
    
    Equi_hat.ri=ri_hat;
    
    Equi_hat.vi=vi_hat;
    
    Equi_hat.hi=hi_hat;
    
    Equi_hat.vhi=vi_hat.*hi_hat;
    
    Equi_hat.pi=pi_hat;
    
    Equi_hat.yi=yi_hat;
    
    
    
    


end