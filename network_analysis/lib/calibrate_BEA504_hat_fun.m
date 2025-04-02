%% Calibrate the model in "Sectoral Effect of Social Distancing"

%% Copyright: Jean-Noel Barrot, Basile Grassi and Julien Sauvagnat
%%AEA Paper and Procedings May 2021
%%Solving the model in BGS (AEA P&P) and performing exercice for various labor supply shocks
%%This version: the AEA P&P
%%Corresponding author: Basile Grassi



function calibrate_BEA504_hat_fun(theta,sigma,epsilon)

%% I %%%%%%%%%%%%%%%%%%%%%% Import DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('%% I %%%%%%%%%%%%%%%%%%%%%% Import DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    tic;
    %final consumption
        C = xlsread('data_deep/WIOT2014.xlsx','2014bis','BM3:BM56');
        G = xlsread('data_deep/WIOT2014.xlsx','2014bis','BN3:BN56');
        I = xlsread('data_deep/WIOT2014.xlsx','2014bis','BO3:BO56');
        pfi = C+G+I
        pfi(isnan(pfi))= 0  
        pfi(pfi<0) = 0


    %IO
        OmegaRaw = xlsread('data_deep/WIOT2014.xlsx','2014bis','C3:BD56');
        OmegaRaw(isnan(OmegaRaw))=0;
        OmegaRaw((OmegaRaw)<0)=0;
        OmegaRaw=OmegaRaw';
        
        pxij = OmegaRaw;
        pmXi = sum(pxij,2);
        
    %total output 
        pyi = sum(pxij,1)' + pfi;

    %number of sectors
        N=length(pxij);

    %normalized the matrix
        Gamma=zeros(size(pxij));
        for i=1:N
            Gamma(i,:) = pxij(i,:) ./ pyi(i);
        end
  
    
    %Labor income share in Value Added 
        %sectoral data
        empcomp = xlsread('data_deep/WIOT2014.xlsx','2014bis','C67:BD67');
        va = xlsread('data_deep/WIOT2014.xlsx','2014bis','C68:BD68');
        wli_over_vhi = (empcomp./va)';     
      

    %%Load sector description
        %Long description and code
        [sector_code_num,sector_info_txt,sector_info]=xlsread('data_deep/WIOT2014.xlsx','Countries','A46:C99');
        
        sector_code=cell(length(sector_info),1);
        for k=1:length(sector_info)
            if isempty(sector_info_txt{k,1})
                sector_code{k} = num2str(sector_info{k,1});
            else
                sector_code{k}=sector_info_txt{k,1};
            end
             
            %sector_code=sector_info_raw(:,1);
        end
        
        sector_descriptionfrench=sector_info_txt(:,2);
        sector_shortdescription=sector_info_txt(:,3);

        
    toc;
    
     %% II %%%%%%%%%%%%%%%%%%%%%% Calibrate Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('%% II %%%%%%%%%%%%%%%%%%%%%% Calibrate Parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    tic;
    
    %Omega,eta,phi,Delta,psi
    
    %%% Omega_ij = share of j in total intermediate inputs cost by i
    % Omega=zeros(size(pxij));
    % for i=1:N
    %     Omega(i,:) = pxij(i,:) ./ pmXi(i);
    % end
    Omega = diag( pmXi.^(-1) ) * pxij;
    %Adjancy matrix 
        
    %%% eta_i = value added share (of revenue) in sector i
    eta = 1 - pmXi./pyi ;
    
    %%% phi_i = share of final demand i in total output of i
    phi = pfi ./ pyi;
    
    %%% Delta_ij =  expenditure on j by i as a share of total production of j
    % Delta=zeros(size(pxij));
    % for i=1:N
    %     for j=1:N
    %         Delta(i,j) = pxij(i,j) ./ pyi(j);
    %     end
    % end
    Delta = pxij * diag( pyi.^(-1) ) ;
    
    %%% psi_i = share of final demand i in to total final demand
    psi = pfi ./sum( pfi);
    
    %%% gamma_i = labor income share in value added in sector i
    gamma = wli_over_vhi;
    
    %%% sector labor income over GDP
    %Labor income relative to GDP
    va = pyi - pmXi;
    %wli_over_GDP
    lambda = wli_over_vhi.*va./sum(pfi);
    %rli_over_GDP
    rho = (1-wli_over_vhi).*va./sum(pfi);
    
    
    
    toc;
    
  
    
%% IV %%%%%%%%%%%%%%%%%%%%%% Networks Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%% IV %%%%%%%%%%%%%%%%%%%%%% Networks Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
tic;

%Outdegree
outdegree = ones(N,1)'* Gamma; 
outdegree=outdegree'

%Centrality
centrality = (pfi ./sum(pfi))'*inv( eye(N) - Gamma );
centrality=centrality'

%Upstreamness
upstreamness = inv( eye(N) - diag( pyi.^(-1) ) * Gamma' * diag( pyi )  )*ones(N,1);

%Domar weights
domar = pyi./sum(pfi)

tablestats = table(sector_code, sector_shortdescription, outdegree, centrality , upstreamness, domar  );
%disp(tablestats)

filename_export = [ 'export/US_BEA504_network_stats.xlsx' ];
delete(filename_export)
writetable(tablestats,filename_export,'Range','A1')

betak=(pfi ./sum(pfi));    
toc;
%% V %%%%%%%%%%%%%%%%%%%%%% Save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%% V %%%%%%%%%%%%%%%%%%%%%% Save  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
tic;



    %%Data
    Data.N = N;
    Data.pfi = pfi; %final final uses, Total-exports/imports 
    Data.pxij=pxij; %transponse of the input output table containing just intermediate goods
    Data.pmXi=pmXi; % Column vector containing the sum of each row of the transpose of the IOtable
    % It contains the sum of intermediate use of each sector
    Data.pyi=pyi; %Total production per sector 
    Data.wli_over_vhi = wli_over_vhi; %empoloyment component over value added
    Data.va=va %value added by each industry
    Data.Gamma=Gamma;% NxN matrix containing the normalized input/output matrix -> Adjency matrix 

    Data.betak=(pfi ./sum(pfi))

    Data.sector_code=sector_code;
    Data.sector_descriptionfrench=sector_descriptionfrench;
    Data.sector_shortdescription=sector_shortdescription;


    %%Calibrated parameters
    Cali.N = N;
    
    Cali.Omega=Omega; %share of j in total intermediate inputs cost by i
    Cali.eta=eta; %value added share (of revenue) in sector i
    Cali.phi=phi; %share of final demand i in total output of i
    Cali.Delta=Delta; %expenditure on j by i as a share of total production of j
    Cali.psi=psi;%share of final demand i in to total final demand of a good i 
    Cali.gamma=gamma; %labor income share in value added in sector i
    Cali.lambda=lambda; %wli_over_GDP
    Cali.rho=rho; %rli_over_GDP
    
    Kappa = inv(eye(N)-Delta);
    
    table_data = table(sector_code, sector_shortdescription, eta  , phi  , psi , gamma, lambda, rho   );
    
    filename_export = [ 'export/US_BEA504_data.xlsx' ];
    delete(filename_export)
    writetable(table_data,filename_export,'Sheet','Para','Range','A1');    
    writetable(table(sector_code,Delta),filename_export,'Sheet','Delta','Range','A1');
    writetable(table(sector_code,Omega),filename_export,'Sheet','Omega','Range','A1');
    writetable(table(sector_code,Kappa),filename_export,'Sheet','Kappa','Range','A1');

    
    
    %%Networks statistics
    Stats.outdegree=outdegree;
    Stats.centrality=centrality;
    Stats.upstreamness=upstreamness;
    Stats.domar=domar;

    %%Saving
    filename = ['data_matlab/US_BEA504_hat'];
    disp(['save using ' filename])
    save(filename,'Cali','Data','Stats') 

    %%Graph%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('data_matlab/US_BEA504_hat')
    
    %Sector Centrality 
    betabark = Stats.centrality;

    %Domar Weight
    Domark= Stats.domar

    %Sector Code
    %code and long description
    sector_code=Data.sector_code;
    sector_longdescription=Data.sector_descriptionfrench;
    sector_shortdescription=Data.sector_shortdescription;

    %threshold above which ones ignore the links
    thresholds = 0.0266; %Set to 0.2 in baseline

    %maximum number of Sink and Sources sectors
    SinkSourceMax=1; %set to 1 in the source code 
    
    tiOmega=Data.Gamma;

    %define the graph
    OmegaPlot= tiOmega.*(tiOmega>thresholds); %get ride of too small edges

    %define the graph object
    G=digraph(OmegaPlot','omitselfloops'); %construct the directed graph matlab object associated with OmegaPlot (get rid of self loop ie diagonale elements)

    %Find the Sinks Sectors (Highest final demand)
    [C,I]=sort(Data.betak*100); 
    SinksRaw=I(end-SinkSourceMax:end)';
    SinksRaw=[ 19 26 18 31  ];
    SinksRaw=[ 19 26   ];
    
    %Find the Sources Sectors (Higehest upstreamness)
    [CC,II]=sort(Cali.gamma*100);
    SourcesRaw=II(end-SinkSourceMax:end)';
    SourcesRaw = [ 1 2 ];
    SourcesRaw = [ 2  27  25 16 29 30];
    SourcesRaw = [ 2  27  25 16 ];
    SourcesRaw = [ 2  27  ];

    %compute the intersection between the top 5 sink and source sectors
    C = intersect(SinksRaw,SourcesRaw); 

    %find the final Sinks and Sources
    Sinks = setdiff(SinksRaw,C);
    Sources = setdiff(SourcesRaw,C);

    if isempty(Sinks)
    error('Number of Sinks/Sources too low: increase SinkSourceMax')
end


disp(['Sinks sectors: ',num2str(Sinks)])
for k=1:length(Sinks)
    disp([num2str(Sinks(k)),' : ', sector_code{Sinks(k)}, ' : ',  sector_longdescription{Sinks(k)} ])
end

disp('--')
disp(['Sources sectors: ',num2str(Sources)])
for k=1:length(Sources)
    disp([num2str(Sources(k)),' : ', sector_code{Sources(k)}, ' : ',  sector_longdescription{Sources(k)} ])
end


%% Defined a colormap bleu-blanc-rouge

ncolor = 1000; % number    
frenchcolormap= zeros(2*ncolor+1,3) ;
frenchcolormap( : , 1 ) = [ones(1,ncolor), linspace(1,0,ncolor+1) ]';
linecolor=linspace(0,1,ncolor+1);
frenchcolormap( : , 2) = [linecolor(1:end-1), linspace(1,0,ncolor+1) ]';
frenchcolormap( : , 3 ) = [ linspace(0,1,ncolor+1), ones(1,ncolor) ]';

frenchcolormapBlue=frenchcolormap( (ncolor+1):end ,:) ;
frenchcolormapRed=frenchcolormap( 1:ncolor ,:) ;


%% Plot the graph (without VA)
figure('units','normalized','outerposition',[0 0 19/21 1]);
%figure('units','normalized','outerposition',[0 0 2/3 3/3]);
%figure();
    
    LWidths = 4*G.Edges.Weight/max(G.Edges.Weight);%the first erm is set to 6 in the source code
    NodeLabel =reshape(sector_code,1,N);
    NodeLabel{2}=' ';
    NodeLabel{27}= 'CONSTRUCTION';
    NodeLabel{44}= 'REAL ESTATE';
    NodeLabel{29}= 'WHOLESALE TRADE';
    NodeLabel{52}= ' ';
    NodeLabel{53}= 'HUMAN HEALTH';
    NodeLabel{51}= 'PUBLIC ADMINISTRATION';
    NodeLabel{50}= 'ADMINISTRATIVE ACTIVITIES';
    NodeLabel{30}= 'RETAIL TRADE';
    NodeLabel{36}= 'FOOD SERVICES';

    %without color of nodes for VA
    p=plot(G,'LineStyle','-','EdgeColor',[0.4, 0.4, 0.4],'LineWidth',LWidths,'MarkerSize',Domark*250,'NodeColor',[  1    0.20    0.2],'Layout','layered','Sinks',Sinks,'Sources',Sources,'NodeLabel',NodeLabel);
    %p=plot(G,'LineStyle','-','EdgeColor',[0.4, 0.4, 0.4],'LineWidth',LWidths,'MarkerSize',Data.betak*250,'NodeColor',[0.4660 0.6740 0.1880],'Layout','layered','Sources',Sources,'NodeLabel',NodeLabel);
   
    %colormap(flipud(frenchcolormap))
    %colormap((frenchcolormapRed))

    %custom
    p.NodeFontWeight='bold';
    
    
    set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[],'Visible','off')

    

    
print(['fig/IO_network_France.eps'],'-depsc','-tiff')   


%fig = gcf;
%fig.PaperUnits = 'centimeters';
%fig.PaperPosition = [0 0 16.4 11.8 ];
print(['fig/IO_network_France.png'],'-dpng','-r600')  

%% heatmap 
%heatmap(sector_shortdescription,sector_shortdescription,Gamma)
    
end
