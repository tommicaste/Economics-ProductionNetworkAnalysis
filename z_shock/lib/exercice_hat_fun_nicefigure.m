%%%%%% Copyright Basile Grassi 2020
%%%%%% Compute effect of a labor supply shocks
%%%%% and with (and without) demand shocks



function exercice_hat_fun_nicefigure(theta,sigma,epsilon,cap, califile ,excelfile,fig)


%% I %%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%% I %%%%%%%%%%%%%%%%%%%%%% Load Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    tic;

    filename = ['data_matlab/' califile ];
    load(filename)
    disp(['load: ' filename])

    %%Data
    sector_code=Data.sector_code;
    sector_shortdescription=Data.sector_shortdescription;

    %%Calibrated parameters
    N=Cali.N;
    
    Omega=Cali.Omega;
    eta=Cali.eta;
    phi=Cali.phi;
    Delta=Cali.Delta;
    psi=Cali.psi;
    gamma=Cali.gamma;
    w=Cali.w;
    toc;


    %%turning off some warning
    warning('off','MATLAB:xlswrite:AddSheet')
    

%% III %%%%%%%%%%%%%%%%%%%%%% Compute Effectt of Shock  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%% III %%%%%%%%%%%%%%%%%%%%%% Compute Effect of Shock  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    tic;

    disp(['Using shocks from', excelfile,'.xlsx' ] )
    
    options = optimset('Display', 'Final');
    
    Alphabet = 'BCDEFGHIJKLMNOPQRSTUVWXYZ';
    
    %load the list of supply shocks
    [crap, header] = xlsread(['data_deep/' excelfile '_laborsupply.xlsx'],'A1:ZZ1');
    confinement = header(2:end);
    
    %load the list of demand shocks
    [crap, header] = xlsread(['data_deep/' excelfile '_demand.xlsx'],'A1:ZZ1');
    demand_header = header(2:end);
    
    if nargin>6 && fig==1
        figure('units','normalized','outerposition',[0 0 0.6 0.6]);
        Legend=cell(length(confinement)+1,1);%  two positions ;
        
    end
    
    for c = 1:length(confinement)
        
        for d = 0:length(demand_header)
        
            disp(['------',  confinement{c} ,'----------'])
            if d >0
                disp(['------ demand:',  demand_header{d} ,'----------'])
            end
            
            %Load the shock
            share_shock =  xlsread(['data_deep/' excelfile '_laborsupply.xlsx'],[Alphabet(c) '2:' Alphabet(c) num2str(2+N-1)]);
            shocks = min(share_shock,cap);
            
            if d>0
                demand_shock = xlsread(['data_deep/' excelfile '_demand.xlsx'],[Alphabet(d) '2:' Alphabet(d) num2str(2+N-1)]);
            end
            
            
            %filename for export
            if d==0
                filename_export = [ 'export/' califile '_' confinement{c} '_theta' , strrep(num2str(theta) ,'.','d') , '_sigma' , strrep(num2str(sigma),'.','d') '_epsilon' strrep(num2str(epsilon),'.','d') , '_N' , num2str(N), '.xlsx' ];
                figname = [ califile '_theta' , strrep(num2str(theta) ,'.','d') , '_sigma' , strrep(num2str(sigma),'.','d') '_epsilon' strrep(num2str(epsilon),'.','d') , '_N' , num2str(N) ];
                disp(['export in the file: ' filename_export])            
            elseif d>0
                filename_export = [ 'export/' califile '_' confinement{c} '_' demand_header{d}  '_theta' , strrep(num2str(theta) ,'.','d') , '_sigma' , strrep(num2str(sigma),'.','d') '_epsilon' strrep(num2str(epsilon),'.','d') , '_N' , num2str(N), '.xlsx' ];
                figname = [ califile '_theta' , strrep(num2str(theta) ,'.','d') , '_sigma' , strrep(num2str(sigma),'.','d') '_epsilon' strrep(num2str(epsilon),'.','d') , '_N' , num2str(N) ];
                disp(['export in the file: ' filename_export])            
            end
            
            %delete if exist
            delete(filename_export)
            
            %Name of results in the structure
            name_struct = strrep(confinement{c},' ','_' );
            name_struct = strrep(name_struct,'(','' );
            name_struct = strrep(name_struct,')','' );

            if d>0
                name_demand_struct = strrep(demand_header{d},' ','_' );
                name_demand_struct = strrep(name_demand_struct,'(','' );
                name_demand_struct = strrep(name_demand_struct,')','' );

                name_struct = [name_struct,'_', name_demand_struct ];
            end
            disp(name_struct)
            
            %chnage in capital stock and productivity
            zi_hat = ones(N,1);
            ki_hat = ones(N,1);
            zi_hat(56)= 1- 0.39; %Gas
            zi_hat(57)= 1- 0.24; %Oil
            
            %demand shocks
            betai_hat =ones(N,1);            
            if d>0
                betai_hat=demand_shock;
            end
            
            
            %Solve for  shocked equilibrium 
            disp('Solve for the effect of the shock under Cobb-Douglas')            
                fun_solve_CD_Shocked= @(lvec)    RES_hat_eq(lvec,zi_hat,(1-shocks) ,ki_hat,betai_hat,1,1,1,N, Omega,eta,gamma,phi,Delta,psi,w) ;
                lvec_eq_CD_shocked = fsolve(fun_solve_CD_Shocked,zeros(2*N+1,1),options); 
                %[RES_shocked,p1MoinsThetai_shocked,pThetaYi_shocked,PSigmaY_shocked,ri_shocked,wi_shocked,vi_shocked,pi_shocked,si_shocked,fi_shocked,yi_shocked,Pmi_shocked,Xi_shocked,xij_shocked,hi_shocked,vai_shocked,rki_shocked,ki_shocked,wli_shocked,li_shocked]  
                [RES_CD,PSigmaY_hat_CD,Equi_CD_hat] =fun_solve_CD_Shocked(lvec_eq_CD_shocked);
                %GDP_CD_hat = PSigmaY_hat_CD;
                GDP_CD_hat= Equi_CD_hat.GDPIta;

            disp('---------------')
            disp('Solve for the effect of the shock')            
                fun_solve_Shocked= @(lvec)    RES_hat_eq(lvec,zi_hat, (1-shocks) ,ki_hat,betai_hat,theta,sigma,epsilon,N, Omega,eta,gamma,phi,Delta,psi,w) ;
                lvec_eq_shocked = fsolve(fun_solve_Shocked,lvec_eq_CD_shocked,options); 
                %[RES_shocked,p1MoinsThetai_shocked,pThetaYi_shocked,PSigmaY_shocked,ri_shocked,wi_shocked,vi_shocked,pi_shocked,si_shocked,fi_shocked,yi_shocked,Pmi_shocked,Xi_shocked,xij_shocked,hi_shocked,vai_shocked,rki_shocked,ki_shocked,wli_shocked,li_shocked]  
                [RES,PSigmaY_hat,Equi_hat] =fun_solve_Shocked(lvec_eq_shocked);
                %GDP_hat = PSigmaY_hat;
                GDP_hat= Equi_hat.GDPIta;
                


            %disp('---------------')
            %disp(['GDP'   ])                
                
                GDP_change = GDP_hat -1;
                GDP_change_CD = GDP_CD_hat-1;
                
                GDP_logdiff = log(GDP_hat );
                GDP_CD_logdiff = log( GDP_CD_hat);                
                
                table_general = table(sigma, theta, epsilon, GDP_hat, GDP_change , GDP_CD_hat, GDP_change_CD, GDP_logdiff, GDP_CD_logdiff );
                writetable(table_general,filename_export,'Sheet','General','Range','A1')
                disp(table_general)
                
                Output.(name_struct).general=table_general;

 
            %disp('-------------------')
            %disp('Sector Domar Weight:----')
                domar_hat = Equi_hat.domar;
                domar_change = domar_hat -1;

                domar_CD_hat = Equi_CD_hat.domar;
                domar_CD_change = domar_CD_hat -1;


                domar_tab = table( sector_code, sector_shortdescription , domar_CD_hat,  domar_CD_change, domar_hat, domar_change);
                writetable(domar_tab,filename_export,'Sheet','Domar','Range','A1')
                
                Output.(name_struct).domar=domar_tab;


%             %disp('-------------------')
%             %disp('Sector Sales:----')

                sales_hat = Equi_hat.pyi;
                sales_change = sales_hat -1;

                sales_CD_hat = Equi_CD_hat.pyi;
                sales_CD_change = sales_CD_hat -1;


                sales_tab = table( sector_code, sector_shortdescription , sales_CD_hat,  sales_CD_change, sales_hat, sales_change);

                writetable(sales_tab,filename_export,'Sheet','Sales','Range','A1')
                
                Output.(name_struct).sales=sales_tab;

                
                
%             %disp('------')
%             %disp('Final Demand:----')                
                
                final_hat = Equi_hat.pfi;
                final_change = final_hat -1;

                final_CD_hat = Equi_CD_hat.pfi;
                final_CD_change = final_CD_hat -1;


                final_tab = table( sector_code, sector_shortdescription , final_CD_hat,  final_CD_change, final_hat, final_change);

                writetable(final_tab,filename_export,'Sheet','Final','Range','A1')
                
                Output.(name_struct).final=final_tab;




            %disp('------')
            %disp('Labor Income :----')                
                laborIncome_hat = Equi_hat.wli;
                laborIncome_change = laborIncome_hat -1;

                laborIncome_CD_hat = Equi_CD_hat.wli;
                laborIncome_CD_change = laborIncome_CD_hat -1;

                laborIncome_tab = table( sector_code, sector_shortdescription , laborIncome_CD_hat,  laborIncome_CD_change, laborIncome_hat, laborIncome_change);

                writetable(laborIncome_tab,filename_export,'Sheet','laborIncome','Range','A1')
                
                
                Output.(name_struct).laborIncome=laborIncome_tab;
                
             %disp('Capital Income :----')                
                capitalIncome_hat = Equi_hat.rki;
                capitalIncome_change = capitalIncome_hat -1;

                capitalIncome_CD_hat = Equi_CD_hat.rki;
                capitalIncome_CD_change = capitalIncome_CD_hat -1;


                capitalIncome_tab = table( sector_code, sector_shortdescription , capitalIncome_CD_hat,  capitalIncome_CD_change, capitalIncome_hat, capitalIncome_change);

                writetable(capitalIncome_tab,filename_export,'Sheet','capitalIncome','Range','A1')
                
                Output.(name_struct).capitalIncome=capitalIncome_tab;

            %disp('------')    
            %disp('Rental Rate :----')                
                rate_hat = Equi_hat.ri;
                rate_change = rate_hat -1;

                rate_CD_hat = Equi_CD_hat.ri;
                rate_CD_change = rate_CD_hat -1;

                rate_tab = table( sector_code, sector_shortdescription , rate_CD_hat ,  rate_CD_change, rate_hat, rate_change);

                writetable(rate_tab,filename_export,'Sheet','rate','Range','A1')
                
                Output.(name_struct).rate=rate_tab;

            %disp('------')    
            %disp('Salaire :----')
                wage_hat = Equi_hat.wi;
                wage_change = wage_hat -1;

                wage_CD_hat = Equi_CD_hat.wi;
                wage_CD_change = wage_CD_hat -1;


                wage_tab = table( sector_code, sector_shortdescription , wage_CD_hat , wage_CD_change, wage_hat, wage_change);

                writetable(wage_tab,filename_export,'Sheet','wage','Range','A1')
                
                Output.(name_struct).wage=wage_tab;


            %disp('------')
            %disp('VA :----')
                va_hat = Equi_hat.vhi;
                va_change = va_hat -1;

                va_CD_hat = Equi_CD_hat.vhi;
                va_CD_change = va_CD_hat -1;


                va_tab = table( sector_code, sector_shortdescription , va_CD_hat,  va_CD_change, va_hat, va_change);

                writetable(va_tab,filename_export,'Sheet','va','Range','A1')

                Output.(name_struct).va=va_tab;



            %disp('------')
            %disp('Prices :----')                
                price_hat = Equi_hat.pi;
                price_change = price_hat -1;

                price_CD_hat = Equi_CD_hat.pi;
                price_CD_change = price_CD_hat -1;
                
                price_Dlog = log(price_hat);


                price_tab = table( sector_code, sector_shortdescription , price_CD_hat,  price_CD_change, price_hat, price_change, price_Dlog);

                writetable(price_tab,filename_export,'Sheet','price','Range','A1')

                Output.(name_struct).va=va_tab;


            %disp('------')
            %disp('Quantity :----')
                quantity_hat = Equi_hat.yi;
                quantity_change = quantity_hat -1;

                quantity_CD_hat = Equi_CD_hat.yi;
                quantity_CD_change = quantity_CD_hat -1;
                
                quantity_Dlog = log(quantity_hat);
                quantity_CD_Dlog = log(quantity_CD_hat);


                quantity_tab = table( sector_code, sector_shortdescription , quantity_CD_hat,  quantity_CD_change, quantity_hat, quantity_change, quantity_Dlog);

                writetable(quantity_tab,filename_export,'Sheet','quantity','Range','A1')
                
                Output.(name_struct).quantity=quantity_tab;
                
                
            if nargin>6 && fig==1
                
                 
             
                 
             %%%%%%%%%%%%%%%%%%%%%%%%    
             %figure;
                if strcmp(confinement{c},'share_kid_constrained')
                    Legend{c}= 'School Closure';
                    color_fig = [ 0 0 1 ];
                elseif strcmp(confinement{c},'share_closednocomputerkid')
                    Legend{c}= 'Administrative + School';
                    color_fig = [0.33 0.34 0.33];
                elseif strcmp(confinement{c},'share_closednocomputer')
                    Legend{c}= 'Administrative Closure';
                    color_fig = [1 0 0];
                else
                    Legend{c+1}=confinement{c};
                end
             
                subplot(121)
                    [f,xi] =ksdensity(quantity_Dlog);
                    h=plot(xi,f,'-','LineWidth',2);
                    %h.Color = color_fig;
                    hold on;

                    
                
                subplot(122)
                    [f,xi] =ksdensity(quantity_CD_Dlog);
                    h=plot(xi,f,'-','LineWidth',2);
                    %h.Color = color_fig;
                    hold on;             
                
                
                
                    
                
                 
            end
        
        
    
    
        end
    end

%theta,sigma,epsilon,cap, califile ,excelfile

if nargin>6 && fig==1

     %load data
    dlq_4 = xlsread('data_deep/Data_quantity.xlsx','B2:B83');
    
    subplot(121)    
    
        xlim([-1, 0.5])
        ylim([0, 55])
        
        %set(gca,'ytick',[])
        %set(gca,'yticklabel',[])
    
        legend(Legend,'Location','NorthWest');
        xlabel('$ \log{} \widehat{y}_i$','Interpreter','latex'); 
        title( [ '($\sigma$ , $\theta$,  $\varepsilon$) = (', num2str(sigma), ' , ', num2str(theta), ' , ',num2str(epsilon) ,')'], 'interpreter', 'latex' )
    
    subplot(122)
    
    
        xlim([-1, 0.5])
        ylim([0, 55])
    
        %set(gca,'ytick',[])
        %set(gca,'yticklabel',[])
    
        xlabel('$ \log{} \widehat{y}_i$','Interpreter','latex'); 
        title( [ 'Cobb-Douglas ($\sigma$ , $\theta$,  $\varepsilon$) = (1,1,1)'], 'interpreter', 'latex' )
    
    
    print(['fig/', figname , '_qty_hist_nice.png'],'-dpng')  
    print(['fig/' figname  '_qty_hist_nice.eps'],'-depsc','-tiff')  
end            
             

Output.params.theta=theta;
Output.params.sigma=sigma;
Output.params.epsilon=epsilon;
Output.params.cap=cap;

filename_matlab = ['data_matlab/' califile '_'  excelfile '_cap' , strrep(num2str(cap) ,'.','d')  '_theta' , strrep(num2str(theta) ,'.','d') , '_sigma' , strrep(num2str(sigma),'.','d') '_epsilon' strrep(num2str(epsilon),'.','d') , '_N' , num2str(N),'.mat'];
disp(['Save the data in ' filename_matlab])
save(filename_matlab,'Output')



end