function fage_out_decade=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
    nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
    dist_scen_start,dist_scen_end,nyout,nages_dec,firstscenyear,lastscenyear)
%Make the calculations of how forest age changes over time.
%To be called from age_class_reconstruction.m
%Track primary and secondary forest from LUHv2h separately and only merge together in the final output array.
%
%T. Pugh
%12.09.19

firstscenind=firstscenyear-year1;
lastscenind=lastscenyear-year1;

%If using a disturbance rate scenario in time, initialise the annual multiplicative array here
if use_dist_scen
    dist_scen=NaN(nyear,1);
    for nn=1:firstscenind-1
        dist_scen(nn)=dist_scen_start;
    end
    for nn=firstscenind:lastscenind
        df=dist_scen_end-dist_scen_start;
        nnn=nn-firstscenind;
        dist_scen(nn)=dist_scen_start+(nnn*(df/(lastscenind-firstscenind)));
    end
    for nn=lastscenind+1:nyear
        dist_scen(nn)=dist_scen_end;
    end
    clear nn
end

%Set weight arrays for likelihood of LUC or wood harvest by age class
if LUC_age_weights==0
    %Do nothing, no age weighting applied
elseif LUC_age_weights==1
    %Increased likelihood of LUC for older classes (arbitrary increase to test sensitivity)
    weights_LUC=1:(4/(nages-1)):5;
elseif LUC_age_weights==2
    %Increased likelihood of LUC for younger classes (arbitrary increase to test sensitivity)
    weights_LUC=5:-(4/(nages-1)):1;
else
    error('LUC_age_class not set to 0, 1 or 2')
end
if LUC_age_weights>0 && length(weights_LUC)>nages
    error('length(weights_LUC)>nages')
end

%Set weight arrays for likelihood of disturbance by age class
if dist_age_weights==0
    %Do nothing, no age weighting applied
elseif dist_age_weights==1
    %Increased likelihood of disturbance for older classes (arbitrary increase to test sensitivity)
    weights_dist=1:(4/(nages-1)):5;
elseif dist_age_weights==2
    %Increased likelihood of disturbance for younger classes (arbitrary increase to test sensitivity)
    weights_dist=5:-(4/(nages-1)):1;
else
    error('dist_age_class not set to 0, 1 or 2')
end
if dist_age_weights>0 && length(weights_dist)>nages
    error('length(weights_dist)>nages')
end

%Settings for desired output years
fage_out_decade=NaN(360,180,nages_dec,nyout);

%Initialise diagnostic arrays to track whether there are any unallocated transitions for the LUH2 data
prim_unallocated=zeros(360,180);
sec_unallocated=zeros(360,180);

%Loop over latitudes to make the array size manageable in memory
llint=20; %index range for latitude bands

for lls=1:llint:180
    fprintf('%d\n',lls)
    
    %Initialise forest fractions at the maximum age (i.e. all old-growth prior to the beginning of the nage-long spin-up period)
    fage_sec=zeros(360,llint,nages);
    fage_sec(:,:,nages)=secdf_init(:,lls:lls+llint-1);
    fage_prim=zeros(360,llint,nages);
    fage_prim(:,:,nages)=primf_init(:,lls:lls+llint-1);
    
    fage=zeros(360,llint,nages,nyear+nspinup);
    
    %Now start the year loop for forest age calculation
    for yy=1:nyear+nspinup
        
        %First increment the ages of forest by one year
        %Secondary forest
        fage_sec(:,:,nages)=fage_sec(:,:,nages)+fage_sec(:,:,nages-1); %Adding to old-growth
        fage_sec(:,:,2:nages-1)=fage_sec(:,:,1:nages-2); %All intermediate ages increment by one age class
        fage_sec(:,:,1)=zeros(360,llint); %First age class must become zero
        %Primary forest
        fage_prim(:,:,nages)=fage_prim(:,:,nages)+fage_prim(:,:,nages-1);
        fage_prim(:,:,2:nages-1)=fage_prim(:,:,1:nages-2);
        fage_prim(:,:,1)=zeros(360,llint);
        
        if inc_luh2
            %Subtract forest loss
            %Remove losses from a randomly-chosen age class until all losses are allocated
            
            for ii=1:360
                for jj=1:llint %Index within current latitude band
                    jjj=lls-1+jj; %Index across all latitude bands
                    
                    %Secondary forest
                    if luh2_forlu_loss_sec_1deg(ii,jjj,yy)>0
                        to_lose=luh2_forlu_loss_sec_1deg(ii,jjj,yy); %Initialise running total of lost forest fraction still to be allocated
                        while to_lose>0.00000001
                            hf=find(fage_sec(ii,jj,:)>0); %Only select the ages with forest cover in order to optimise the calculation
                            if ~isempty(hf)
                                sec_carryover=0;
                                
                                if LUC_age_weights==0
                                    %Randomly choose age class to remove forest from, weighted by fraction of forest in age class to avoid abnormally low probabilities for old growth
                                    rr=randsample(length(hf),1,true,fage_sec(ii,jj,hf));
                                elseif LUC_age_weights==1 || LUC_age_weights==2
                                    %Randomly choose age class with modified likelihood by age class
                                    rr=randsample(length(hf),1,true,weights_LUC(hf).*squeeze(fage_sec(ii,jj,hf))');
                                end
                                
                                if fage_sec(ii,jj,hf(rr)) > to_lose
                                    fage_sec(ii,jj,hf(rr))=fage_sec(ii,jj,hf(rr))-to_lose;
                                    break
                                else
                                    to_lose=to_lose-fage_sec(ii,jj,hf(rr));
                                    fage_sec(ii,jj,hf(rr))=0;
                                end
                            else
                                %If there is no secondary forest to lose then ignore this loss
                                sec_unallocated(ii,jjj)=sec_unallocated(ii,jjj)+to_lose;
                                sec_carryover=to_lose; %Attempt to carryover unallocated losses to primary forest
                                break
                            end
                        end
                        clear to_lose rr
                    end
                    
                    %Primary forest
                    if luh2_forlu_loss_prim_1deg(ii,jjj,yy)>0
                        to_lose=luh2_forlu_loss_prim_1deg(ii,jjj,yy)+sec_carryover; %Initialise running total of lost forest fraction still to be allocated
                        while to_lose>0.00000001
                            hf=find(fage_prim(ii,jj,:)>0);
                            if ~isempty(hf)
                                if LUC_age_weights==0
                                    %Randomly choose age class to remove forest from, weighted by fraction of forest in age class to avoid abnormally low probabilities for old growth
                                    rr=randsample(length(hf),1,true,fage_prim(ii,jj,hf));
                                elseif LUC_age_weights==1 || LUC_age_weights==2
                                    %Randomly choose age class, with modified likelihood by age class
                                    rr=randsample(length(hf),1,true,weights_LUC(hf).*squeeze(fage_prim(ii,jj,hf))');
                                end
                                
                                if fage_prim(ii,jj,hf(rr)) > to_lose
                                    fage_prim(ii,jj,hf(rr))=fage_prim(ii,jj,hf(rr))-to_lose;
                                    break
                                else
                                    to_lose=to_lose-fage_prim(ii,jj,hf(rr));
                                    fage_prim(ii,jj,hf(rr))=0;
                                end
                            else
                                %If there is no primary forest to lose then ignore this loss
                                prim_unallocated(ii,jj)=prim_unallocated(ii,jj)+to_lose;
                                break
                            end
                        end
                        clear to_lose rr
                    end
                end
                clear jj jjj
            end
            clear ii
            
            %Add forest gain to youngest secondary age class
            fage_sec(:,:,1)=fage_sec(:,:,1)+luh2_forlu_gain_1deg(:,lls:lls+llint-1,yy);
        end
        
        if inc_dist
            %Carry out random disturbance (equal probability across all ages of secondary and primary forest) and add
            %to youngest age class of secondary forest.
            %Disturb a fixed fraction per year defined by 1/distint
            
            if use_dist_scen
                %Modify the disturbance rate
                if yy<=nspinup
                    distrate=(1./distint(:,lls:lls+llint-1))*dist_scen(1);
                else
                    distrate=(1./distint(:,lls:lls+llint-1))*dist_scen(yy-nspinup);
                end
            else
                distrate=(1./distint(:,lls:lls+llint-1));
            end
            
            if dist_age_weights==0
                %If applying a uniform disturbance then can do this as a simple multiplication across the arrays
                
                %Secondary forest
                frac_dist_sec=fage_sec(:,:,2:nages).*repmat(distrate,[1 1 nages-1]);
                fage_sec(:,:,2:nages)=fage_sec(:,:,2:nages)-frac_dist_sec; %Do not disturb youngest age class
                fage_sec(:,:,1)=fage_sec(:,:,1)+sum(frac_dist_sec,3);
                clear frac_dist_sec
                %Primary forest
                frac_dist_prim=fage_prim(:,:,2:nages).*repmat(distrate,[1 1 nages-1]);
                fage_prim(:,:,2:nages)=fage_prim(:,:,2:nages)-frac_dist_prim;
                fage_prim(:,:,1)=fage_prim(:,:,1)+sum(frac_dist_prim,3);
                clear frac_dist_prim
                
            elseif dist_age_weights==1 || dist_age_weights==2
                %Carry out disturbance with likelihood weighted by age class
                
                %Calculate total area to disturb
                to_dist_prim=sum(fage_prim(:,:,2:nages).*repmat(distrate,[1 1 nages-1]),3); 
                to_dist_sec=sum(fage_sec(:,:,2:nages).*repmat(distrate,[1 1 nages-1]),3);
                
                %Remove from random age class until all disturbances are allocated
                for ii=1:360
                    for jj=1:llint
                        
                        if sum(fage_prim(ii,jj,:))>0
                            %Primary forest
                            hf=find(fage_prim(ii,jj,:)>0);
                            if ~isempty(hf)
                                while to_dist_prim(ii,jj)>0.00000001
                                    %Randomly choose age class, with increased likelihood for older classes, weighted by fraction of forest in age class to avoid abnormally low probabilities for old growth
                                    rr=randsample(length(hf),1,true,weights_dist(hf).*squeeze(fage_prim(ii,jj,hf))');
                                    
                                    if fage_prim(ii,jj,hf(rr)) > to_dist_prim(ii,jj)
                                        fage_prim(ii,jj,hf(rr))=fage_prim(ii,jj,hf(rr))-to_dist_prim(ii,jj);
                                        fage_prim(ii,jj,1)=fage_prim(ii,jj,1)+to_dist_prim(ii,jj);
                                        break
                                    else
                                        to_dist_prim(ii,jj)=to_dist_prim(ii,jj)-fage_prim(ii,jj,hf(rr));
                                        fage_prim(ii,jj,1)=fage_prim(ii,jj,1)+fage_prim(ii,jj,hf(rr));
                                        fage_prim(ii,jj,hf(rr))=0;
                                    end
                                end
                            end
                            clear rr hf
                        end
                        
                        if sum(fage_sec(ii,jj,:))>0
                            %Secondary forest
                            hf=find(fage_sec(ii,jj,:)>0);
                            if ~isempty(hf)
                                while to_dist_sec(ii,jj)>0.00000001
                                    %Randomly choose age class, with increased likelihood for older classes, weighted by fraction of forest in age class to avoid abnormally low probabilities for old growth
                                    rr=randsample(length(hf),1,true,weights_dist(hf).*squeeze(fage_sec(ii,jj,hf))');
                                    
                                    if fage_sec(ii,jj,hf(rr)) > to_dist_sec(ii,jj)
                                        fage_sec(ii,jj,hf(rr))=fage_sec(ii,jj,hf(rr))-to_dist_sec(ii,jj);
                                        fage_sec(ii,jj,1)=fage_sec(ii,jj,1)+to_dist_sec(ii,jj);
                                        break
                                    else
                                        to_dist_sec(ii,jj)=to_dist_sec(ii,jj)-fage_sec(ii,jj,hf(rr));
                                        fage_sec(ii,jj,1)=fage_sec(ii,jj,1)+fage_sec(ii,jj,hf(rr));
                                        fage_sec(ii,jj,hf(rr))=0;
                                    end
                                end
                            end
                            clear rr hf
                        end
                    end
                end
                clear to_dist_prim to_dist_sec
            end
                
        end
        
        %Sum both primary and secondary forest in output array
        fage(:,:,:,yy)=fage_sec+fage_prim;
    end
    
    %---
    %Postprocessing of fage
    
    %Select desired years
    output_years_ind=output_years-year1+nspinup;
    fage_out=NaN(360,llint,nages,nyout);
    for nn=1:nyout
        fage_out(:,:,:,nn)=fage(:,:,:,output_years_ind(nn));
    end
    clear nn
    
    %Group age classes to decadal timesteps
    for nn=1:nyout
        dd=0;
        for aa=1:10:nages
            dd=dd+1;
            fage_out_decade(:,lls:lls+llint-1,dd,nn)=sum(fage_out(:,:,aa:aa+9,nn),3);
        end
        clear dd aa
    end
    clear nn
    
end
clear lls llint
clear fage fage_out fage_prim fage_sec