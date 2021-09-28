function []=TSA_EQ_fit(data,EQ_UTC,drop_ifgidx,model_type,aps_flag,polygon,wgt_flag)

if nargin>1
    %   Set parameters
    if nargin<7
        wgt_flag='semi';
    end
    
    if nargin<6
        lon1=10.2;
        lon2=11;
        lat1=20.2;
        lat2=19.2;
        
        polygon=[lon1 lon1 lon2 lon2;lat1 lat2 lat2 lat1];
    end

    if nargin<5
        aps_flag='none';
    end

    if nargin<4
        model_type='cidp';
    end
    
    if nargin<3
        drop_ifgidx=[];
    end

    %   Load uph and latlon
    if ~strcmp(aps_flag,'none')
        if contains(model_type,'_mle')
            savefilename=['TSA_mlefit_' aps_flag];
            mleresultname=['mle_result_' aps_flag];
        else
            savefilename=['TSA_fit_' aps_flag];
        end
        uph=data.ifg_aps;
    else
        if contains(model_type,'_mle')
            savefilename='TSA_mlefit';
            mleresultname='mle_result';
        else
            savefilename='TSA_fit';
        end
        uph=data.ifg;
    end

    idx=inpolygon(data.lon,data.lat,polygon(1,:),polygon(2,:));
    lon=data.lon(idx);
    lat=data.lat(idx);
    
    LOS=-uph(idx,:)*1000*0.0555/4/pi;  
    if exist('ref_aps.mat','file') && ~strcmp(aps_flag,'none')
        ref_aps=load('ref_aps.mat');
        LOS=LOS-repmat(ref_aps.ref_LOS',size(LOS,1),1);
    elseif exist('ref.mat','file') && strcmp(aps_flag,'none')
        ref=load('ref.mat');
        LOS=LOS-repmat(ref.ref_LOS',size(LOS,1),1);
    else
        disp('Reference pixels are not selected!');
    end
    LOS(:,drop_ifgidx)=[];
    ym=LOS';

    %   Get SLC acquisition time
    load('slclist');
    slclist(drop_ifgidx)=[];
    n_slc=size(slclist,1);
    slcday=zeros(n_slc,1);
    sat_UTC=load('./parms_aps.mat','UTC_sat');
    sat_UTC=sat_UTC.UTC_sat;
    for i=1:n_slc
        slcday(i)=datenum([slclist{i}(1:8) sat_UTC],'yyyymmddHH:MM');
    end

    %   Load bp and time
    bp=data.bp;
    bp(drop_ifgidx)=[];
    t=data.day;
    t(drop_ifgidx)=[];
    [n_obs,n_pixel]=size(ym);
    yday=365.25;
    year=t/yday;
    year=year+str2num(slclist{1}(1:4))+(str2num(slclist{1}(5:6))-1)/12 ...
        +str2num(slclist{1}(7:8))/yday;   
    
    %   Multi-earthquakes
    n_eq=size(EQ_UTC,1);
    eq_idx=zeros(n_eq,1);
    eq_day=datenum(EQ_UTC,'yyyymmddTHHMM');
    for i=1:n_eq
        eq_idx(i)=sum(slcday<eq_day(i))+1;
    end
    
    %   weight matrix
    if strcmp(wgt_flag,'none')
        P=eye(n_obs);
    else
        % using the semi-varigram fitting results as weight matrix
        load('semi_fit','semi','semi_aps');
        if strcmp(aps_flag,'none')
            cvar=semi(:,2);
        else
            cvar=semi_aps(:,2);
        end
        cvar(drop_ifgidx)=[];
        cvar(1)=mean(cvar(2:end)); % as first epoch is manually set to 0
        cvar=cvar*1e6;             % convert sill unit from m^2 to mm^2
        P=diag(1./cvar);
    end

    %   Modelling
    rms=zeros(1,n_pixel);
    % constant ratio for postseismic decay 1/tau, which applys for the case in Iran
    pcst=0.1784;    
    switch(model_type)
        % Coseismic and interseismic
        case {'ci'}
            A=zeros(n_obs,2+n_eq);
            for i=1:n_eq
                A(eq_idx(i):end,i)=1;
            end
            A(:,end-1)=t;
            A(:,end)=1;

            x=(A'*P*A)\A'*P*ym;
            y=A*x;
            
            Qmm=inv(A'*P*A);
            EQ_def=x(1:n_eq,:);
            sigma=zeros(2+n_eq,n_pixel);
            for i=1:n_pixel
                rms(i)=(y(:,i)-ym(:,i))'*P*(y(:,i)-ym(:,i));
                sigma(:,i)=sqrt(diag(Qmm)*rms(i)/n_obs);
                rms(i)=sqrt(rms(i)/sum(P,'all'));
            end
            mae=sum(abs(y-ym))/n_obs;
            
        % Coseismic, interseismic, and DEM error
        case {'cid'}
            A=zeros(n_obs,3+n_eq);
            for i=1:n_eq
                A(eq_idx(i):end,i)=1;
            end
            A(:,end-2)=t;
            A(:,end-1)=1;
            A(:,end)=bp;

            x=(A'*P*A)\A'*P*ym;
            y=A*x;
            
            Qmm=inv(A'*P*A);
            EQ_def=x(1:n_eq,:);
            sigma=zeros(3+n_eq,n_pixel);
            for i=1:n_pixel
                rms(i)=(y(:,i)-ym(:,i))'*P*(y(:,i)-ym(:,i));
                sigma(:,i)=sqrt(diag(Qmm)*rms(i)/n_obs);
                rms(i)=sqrt(rms(i)/sum(P,'all'));
            end
            mae=sum(abs(y-ym))/n_obs;
            
        % Coseismic, interseismic, DEM error, and seasonal effect
        case {'cids'}
            A=zeros(n_obs,7+n_eq);
            for i=1:n_eq
                A(eq_idx(i):end,i)=1;
            end
            A(:,n_eq+1)=year;
            A(:,n_eq+2)=sin(2*pi*year);
            A(:,n_eq+3)=cos(2*pi*year);
            A(:,n_eq+4)=sin(4*pi*year);
            A(:,n_eq+5)=cos(4*pi*year);
            A(:,n_eq+6)=1;
            A(:,n_eq+7)=bp;
            
            x=(A'*P*A)\A'*P*ym;
            y=A*x;
            
            Qmm=inv(A'*P*A);
            EQ_def=x(1:n_eq,:);
            sigma=zeros(7+n_eq,n_pixel);
            for i=1:n_pixel
                rms(i)=(y(:,i)-ym(:,i))'*P*(y(:,i)-ym(:,i));
                sigma(:,i)=sqrt(diag(Qmm)*rms(i)/n_obs);
                rms(i)=sqrt(rms(i)/sum(P,'all'));
            end
            mae=sum(abs(y-ym))/n_obs;
            
        % Coseismic, interseismic, DEM error, and postseismic
        case {'cidp'}
            n_var=n_eq+4;
            
            A=zeros(n_obs,n_var);
            for i=1:n_eq
                A(eq_idx(i):end,i)=1;
            end
            A(eq_idx(1):end,n_eq+1)=log(1+pcst*(t(eq_idx(1):end) ...
            	-t(eq_idx(1))+slcday(eq_idx(1)+1)-eq_day(1)));
            
            A(:,n_eq+2)=t;
            A(:,n_eq+3)=1;
            A(:,n_eq+4)=bp;
            
            x=(A'*P*A)\A'*P*ym;
            y=A*x;
            
            Qmm=inv(A'*P*A);
            EQ_def=x(1:n_eq,:);
            sigma=zeros(4+n_eq,n_pixel);
            for i=1:n_pixel
                rms(i)=(y(:,i)-ym(:,i))'*P*(y(:,i)-ym(:,i));
                sigma(:,i)=sqrt(diag(Qmm)*rms(i)/n_obs);
                rms(i)=sqrt(rms(i)/sum(P,'all'));
            end
            mae=sum(abs(y-ym))/n_obs;

        %   Apply Maximum likelihood estimation (MLE)
        case {'cidp_mle'}
            rng shuffle
            %% set these parameters manually with caution, and must use 'semi' weight matrix for this estimation
            pidx=[1];   % the index of earthquakes that has postseismic deformation model  
            lbtau=0.01;  % lower boundary of parameter 1/tau
            ubtau=5.0;   % upper boundary of parameter 1/tau
            
            % random walk size for: C(may have multiple values depends on n_eq),A,1/tau,V,b,thita
            %   since 1/tau is set normal distributed ([lbtau,ubtau]) thus its random walk value will be ignored.
            step=[5;5;1;5/yday;5;1e-3]; 
            n_iter=1e6;  % times of iteration
            burn_in=0.2; % ratio for burn_in
            %% 
            %   Set up parameter matrix
            n_pidx=size(pidx,1);
            n_var=3+n_eq+2*n_pidx;
            
            A=zeros(n_obs,n_eq+n_pidx+3);
            for i=1:n_eq
                A(eq_idx(i):end,i)=1;
            end
            for i=1:n_pidx
                A(eq_idx(pidx(i)):end,n_eq+i)=t(eq_idx(pidx(i)):end) ...
                    -t(eq_idx(pidx(i)))+slcday(eq_idx(pidx(i)))-eq_day(pidx(i));
            end
            A(:,n_eq+n_pidx+1)=t;
            A(:,n_eq+n_pidx+2)=1;
            A(:,n_eq+n_pidx+3)=bp;
           
            %   Set up modelfun
            funstr='@(b,x)';
            for i=1:n_eq
                funstr=[funstr 'x(:,' num2str(i) ')*b(' num2str(i) ')+'];
            end
            for i=1:n_pidx
                funstr=[funstr 'x(:,' num2str(pidx(i)) ')*b(' num2str(n_eq+2*i-1) ...
                    ').*log(1+x(:,' num2str(n_eq+i) ')*b(' num2str(n_eq+2*i) '))+'];
            end
            funstr=[funstr 'x(:,' num2str(n_eq+n_pidx+1) ')*b(' num2str(n_var-2) ')+x(:,' ...
                num2str(n_eq+n_pidx+2) ')*b(' num2str(n_var-1) ')+x(:,' num2str(n_eq+n_pidx+3) ')*b(' num2str(n_var) ')'];
            modelfun=eval(funstr);
            
            %   First use a simpler linear model to provide intial result
            B=zeros(n_obs,3+n_eq);
            B(:,1:n_eq)=A(:,1:n_eq);
            B(:,end-2:end)=A(:,end-2:end);
            x0=(B'*P*B)\B'*P*ym;
            
            %   Prepare for inversing pixels one by one
            n_pixel=size(ym,2); 
            x=zeros(n_var,n_pixel);
            y=zeros(n_obs,n_pixel);
            %   Set values that do not in the loop
            m_all=zeros(n_var,n_pixel);
            m_all(1:n_eq,:)=x0(1:n_eq,:);
            m_all(n_eq+1:2:n_eq+2*n_pidx,:)=0;
            m_all(n_eq+2:2:n_eq+2*n_pidx,:)=(lbtau+ubtau)/2;
            m_all(end-2:end,:)=x0(end-2:end,:);
                
            m_mean=zeros(n_var,n_pixel);
            m_std=zeros(n_var,n_pixel);
            m_median=zeros(n_var,n_pixel);
            ratio=zeros(n_pixel,1);
            
            disp(['There are total ' num2str(n_pixel) ' pixels need to be processed!']);
            for i=1:n_pixel
                if i/1000==floor(i/1000)
                    disp([num2str(i) ' pixels have been processed.']);
                end
                
                %   MLE
                m=m_all(:,i);
                j_accept=0;
                m_keep=zeros(n_var,n_iter);
                log_likelihood_keep=zeros(n_iter,1);
                
                for j=1:n_iter
                    m_trial=m+step.*(2*rand(n_var,1)-1);
                    m_trial(n_eq+2:2:n_eq+2*n_pidx)=(ubtau-lbtau)*rand(n_pidx)+lbtau;
                    
                    y_trial=modelfun(m_trial,A);
                    r=ym(:,i)-y_trial;
                    log_likelihood_trial=-(r'*P*r)/2;

                    y_current=modelfun(m,A);
                    r=ym(:,i)-y_current;
                    log_likelihood_current=-(r'*P*r)/2;

                    ratio_likelihoods=log_likelihood_trial-log_likelihood_current;

                    if exp(ratio_likelihoods)>rand(1,1)
                        m=m_trial;
                        j_accept=j_accept+1;
                    end
                    m_keep(:,j)=m;
                    log_likelihood_keep(j)=log_likelihood_current;
                end

                x(:,i)=m_keep(:,end);
                y(:,i)=modelfun(m_keep(:,end),A);
                ratio(i)=j_accept/n_iter;
                
                m_keep=m_keep(:,n_iter*burn_in+1:end);
                m_mean(:,i)=mean(m_keep,2);
                m_std(:,i)=std(m_keep,0,2);
                m_median(:,i)=median(m_keep,2);
            end
            rms=sqrt(sum((y-ym).^2)/n_obs);
            mae=sum(abs(y-ym))/n_obs;
            EQ_def=x(1:n_eq,:);
            
            save(mleresultname,'m_mean','m_std','m_median','rms','mae');
    end
    
    save(savefilename,'year','x','y','ym','EQ_def','rms','mae','lon','lat','slcday','eq_day','eq_idx','bp','model_type'); 
    if ~contains(model_type,'_mle')
        save(savefilename,'sigma','-append')
    end
else
    load(data);
    if contains(data,'_gacos')
        aps_flag='gacos';
    elseif contains(data,'_powerlaw')
        aps_flag='powerlaw';
    elseif contains(data,'_linear')
        aps_flag='linear';
    else
        aps_flag='none';
    end
    savefilename=data;
end

n_eq=size(eq_idx,1);
for i=1:n_eq
   figure;
   scatter(lon(:),lat(:),3,'filled','cdata',EQ_def(i,:));
   colorbar;
   title(['Coseismic Deformatioin Field of Earthquake ' num2str(i) ' (unit: mm)']);
   % set the color scale manually here for better plotting
   switch(i)
       case 1
           set(gca,'clim',[-200 200])
       case 2
           set(gca,'clim',[-30 30])
       case 3
           set(gca,'clim',[-30 30])
       case 4
           set(gca,'clim',[-30 30])
   end
   
   mButton=uicontrol('Style', 'pushbutton', 'Callback',...
        {@TSA_EQ_pixel,aps_flag,savefilename,model_type},'String','TSA_EQ', 'Position', [0 0 80 20]);
end
disp('Done.');

end