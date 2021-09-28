function[]=TSA_EQ_pixel(src,evt,aps_flag,fitfilename,model_type)
% radius of selecting pixels, unit is degree
radius=0.01;
% whether to subtract the DEM errors during plotting for better illustration
subtractDEM=1;
% days in a year
yday=365.25;

fitdata=load(fitfilename);
[lon0,lat0]=ginput(1);
disp(['Selected point coordinates (lon,lat):' num2str(lon0),', ', num2str(lat0) ]);

circle=linspace(0.2*pi,114);
x=radius*cos(circle)+lon0;
y=radius*sin(circle)+lat0;

in=inpolygon(fitdata.lon,fitdata.lat,x,y);
n_ps=sum(in);

if sum(in)==0
    disp('No points found. Please make a new selection or increase the radius.');
    return
end

figure;
year=fitdata.year;
y=fitdata.y(:,in);
ym=fitdata.ym(:,in);
% choose a single pixel to plot
y=y(:,1);
ym=ym(:,1);
x=fitdata.x(:,in);
x=x(:,1);
coseis=fitdata.EQ_def(:,in);
coseis=coseis(:,1);
if contains(fitfilename,'_mlefit')
    if strcmp(aps_flag,'gacos')
        mle_res=load('mle_result_gacos');
    elseif strcmp(aps_flag,'powerlaw')
        mle_res=load('mle_result_powerlaw');
    elseif strcmp(aps_flag,'linear')
        mle_res=load('mle_result_linear');
    else
        mle_res=load('mle_result');
    end
    se=mle_res.m_std(:,in);
    se=se(:,1);
else
    se=fitdata.sigma(:,in);
    se=se(:,1);
end
if subtractDEM==1
    y=y-fitdata.bp*x(end);
end

plot(year,ym,'o',year,y,'*');
set(gcf, 'Position', get(0, 'Screensize'));
title(['Fit Result for Pixel(' num2str(lon0) ',' num2str(lat0) ')']);
xlabel('Year');
ylabel('LOS(mm)');
dim=[.15 .15 .75 .75];

n_var=max(size(x));
n_eq=size(fitdata.eq_idx,1);
str=cell(0);

ratio=tinv(0.95,size(y,1)-size(x,1));
count=0;
if contains(fitfilename,'_mlefit')
    for i=1:n_eq
        str{end+1}=['Coseis_' num2str(i) ': ' num2str(coseis(i)) ' mm'];
        count=count+1;
    end
    switch(model_type)
        case {'cidp_mle'}
            str{end+1}=['A: ' num2str(x(count+1)) ' mm'];
            str{end+1}=['1/tau: ' num2str(x(count+2)) ' day^-^1'];
        case {'cids_mle'}
            str{end+1}=['I: ' num2str(x(count+1))];
            str{end+1}=['phi: ' num2str(x(count+2))];
    end
    
    if size(x,1)==n_eq+5
        str{end+1}=['V: ' num2str(x(count+3)*yday) char(177) num2str(ratio*se(count+3)) ' mm/yr'];
        count=count+1;
    end
    str{end+1}=['Offset: ' num2str(x(count+3)) char(177) num2str(ratio*se(count+3)) ' mm'];
    str{end+1}=['thita: ' num2str(x(count+4)) char(177) num2str(ratio*se(count+4))]; 
else
    for i=1:n_eq
        str{end+1}=['Coseis_' num2str(i) ': ' num2str(coseis(i)) char(177) num2str(ratio*se(i)) ' mm'];
        count=count+1;
    end
    switch(model_type)
        case {'cid'}
            str{end+1}=['V: ' num2str(x(count+1)*yday) char(177) num2str(yday*ratio*se(count+1)) ' mm/yr'];
            count=count+1;
        case {'cdp'}
            str{end+1}=['A: ' num2str(x(count+1)) char(177) num2str(ratio*se(count+1)) ' mm'];
            count=count+1;
        case {'cidp'}
            str{end+1}=['A: ' num2str(x(count+1)) char(177) num2str(ratio*se(count+1)) ' mm'];
            count=count+1;
            str{end+1}=['V: ' num2str(x(count+1)*yday) char(177) num2str(yday*ratio*se(count+1)) ' mm/yr'];
            count=count+1;
        case {'cids'}
            str{end+1}=['V: ' num2str(x(count+1)*yday) char(177) num2str(yday*ratio*se(count+1)) ' mm/yr'];
            count=count+1;
            str{end+1}=['I: ' num2str(x(count+1)) char(177) num2str(ratio*se(count+1)) ' mm'];
            count=count+1;
    end
    str{end+1}=['Offset: ' num2str(x(count+1)) ' mm'];
    str{end+1}=['thita: ' num2str(x(count+2))];
end

annotation('textbox',dim,'String',str,'FitBoxToText','on');

hold on
eq_idx=fitdata.eq_idx;
ymax=max(ym);
ymin=min(ym);
yline=[ymin-10:10:ymax+10];
n_line=size(yline,2);
Legend=cell(3,1);
Legend{1}='Observations';
Legend{2}='Fitting';
Legend{3}='Earthquakes';
delta=3/yday;
for i=1:n_eq
    plot(repmat(year(eq_idx(i))-delta,1,n_line),yline,'--','Color','green');
end
legend(Legend);
end