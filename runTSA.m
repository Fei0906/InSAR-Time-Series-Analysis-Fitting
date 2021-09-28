%% InSAR time series analysis fitting code for sharing

% choose which tropospheric correction you use. if you do not want to apply
% any corrections, set it to 'none'
aps_flag='gacos';

% drop ifgs that you do not wish to keep due to high atmospheric noise or
% unwrapping errors
drop_ifgidx=[];

% read data which contains the following information:
% bp: the perpendicular baseline, reference to the first image
% day: the number of days, reference to the first image
% ifg: unwrapped interferograms (time series data)
% ifg_aps: unwrapped interferograms with tropospheric correction
% lat, lon: lat and lon of pixels
data=load('TSA_ifg_ds');

% Polygon that you wish to run, which define the AOI
lon0=45.30;
lon1=46.50;
lat0=35.25;
lat1=33.70;
polygon=[lon0 lon0 lon1 lon1;lat0 lat1 lat1 lat0];

% the earthquake event time (format: yyyymmddThhmm)
EQ_UTC=['20171112T1818';'20180825T2213';'20181125T1637';'20190106T1341'];

% used model: coseismic, interseismic, DEM errors, and postseismic
% deformations
model_type='cidp';

% variance for individual interferograms. Here set to semi, which uses the
% sill values of semi-varigram fitting as the covariance matrix.
% Note the off-diagonal elements are zeros.
cov_matrix='semi';

% call the function to run
TSA_EQ_fit(data,EQ_UTC,drop_ifgidx,model_type,aps_flag,polygon,cov_matrix);

