%function do_make_isomip_output(hisname,grdname,outname) 
% read in file
if 1 % ocean3
hisname_format = 'ocean_his_';
grdname_format = 'ocean_his_';
outname = 'Ocean3_COM_ROMSUTAS_FISOC.nc'
output_vec=[1:100]; % vector of years
elseif 0 % ocean4
hisname_format = 'ocean_his_';
grdname_format = 'ocean_his_';
outname = 'Ocean4_COM_ROMSUTAS_FISOC.nc'
output_vec=[1:50]; % vector of years
end
Vtransform=2;
Vstretching=4;
theta_s=4.0;
theta_b=0.9;
hc=100;

% define dims
Lpinfo=ncinfo([grdname_format,'0001.nc'],'lon_rho'); Lp = Lpinfo.Size(1);
Mpinfo=ncinfo([grdname_format,'0001.nc'],'lat_rho'); Mp = Mpinfo.Size(2);
L=Lp-1;
M=Mp-1;
time_length=length(output_vec)*12+1; %!!!!!!!!!!!!!!
disp('note: since ROMS is also keeping the initial conditions as the t=1 in the first file, it means that it will be number_of_years*12+1 months of data')
Z_interp=-717.5:5:-2.5; %m
Zp=length(Z_interp);

%% define arrays of compiled data (preallocate for memory savings)
ocean_time_big=nan(time_length,1);
meanMeltRate_big=nan(time_length,1);
meanMeltFlux_big=nan(time_length,1);
sumOfVolume_big=nan(time_length,1);
meanTemperature_big=nan(time_length,1);
meanSalinity_big=nan(time_length,1);
iceDraft_big=nan(Lp,Mp,time_length);
bathymetry_big=nan(Lp,Mp,time_length);
meltRate_big=nan(Lp,Mp,time_length);
Ustar_big=nan(Lp,Mp,time_length);
Tstar_big=nan(Lp,Mp,time_length);
Sstar_big=nan(Lp,Mp,time_length);
v_top_big=nan(Lp,Mp,time_length);
u_top_big=nan(Lp,Mp,time_length);
baroSF_big=nan(Lp,Mp,time_length);
otSF_big=nan(Mp,Zp,time_length);
T_bot_big=nan(Lp,Mp,time_length);
S_bot_big=nan(Lp,Mp,time_length);
tempXZ_big=nan(Zp,Mp,time_length);
saltXZ_big=nan(Zp,Mp,time_length);
tempYZ_big=nan(Zp,Lp,time_length);
saltYZ_big=nan(Zp,Lp,time_length);


fillstart=0;
fillend=0; 

for yy=1:length(output_vec);  % loop over hte number of years


if yy<10,yy_suffix = ['000',num2str(output_vec(yy))];
elseif yy<100 & yy>9,yy_suffix = ['00',num2str(output_vec(yy))];
elseif yy<1000,yy_suffix = ['0',num2str(output_vec(yy))];
end

hisname = [hisname_format,yy_suffix,'.nc'];
grdname = [grdname_format,yy_suffix,'.nc'];
ocean_time=ncread(hisname,'ocean_time'); %load time
disp(['doing ' hisname ])

% Get vertical coords
h=repmat(flipud(ncread(grdname,'h')),[1 1 length(ocean_time)]);
iceDraft=flipud(ncread(grdname,'draft'));
lat_rho=ncread(grdname,'lat_rho');
y_rho=ncread('../Ocean1/isomip_plus_ocean1.nc','y_rho');
x_rho=ncread('../Ocean1/isomip_plus_ocean1.nc','x_rho');
lon_rho=ncread(grdname,'lon_rho');
lat_u=ncread(grdname,'lat_u');
lon_u=ncread(grdname,'lon_u');
lat_v=ncread(grdname,'lat_v');
lon_v=ncread(grdname,'lon_v');
mask_rho=flipud(ncread(grdname,'wetdry_mask_rho'));
%mask_zice=flipud(ncread(grdname,'mask_zice'));
pm=flipud(ncread(grdname,'pm'));
pn=flipud(ncread(grdname,'pn'));
Cs_w=ncread(hisname,'Cs_w');
Cs_r=ncread(hisname,'Cs_r');
mask_zice = mask_rho; mask_draft = (iceDraft<0); mask_zice = mask_zice.*mask_draft;

N = length(Cs_r);
dx = repmat(1./pm,[1 1 length(ocean_time)]);
dy = repmat(1./pn,[1 1 length(ocean_time)]);
mask_closed=ones(size(mask_zice)); mask_closed([1,end],:,:)=NaN; mask_closed(:,[1,end],:)=NaN;
mask_zice(mask_zice==0)=NaN; mask_zice=mask_zice.*mask_closed;
mask_rho(mask_rho==0)=NaN; mask_rho=mask_rho.*mask_closed;
Area_ice = dx.*dy.*mask_zice;
Area_water=dx.*dy.*mask_rho;
rho_i=918;%kgm^-3
Cd=2.5e-3; %dimensionless
u_t=0.01;%m/s
addpath('/ds/projects/iomp/matlab_scripts/ROMS_MATLAB/utility/')
%Zw=nan(size(lat_rho,1),size(lat_rho,2),N+1);
%Z=nan(size(lat_rho,1),size(lat_rho,2),N);
%for jj=1:size(lat_rho,1)
%[z,~,~]=scoord(h,iceDraft,lon_rho,lat_rho,Vtransform,Vstretching,4,0.9,20,N,0,1,jj,0); %z is 3d of depths of each cell.
%[zw,~,~]=scoord(h,iceDraft,lon_rho,lat_rho,Vtransform,Vstretching,4,0.9,20,N,1,1,jj,0); %z is 3d of depths of each cell.
%Z(jj,:,:)=z;
%Zw(jj,:,:)=zw;
%end
Z=nan(length(ocean_time),size(h,1),size(h,2),N);
Zw=nan(length(ocean_time),size(h,1),size(h,2),N+1);
for mini_tt=1:length(ocean_time)
Z(mini_tt,:,:,:)=set_depth_nozice(Vtransform, Vstretching,theta_s, theta_b, hc, N,1, squeeze(h(:,:,mini_tt)),squeeze(iceDraft(:,:,mini_tt)),0); 
Zw(mini_tt,:,:,:)=set_depth_nozice(Vtransform, Vstretching,theta_s, theta_b, hc, N,5, squeeze(h(:,:,mini_tt)),squeeze(iceDraft(:,:,mini_tt)),0);
end
rmpath('/ds/projects/iomp/matlab_scripts/ROMS_MATLAB/utility/')

dx = permute(dx,[3 1 2]);
dy = permute(dy,[3 1 2]);
dz = Zw(:,:,:,2:end)-Zw(:,:,:,1:end-1);
dxdydz=bsxfun(@times,bsxfun(@times,dz,dx),dy);
dxdydz_masked=bsxfun(@times,dxdydz,permute(mask_rho,[3 1 2]));
dx = permute(dx,[2 3 1]);
dy = permute(dy,[2 3 1]);
dz = permute(dz,[2 3 4 1]);
Zw = permute(Zw,[2 3 4 1]);
Z = permute(Z,[2 3 4 1]);
dxdydz = permute(dxdydz,[2 3 4 1]);

mask_rho_3d = permute(repmat(permute(mask_rho,[3 1 2]),[1 1 1 N]),[2 3 4 1]);
% load data from output file
% flipud because the lat dirn in the model output is pointing in the opposite direction to that expected by ISOMIP+
m = flipud(double(ncread(hisname,'m')));
zeta = flipud(ncread(hisname,'zeta'));
temp = flipud(ncread(hisname,'temp'));
salt = flipud(ncread(hisname,'salt'));
u = flipud(ncread(hisname,'u'));
v = flipud(ncread(hisname,'v'));
Tb = flipud(double(ncread(hisname,'Tb')));
Sb = flipud(double(ncread(hisname,'Sb')));
ubar=flipud(ncread(hisname,'ubar'));
vbar=flipud(ncread(hisname,'vbar'));
% make requested data
meanMeltRate = squeeze(nansum(nansum(bsxfun(@times,m,Area_ice),2),1)) ./ squeeze(nansum(nansum(Area_ice,2),1));
meanMeltFlux = squeeze(nansum(nansum(bsxfun(@times,m,Area_ice*rho_i),2),1));
sumOfVolume = squeeze(nansum(nansum(bsxfun(@rdivide,bsxfun(@plus,h,zeta),(pm.*pn)),2),1));
meanTemperature = squeeze(nansum(nansum(nansum(bsxfun(@times,temp,bsxfun(@times,dxdydz,mask_rho_3d)),3),2),1)./nansum(nansum(nansum(bsxfun(@times,dxdydz,mask_rho_3d),3),2),1));
meanSalinity = squeeze(nansum(nansum(nansum(bsxfun(@times,salt,bsxfun(@times,dxdydz,mask_rho_3d)),3),2),1)./nansum(nansum(nansum(bsxfun(@times,dxdydz,mask_rho_3d),3),2),1));
meltRate=bsxfun(@times,m,mask_zice);
%u*
u_mod = (u(1:end-1,2:end-1,:,:)+u(2:end,2:end-1,:,:))/2;
v_mod = (v(2:end-1,1:end-1,:,:)+v(2:end-1,2:end,:,:))/2;
u_mod(end+1,:,:,:)=NaN; u_mod(:,end+1,:,:)=NaN; 
v_mod(end+1,:,:,:)=NaN;v_mod(:,end+1,:,:)=NaN;
u_mod(2:end+1,:,:,:)=u_mod;u_mod(1,:,:,:)=NaN;u_mod(:,2:end+1,:,:)=u_mod;u_mod(:,1,:,:)=NaN;
v_mod(2:end+1,:,:,:)=v_mod;v_mod(1,:,:,:)=NaN;v_mod(:,2:end+1,:,:)=v_mod;v_mod(:,1,:,:)=NaN;
Ustar = sqrt(Cd)*sqrt(sqrt(squeeze(u_mod(:,:,N,:)).^2 + squeeze(v_mod(:,:,N,:)).^2).^2+u_t^2); disp('ustar computed post-results from sqrt(cd)*sqrt(u^2+u_t^2)') 
Ustar = bsxfun(@times,Ustar,mask_zice);
%Tb is in situ -> convert to pot
addpath(genpath('/ds/projects/iomp/matlab_scripts/GSW'))
p_top = gsw_p_from_z(squeeze(Zw(:,:,N,:)),repmat(lat_rho,[1 1 length(ocean_time)]));
disp('Changing basal temperature to potential temperature')
clear TbPot
for tt = 1:size(temp,4)
TbPot(:,:,tt) = gsw_pt0_from_t(gsw_SA_from_SP(squeeze(Sb(:,:,tt)),p_top(:,:,tt),lon_rho,lat_rho),squeeze(Tb(:,:,tt)),p_top(:,:,tt));
if ~rem(tt,12)
disp([num2str(tt/size(temp,4)*100) '%'])
end
end
rmpath(genpath('/ds/projects/iomp/matlab_scripts/GSW'))
Tstar = bsxfun(@times,squeeze(temp(:,:,N,:))-TbPot,mask_zice); 
Sstar = bsxfun(@times,squeeze(salt(:,:,N,:))-Sb,mask_zice);
u_top = squeeze(u_mod(:,:,N,:));
v_top = squeeze(v_mod(:,:,N,:));
%psi = -cumsum(ubarMod2.*(h_Mod+zice_Mod).*dy(2:end-1,2:end-1),2);
ubar_mod = 1/2*(ubar(1:end-1,2:end-1,:)+ubar(2:end,2:end-1,:)); ubar_mod(isnan(ubar_mod))=0;
vbar_mod = 1/2*(vbar(1:end-1,2:end-1,:)+vbar(2:end,2:end-1,:)); vbar_mod(isnan(vbar_mod))=0;

baroSF = -cumsum(bsxfun(@times,ubar_mod,(h(2:end-1,2:end-1,:)+iceDraft(2:end-1,2:end-1,:)).*dy(2:end-1,2:end-1,:)),2);
baroSF(end+1,:,:)=NaN; baroSF(:,end+1,:)=NaN;
baroSF(2:end+1,:,:)=baroSF;baroSF(1,:,:)=NaN;baroSF(:,2:end+1,:)=baroSF;baroSF(:,1,:)=NaN;

if 0 %calculate OTSF
% dx*dz*v integrate into page, cumsum bottom-top
disp('interpolating v-velocity to z-levels')
[X_isurf,Y_isurf,Z_isurf]=meshgrid(squeeze(x_rho(1,:)),squeeze(y_rho(:,1)),Z_interp);
v_mod_i = nan(size(X_isurf,1),size(X_isurf,2),size(X_isurf,3),size(v_mod,4));
X_RHO = repmat(x_rho,[1 1 N]);
Y_RHO = repmat(y_rho,[1 1 N]);
for tt=1:size(v_mod,4)
v_mod_i(3:end-2,:,:,tt) = griddata(squeeze(X_RHO(3:end-2,:,:)),squeeze(Y_RHO(3:end-2,:,:)),squeeze(Z(3:end-2,:,:)),squeeze(v_mod(3:end-2,:,:,tt)),squeeze(X_isurf(3:end-2,:,:)),squeeze(Y_isurf(3:end-2,:,:)),squeeze(Z_isurf(3:end-2,:,:)));
v_mod_i(2,:,:,tt) = griddata(squeeze(X_RHO(2,:,:)),squeeze(Z(2,:,:)),squeeze(v_mod(2,:,:,tt)),squeeze(X_isurf(end-1,:,:)),squeeze(Z_isurf(end-1,:,:)));
v_mod_i(end-1,:,:,tt) = griddata(squeeze(X_RHO(end-1,:,:)),squeeze(Z(end-1,:,:)),squeeze(v_mod(end-1,:,:,tt)),squeeze(X_isurf(end-1,:,:)),squeeze(Z_isurf(end-1,:,:)));
if ~rem(tt,12)
disp([num2str(tt/size(v_mod,4)*100) '%'])
end
end

disp('masking interpolation within convex hull, but not within data')
%remove bad interp within convex hull
IN.maskICE=nan(size(v_mod_i,1),size(v_mod_i,2),size(v_mod_i,3));
disp('making ice mask for otsf')
for jj=1:size(v_mod,1)
IN.maskICE(jj,:,:) = double(inpolygon(squeeze(X_isurf(jj,:,:)),squeeze(Z_isurf(jj,:,:)),...
[squeeze(X_RHO(jj,:,N)),squeeze(X_RHO(jj,end,:))',fliplr(squeeze(X_RHO(jj,:,1))),fliplr(squeeze(X_RHO(jj,1,:))')],...
[squeeze(Z(jj,:,N)),squeeze(Z(jj,end,:))',fliplr(squeeze(Z(jj,:,1))),fliplr(squeeze(Z(jj,1,:))')]));
if ~rem(jj,5)
disp([num2str(jj/size(v_mod,1)*100) '%'])
end
end
IN.maskICE(IN.maskICE==0)=NaN;
v_mod_I = bsxfun(@times,v_mod_i,IN.maskICE); v_mod_I(isnan(v_mod_I))=0;
vbarx2= squeeze(nansum(bsxfun(@times,bsxfun(@times,v_mod_I,dx),(Z_interp(2)-Z_interp(1))*ones(size(dx))),1));
otSF = cumsum(vbarx2,2);
else
otSF = zeros(size(x_rho,2),length(Z_interp),length(ocean_time));
end

T_bot = squeeze(temp(:,:,1,:));
S_bot = squeeze(salt(:,:,1,:));

%interpolation of the form:
X_trans1=repmat(squeeze(x_rho(20,:))',[1 N length(ocean_time)]);
Z_trans1=squeeze(Z(20,:,:,:));
X_trans2=repmat(squeeze(x_rho(21,:))',[1 N length(ocean_time)]);
Z_trans2=squeeze(Z(21,:,:,:));
Y_trans1=repmat(squeeze(y_rho(:,100)),[1 N length(ocean_time)]);
Z_trans3=squeeze(Z(:,100,:,:));
Y_trans2=repmat(squeeze(y_rho(:,101)),[1 N length(ocean_time)]);
Z_trans4=squeeze(Z(:,101,:,:));
[X_isurf,Z_isurf]=meshgrid(squeeze(x_rho(20,:)),Z_interp);
tempXZ1=nan(size(X_isurf,1),size(X_isurf,2),size(temp,4));
tempXZ2=nan(size(X_isurf,1),size(X_isurf,2),size(temp,4));
saltXZ1=nan(size(X_isurf,1),size(X_isurf,2),size(temp,4));
saltXZ2=nan(size(X_isurf,1),size(X_isurf,2),size(temp,4));
[Y_isurf,Z2_isurf]=meshgrid(squeeze(y_rho(:,100)),Z_interp);
tempYZ1=nan(size(Y_isurf,1),size(Y_isurf,2),size(temp,4));
tempYZ2=nan(size(Y_isurf,1),size(Y_isurf,2),size(temp,4));
saltYZ1=nan(size(Y_isurf,1),size(Y_isurf,2),size(temp,4));
saltYZ2=nan(size(Y_isurf,1),size(Y_isurf,2),size(temp,4));
disp('interpolating transect data')
for tt=1:size(temp,4)
tempXZ1(:,:,tt) = griddata(squeeze(X_trans1(:,:,tt)),squeeze(Z_trans1(:,:,tt)),squeeze(temp(20,:,:,tt)),X_isurf,Z_isurf);
tempXZ2(:,:,tt) = griddata(squeeze(X_trans2(:,:,tt)),squeeze(Z_trans2(:,:,tt)),squeeze(temp(21,:,:,tt)),X_isurf,Z_isurf);
saltXZ1(:,:,tt) = griddata(squeeze(X_trans1(:,:,tt)),squeeze(Z_trans1(:,:,tt)),squeeze(salt(20,:,:,tt)),X_isurf,Z_isurf);
saltXZ2(:,:,tt) = griddata(squeeze(X_trans2(:,:,tt)),squeeze(Z_trans2(:,:,tt)),squeeze(salt(21,:,:,tt)),X_isurf,Z_isurf);
tempYZ1(:,:,tt) = griddata(squeeze(Y_trans1(:,:,tt)),squeeze(Z_trans3(:,:,tt)),squeeze(temp(:,100,:,tt)),Y_isurf,Z2_isurf);
tempYZ2(:,:,tt) = griddata(squeeze(Y_trans2(:,:,tt)),squeeze(Z_trans4(:,:,tt)),squeeze(temp(:,101,:,tt)),Y_isurf,Z2_isurf);
saltYZ1(:,:,tt) = griddata(squeeze(Y_trans1(:,:,tt)),squeeze(Z_trans3(:,:,tt)),squeeze(salt(:,100,:,tt)),Y_isurf,Z2_isurf);
saltYZ2(:,:,tt) = griddata(squeeze(Y_trans2(:,:,tt)),squeeze(Z_trans4(:,:,tt)),squeeze(salt(:,101,:,tt)),Y_isurf,Z2_isurf);
if ~rem(tt,12)
disp([num2str(tt/size(temp,4)*100) '%'])
end
end

%remove interp within convex hull
IN.maskXZ1 = double(inpolygon(X_isurf,Z_isurf,[X_trans1(:,N);X_trans1(end,:)';flipud(X_trans1(:,1));flipud(X_trans1(1,:)')],[Z_trans1(:,N);Z_trans1(end,:)';flipud(Z_trans1(:,1));flipud(Z_trans1(1,:)')])); IN.maskXZ1(IN.maskXZ1==0)=NaN;
IN.maskXZ2 = double(inpolygon(X_isurf,Z_isurf,[X_trans2(:,N);X_trans2(end,:)';flipud(X_trans2(:,1));flipud(X_trans2(1,:)')],[Z_trans2(:,N);Z_trans2(end,:)';flipud(Z_trans2(:,1));flipud(Z_trans2(1,:)')])); IN.maskXZ2(IN.maskXZ2==0)=NaN;
IN.maskYZ1 = double(inpolygon(Y_isurf,Z2_isurf,[Y_trans1(:,N);Y_trans1(end,:)';flipud(Y_trans1(:,1));flipud(Y_trans1(1,:)')],[Z_trans3(:,N);Z_trans3(end,:)';flipud(Z_trans3(:,1));flipud(Z_trans3(1,:)')])); IN.maskYZ1(IN.maskYZ1==0)=NaN;
IN.maskYZ2 = double(inpolygon(Y_isurf,Z2_isurf,[Y_trans2(:,N);Y_trans2(end,:)';flipud(Y_trans2(:,1));flipud(Y_trans2(1,:)')],[Z_trans4(:,N);Z_trans4(end,:)';flipud(Z_trans4(:,1));flipud(Z_trans4(1,:)')])); IN.maskYZ2(IN.maskYZ2==0)=NaN;

tempXZ=( bsxfun(@times,tempXZ1,IN.maskXZ1) + bsxfun(@times,tempXZ2,IN.maskXZ2) )/2;
saltXZ=( bsxfun(@times,saltXZ1,IN.maskXZ1) + bsxfun(@times,saltXZ2,IN.maskXZ2) )/2;
tempYZ=( bsxfun(@times,tempYZ1,IN.maskYZ1) + bsxfun(@times,tempYZ2,IN.maskYZ2) )/2;
saltYZ=( bsxfun(@times,saltYZ1,IN.maskYZ1) + bsxfun(@times,saltYZ2,IN.maskYZ2) )/2;

% fill BIG variables
fillstart = fillend+1; %needs to be ini at 0
fillend = fillstart+length(ocean_time)-1; %also needs to be ini at 0.

ocean_time_big(fillstart:fillend) = ocean_time;
meanMeltRate_big(fillstart:fillend)=meanMeltRate;
meanMeltFlux_big(fillstart:fillend)=meanMeltFlux;
sumOfVolume_big(fillstart:fillend)=sumOfVolume;
meanTemperature_big(fillstart:fillend)=meanTemperature;
meanSalinity_big(fillstart:fillend)=meanSalinity;
iceDraft_big(:,:,fillstart:fillend)=iceDraft;
bathymetry_big(:,:,fillstart:fillend)=h;
meltRate_big(:,:,fillstart:fillend)=meltRate;
Ustar_big(:,:,fillstart:fillend)=Ustar;
Tstar_big(:,:,fillstart:fillend)=Tstar;
Sstar_big(:,:,fillstart:fillend)=Sstar;
v_top_big(:,:,fillstart:fillend)=v_top;
u_top_big(:,:,fillstart:fillend)=u_top;
baroSF_big(:,:,fillstart:fillend)=baroSF;
otSF_big(:,:,fillstart:fillend)=otSF;
T_bot_big(:,:,fillstart:fillend)=T_bot;
S_bot_big(:,:,fillstart:fillend)=S_bot;
tempXZ_big(:,:,fillstart:fillend)=tempXZ;
saltXZ_big(:,:,fillstart:fillend)=saltXZ;
tempYZ_big(:,:,fillstart:fillend)=tempYZ;
saltYZ_big(:,:,fillstart:fillend)=saltYZ;

%fillstart=fillend+1
%if yy+1<10,yy_suffix = ['000',num2str(output_vec(yy+1))];
%elseif yy+1<100 & yy+1>9,yy_suffix = ['00',num2str(output_vec(yy+1))];
%elseif yy+1<1000,yy_suffix = ['0',num2str(output_vec(yy+1))];
%end
%fillend=fillstart+length(ncread([hisname_format,num2str(yy_suffix),'.nc'],'ocean_time')) - 1;



end  %stop looping through files



disp(' ')
disp([' Creating the file : ',outname])
disp(' ')

id = netcdf.create(outname, 'clobber');

% define dims
nx_dim = netcdf.defDim(id, 'nx', Mp);
ny_dim = netcdf.defDim(id, 'ny', Lp);
nz_dim = netcdf.defDim(id, 'nz', Zp);
nTime_dim=netcdf.defDim(id, 'nTime', time_length);

%define vars
x_id = netcdf.defVar(id, 'x', 'double', nx_dim);
netcdf.putAtt(id, x_id, 'long_name', 'cell centre in x-direction');
netcdf.putAtt(id, x_id, 'units', 'meter');
y_id = netcdf.defVar(id, 'y', 'double', ny_dim);
netcdf.putAtt(id, y_id, 'long_name', 'cell centre in y-direction');
netcdf.putAtt(id, y_id, 'units', 'meter');
z_id = netcdf.defVar(id, 'z', 'double', nz_dim);
netcdf.putAtt(id, z_id, 'long_name', 'cell centre in z-direction');
netcdf.putAtt(id, z_id, 'units', 'meter');
time_id = netcdf.defVar(id, 'time', 'double', nTime_dim);
netcdf.putAtt(id, time_id, 'long_name', 'time from the start of the simulation');
netcdf.putAtt(id, time_id, 'units', 'second');


meanMeltRate_id = netcdf.defVar(id, 'meanMeltRate', 'double', nTime_dim);
netcdf.putAtt(id, meanMeltRate_id, 'long_name', 'water equivalent melt rate, positive for melting and negative for freezing, averaged over the ice-shelf base');
netcdf.putAtt(id, meanMeltRate_id, 'units', 'm second^-1');
disp('!!!! check for water equiv')

totalMeltFlux_id = netcdf.defVar(id, 'totalMeltFlux', 'double', nTime_dim);
netcdf.putAtt(id, totalMeltFlux_id, 'long_name', 'the total mass flux of freshwater across the ice-ocean interface, positive for melting and negative for freezing');
netcdf.putAtt(id, totalMeltFlux_id, 'units', 'kg second^-1');


totalOceanVolume_id = netcdf.defVar(id, 'totalOceanVolume', 'double', nTime_dim);
netcdf.putAtt(id, totalOceanVolume_id, 'long_name', 'the total volume of the ocean');
netcdf.putAtt(id, totalOceanVolume_id, 'units', 'm^3');

meanTemperature_id = netcdf.defVar(id, 'meanTemperature', 'double', nTime_dim);
netcdf.putAtt(id, meanTemperature_id, 'long_name', 'the potential temperature averaged over the ocean volume');
netcdf.putAtt(id, meanTemperature_id, 'units', 'Celsius');

meanSalinity_id = netcdf.defVar(id, 'meanSalinity', 'double', nTime_dim);
netcdf.putAtt(id, meanSalinity_id, 'long_name', 'the salinity averaged over the ocean volume');
netcdf.putAtt(id, meanSalinity_id, 'units', 'PSU');

iceDraft_id = netcdf.defVar(id, 'iceDraft', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, iceDraft_id, 'long_name', 'the elevation of the ice-ocean interface');
netcdf.putAtt(id, iceDraft_id, 'units', 'm');
netcdf.putAtt(id, iceDraft_id, 'missing_value', -1e34);
netcdf.putAtt(id, iceDraft_id, '_FillValue', -1e34);

bathymetry_id = netcdf.defVar(id, 'bathymetry', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, bathymetry_id, 'long_name', 'the elevation of the bathymetry');
netcdf.putAtt(id, bathymetry_id, 'units', 'm');
netcdf.putAtt(id, bathymetry_id, 'missing_value', -1e34);
netcdf.putAtt(id, bathymetry_id, '_FillValue', -1e34);

meltRate_id = netcdf.defVar(id, 'meltRate', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, meltRate_id, 'long_name', 'the melt rate, positive for melting and negative for freezing');
netcdf.putAtt(id, meltRate_id, 'units', 'm s^-1');
netcdf.putAtt(id, meltRate_id, 'missing_value', -1e34);
netcdf.putAtt(id, meltRate_id, '_FillValue', -1e34);

frictionVelocity_id = netcdf.defVar(id, 'frictionVelocity', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, frictionVelocity_id, 'long_name', 'the friction velocity, used in melt calculations');
netcdf.putAtt(id, frictionVelocity_id, 'units', 'm s^-1');
netcdf.putAtt(id, frictionVelocity_id, 'missing_value', -1e34);
netcdf.putAtt(id, frictionVelocity_id, '_FillValue', -1e34);

thermalDriving_id = netcdf.defVar(id, 'thermalDriving', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, thermalDriving_id, 'long_name', 'the thermal driving used in melt calculations, calculated as the difference between theta in the boundary layer and freezing theta at the interface.');
netcdf.putAtt(id, thermalDriving_id, 'units', 'Celsius');
netcdf.putAtt(id, thermalDriving_id, 'missing_value', -1e34);
netcdf.putAtt(id, thermalDriving_id, '_FillValue', -1e34);

halineDriving_id = netcdf.defVar(id, 'halineDriving', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, halineDriving_id, 'long_name', 'the haline driving used in melt calculations, calculated as the difference between salinity in the boundary layer and salinity at the interface.');
netcdf.putAtt(id, halineDriving_id, 'units', 'PSU');
netcdf.putAtt(id, halineDriving_id, 'missing_value', -1e34);
netcdf.putAtt(id, halineDriving_id, '_FillValue', -1e34);


uBoundaryLayer_id = netcdf.defVar(id, 'uBoundaryLayer', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, uBoundaryLayer_id, 'long_name', 'the u-components of velocity in the boundary layer that were used to compute friction velocity');
netcdf.putAtt(id, uBoundaryLayer_id, 'units', 'm s^-1');
netcdf.putAtt(id, uBoundaryLayer_id, 'missing_value', -1e34);
netcdf.putAtt(id, uBoundaryLayer_id, '_FillValue', -1e34);

vBoundaryLayer_id = netcdf.defVar(id, 'vBoundaryLayer', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, vBoundaryLayer_id, 'long_name', 'the v-components of velocity in the boundary layer that were used to compute friction velocity');
netcdf.putAtt(id, vBoundaryLayer_id, 'units', 'm s^-1');
netcdf.putAtt(id, vBoundaryLayer_id, 'missing_value', -1e34);
netcdf.putAtt(id, vBoundaryLayer_id, '_FillValue', -1e34);

barotropicStreamfunction_id = netcdf.defVar(id, 'barotropicStreamfunction', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, barotropicStreamfunction_id, 'long_name', 'the barotropic streamfunction');
netcdf.putAtt(id, barotropicStreamfunction_id, 'units', 'm^3 s^-1');
netcdf.putAtt(id, barotropicStreamfunction_id, 'missing_value', -1e34);
netcdf.putAtt(id, barotropicStreamfunction_id, '_FillValue', -1e34);

overturningStreamfunction_id = netcdf.defVar(id, 'overturningStreamfunction', 'double', [nx_dim nz_dim nTime_dim]); 
netcdf.putAtt(id, overturningStreamfunction_id, 'long_name', 'the overturning streamfunction');
netcdf.putAtt(id, overturningStreamfunction_id, 'units', 'm^3 s^-1');
netcdf.putAtt(id, overturningStreamfunction_id, 'missing_value', -1e34);
netcdf.putAtt(id, overturningStreamfunction_id, '_FillValue', -1e34);

bottomTemperature_id = netcdf.defVar(id, 'bottomTemperature', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, bottomTemperature_id, 'long_name', 'the potential temperature in the bottom-most cell in each ocean column');
netcdf.putAtt(id, bottomTemperature_id, 'units', 'Celsius');
netcdf.putAtt(id, bottomTemperature_id, 'missing_value', -1e34);
netcdf.putAtt(id, bottomTemperature_id, '_FillValue', -1e34);

bottomSalinity_id = netcdf.defVar(id, 'bottomSalinity', 'double', [nx_dim ny_dim nTime_dim]);
netcdf.putAtt(id, bottomSalinity_id, 'long_name', 'the salinity in the bottom-most cell in each ocean column');
netcdf.putAtt(id, bottomSalinity_id, 'units', 'PSU');
netcdf.putAtt(id, bottomSalinity_id, 'missing_value', -1e34);
netcdf.putAtt(id, bottomSalinity_id, '_FillValue', -1e34);

temperatureXZ_id = netcdf.defVar(id, 'temperatureXZ', 'double', [nx_dim nz_dim nTime_dim]);
netcdf.putAtt(id, temperatureXZ_id, 'long_name', 'potential temperature transect in the x-z plane through the centre of the domain');
netcdf.putAtt(id, temperatureXZ_id, 'units', 'Celsius');
netcdf.putAtt(id, temperatureXZ_id, 'missing_value', -1e34);
netcdf.putAtt(id, temperatureXZ_id, '_FillValue', -1e34);

salinityXZ_id = netcdf.defVar(id, 'salinityXZ', 'double', [nx_dim nz_dim nTime_dim]);
netcdf.putAtt(id, salinityXZ_id, 'long_name', 'salinity transect in the x-z plane through the centre of the domain, y=40km');
netcdf.putAtt(id, salinityXZ_id, 'units', 'PSU');
netcdf.putAtt(id, salinityXZ_id, 'missing_value', -1e34);
netcdf.putAtt(id, salinityXZ_id, '_FillValue', -1e34);

temperatureYZ_id = netcdf.defVar(id, 'temperatureYZ', 'double', [ny_dim nz_dim nTime_dim]);
netcdf.putAtt(id, temperatureYZ_id, 'long_name', 'potential temperature transect in the y-z plane through the centre of the domain, x=520km');
netcdf.putAtt(id, temperatureYZ_id, 'units', 'Celsius');
netcdf.putAtt(id, temperatureYZ_id, 'missing_value', -1e34);
netcdf.putAtt(id, temperatureYZ_id, '_FillValue', -1e34);

salinityYZ_id = netcdf.defVar(id, 'salinityYZ', 'double', [ny_dim nz_dim nTime_dim]);
netcdf.putAtt(id, salinityYZ_id, 'long_name', 'salinity transect in the y-z plane through the centre of the domain, x=520km');
netcdf.putAtt(id, salinityYZ_id, 'units', 'PSU');
netcdf.putAtt(id, salinityYZ_id, 'missing_value', -1e34);
netcdf.putAtt(id, salinityYZ_id, '_FillValue', -1e34);

%add global data
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', 'ISOMIP results: Regional Ocean Modelling System');
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'date', date);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'clim_file', outname);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'grd_file', grdname);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'type', 'RESULTS file');
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', 'ROMS');
netcdf.endDef(id);

%fill in variables
netcdf.putVar(id, x_id, x_rho(1,:));
netcdf.putVar(id, y_id, y_rho(:,1));
netcdf.putVar(id, z_id, flip(Z_interp)); %flip dim to start at -2.5 and go through to -717.5
netcdf.putVar(id, time_id, ocean_time_big);
netcdf.putVar(id, meanMeltRate_id, meanMeltRate_big);
netcdf.putVar(id, totalMeltFlux_id, meanMeltFlux_big);
netcdf.putVar(id, totalOceanVolume_id, sumOfVolume_big);
netcdf.putVar(id, meanTemperature_id, meanTemperature_big);
netcdf.putVar(id, meanSalinity_id, meanSalinity_big);
netcdf.putVar(id, iceDraft_id, permute(iceDraft_big,[2 1 3]));
netcdf.putVar(id, bathymetry_id, permute(bathymetry_big,[2 1 3]));
netcdf.putVar(id, meltRate_id, permute(meltRate_big,[2 1 3]));
netcdf.putVar(id, frictionVelocity_id, permute(Ustar_big,[2 1 3]));
netcdf.putVar(id, thermalDriving_id, permute(Tstar_big,[2 1 3]));
netcdf.putVar(id, halineDriving_id, permute(Sstar_big,[2 1 3]));
netcdf.putVar(id, uBoundaryLayer_id, permute(v_top_big,[2 1 3])); %note that v and u are swapped here. this is because in this setup the left-right direction (=u velocity) is actually rotated CW by 90deg, and so points up-down in the requested output. So, the simulatd u is actually the output v. (and vice versa). Also note that as u-vel from ROMS is positive to the right (by convention), AND we switch the u and v velocities, that means that the previous positive u is pointing downwards. However, the requested output convention is positive upwards, meaning that the model u needs to become the output v and the sign flipped (*-1). For example, flow that moves into the cavity on the right hand side, moves clockwise and exits, would be negative u flow in the model raw output, but in the isomip+ requested output, it will be positive v flow.
netcdf.putVar(id, vBoundaryLayer_id, permute(-u_top_big,[2 1 3]));
netcdf.putVar(id, barotropicStreamfunction_id, permute(baroSF_big,[2 1 3]));
netcdf.putVar(id, overturningStreamfunction_id, otSF_big);
netcdf.putVar(id, bottomTemperature_id, permute(T_bot_big,[2 1 3]));
netcdf.putVar(id, bottomSalinity_id, permute(S_bot_big,[2 1 3]));
netcdf.putVar(id, temperatureXZ_id, flip(permute(tempXZ_big,[2 1 3]),2));
netcdf.putVar(id, salinityXZ_id, flip(permute(saltXZ_big,[2 1 3]),2));
netcdf.putVar(id, temperatureYZ_id, flip(permute(tempYZ_big,[2 1 3]),2));
netcdf.putVar(id, salinityYZ_id, flip(permute(saltYZ_big,[2 1 3]),2));

netcdf.close(id);


