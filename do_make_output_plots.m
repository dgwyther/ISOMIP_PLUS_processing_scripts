%define output file
filename='/home/ubuntu/IceOceanVolume/ISOMIP_PLUS/Ocean3/MakeOSF/Ocean3_COM_ROMSUTAS_FISOC_OSF.nc'
%load data
x=ncread(filename,'x');
y=ncread(filename,'y');
z=ncread(filename,'z');
[X,Y,Z]=ndgrid(x,y,z);
time=ncread(filename,'time');
meanMeltRate=ncread(filename,'meanMeltRate');
totalMeltFlux=ncread(filename,'totalMeltFlux');
meltRate=ncread(filename,'meltRate');
barotropicStreamfunction=ncread(filename,'barotropicStreamfunction');
overturningStreamfunction=ncread(filename,'overturningStreamfunction');
temperatureXZ=ncread(filename,'temperatureXZ');
salinityXZ=ncread(filename,'salinityXZ');
temperatureYZ=ncread(filename,'temperatureYZ');
salinityYZ=ncread(filename,'salinityYZ');
uBoundaryLayer=ncread(filename,'uBoundaryLayer');
mask_ocean = double((uBoundaryLayer<1e30)); mask_ocean(mask_ocean==0)=NaN;
% plots
disp('doing plots')
% transect plots
figure('pos',[1233          63        1662        1695])
set(0,'defaultAxesFontSize',10)
smplot(4,2,1,'axis','on')
plot(time/(60*60*24*365),meanMeltRate*(60*60*24*365))
yyaxis right
plot(time/(60*60*24*365),totalMeltFlux*(1e-12*60*60*24*365))
ntitle('mean melt rate (left axis) and total mass loss (right)','fontsize',10)
smplot(4,2,2,'axis','on')
flat_pcolor(X(:,:,1),Y(:,:,1),meltRate(:,:,end)*(60*60*24*365)); colorbar('location','east'),cmocean('curl','pivot',0),caxis([-10 70])
ntitle('melt rate','fontsize',10)
smplot(4,2,3,'axis','on')
flat_pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),temperatureXZ(:,:,end)); caxis([-2.3 1]),cmocean('thermal')
ntitle('XZ transect temp','fontsize',10)
smplot(4,2,4,'axis','on')
flat_pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),salinityXZ(:,:,end)); caxis([33.6 34.5]),cmocean('haline')
ntitle('XZ transect salt','fontsize',10)
smplot(4,2,5,'axis','on')
flat_pcolor(squeeze(Y(1,:,:)),squeeze(Z(1,:,:)),temperatureYZ(:,:,end)); caxis([-2.3 1]),cmocean('thermal'),colorbar('location','south')
ntitle('YZ transect temp','fontsize',10)
smplot(4,2,6,'axis','on')
flat_pcolor(squeeze(Y(1,:,:)),squeeze(Z(1,:,:)),salinityYZ(:,:,end)); caxis([33.6 34.5]),cmocean('haline'),colorbar('location','south')
ntitle('YZ transect salt','fontsize',10)
smplot(4,2,7,'axis','on')
flat_pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),overturningStreamfunction(:,:,end)); cmocean('balance','pivot',0),colorbar('location','west');
ntitle('overturning streamfunction','fontsize',10)
smplot(4,2,8,'axis','on')
flat_pcolor(X(:,:,1),Y(:,:,1),barotropicStreamfunction(:,:,end).*mask_ocean(:,:,end)); cmocean('balance','pivot',0),colorbar('location','west');
ntitle('barotropic streamfunction','fontsize',10)
plotname=[filename,'.summaryplots.png']
export_fig(plotname,'-png','-transparent')
