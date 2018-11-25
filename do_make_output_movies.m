filename='/home/ubuntu/IceOceanVolume/ISOMIP_PLUS/Ocean3/MakeOSF/Ocean3_COM_ROMSUTAS_FISOC_OSF.nc'
PlotName='Ocean3_movie'
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

spy=(365*24*60*60);

%%
disp('start plots')

v = VideoWriter(PlotName);
open(v)
fig=figure('visible','off','position',[700 1000 900 740]),
for tt=1:length(time)

set(0,'defaultAxesFontSize',10)
smplot(3,2,1,'axis','on')
plot(time/(60*60*24*365),meanMeltRate*(60*60*24*365))
hold on,plot(time(tt)/(60*60*24*365),meanMeltRate(tt)*(60*60*24*365),'.')
ntitle('melt rate','fontsize',10)
smplot(3,2,2)
flat_pcolor(X(:,:,1),Y(:,:,1),meltRate(:,:,tt)*(60*60*24*365)); colorbar('location','east'),caxis([-10 70]),cmocean('curl','pivot',0)
ntitle('melt rate','fontsize',10)
set(gca, 'XTickLabel', []),set(gca,'YTickLabel',[])
smplot(3,2,3)
flat_pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),temperatureXZ(:,:,tt)); caxis([-2.3 1]),cmocean('thermal'),colorbar('location','west')
ntitle('XZ transect temp','fontsize',10)
set(gca, 'XTickLabel', []),set(gca,'YTickLabel',[])
smplot(3,2,4)
flat_pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),salinityXZ(:,:,tt)); caxis([33.6 34.5]),cmocean('haline'),colorbar('location','west')
ntitle('XZ transect salt','fontsize',10)
set(gca, 'XTickLabel', []),set(gca,'YTickLabel',[])
smplot(3,2,5)
flat_pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),overturningStreamfunction(:,:,tt)); cmocean('balance','pivot',0),colorbar('location','west');
ntitle('overturning streamfunction','fontsize',10),caxis([-3e5 3e5]),cmocean('balance','pivot',0)
set(gca, 'XTickLabel', []),set(gca,'YTickLabel',[])
smplot(3,2,6)
flat_pcolor(X(:,:,1),Y(:,:,1),barotropicStreamfunction(:,:,tt).*mask_ocean(:,:,tt)); cmocean('balance','pivot',0),colorbar('location','west');
ntitle('barotropic streamfunction','fontsize',10),caxis([-5e5 1e6]),cmocean('balance','pivot',0)
set(gca, 'XTickLabel', []),set(gca,'YTickLabel',[])
set(gca,'nextplot','replacechildren');

   frame = getframe(fig);
   writeVideo(v,frame);
disp(['done frame ' num2str(tt)])
clf
end
close(v)

