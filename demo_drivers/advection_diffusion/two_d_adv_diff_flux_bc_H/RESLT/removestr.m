clearvars;
file=['soln17.dat'];
newfile='soln00.dat';
readdata=fileread(file);
readdata=regexprep(readdata,'ZONE I=5, J=5','');

filewrite=fopen(newfile,'w');
fwrite(filewrite,readdata);
fclose(filewrite);

%%

% aa=load('soln.dat');
% 
% rangevec=floor(linspace(1,length(aa),500));
% [xx,yy]=meshgrid(aa(rangevec,1),aa(rangevec,2));
% zz=griddata(aa(rangevec,1),aa(rangevec,2),aa(rangevec,3),xx,yy);
% pcolor(xx,yy,zz)
% colorbar;

%%
clearvars;
aa=load('soln00.dat');

rangevec=1:length(aa);

xxr=0:0.01:2;
yyr=0:0.1:100;

[xx,yy]=meshgrid(xxr,yyr);
zz=griddata(aa(rangevec,1),aa(rangevec,2),aa(rangevec,3),xx,yy);

% subplot(122);
figure;
h=pcolor(xx,yy,zz);
set(h, 'EdgeColor', 'none');
colorbar;
set(gcf,'Position',[680 103 235 775])
set(gca, 'YDir','reverse')
%%
% 
% xr=0:0.1:2;
% 
% yr=(1-(xr-1).^2); plot(xr,yr)
%%
% r=0:0.01:1; z=0:0.01:10; [rr,zz]=meshgrid(r,z); % u=exp(-4*zz.*zz).*exp(
% -(rr-0.1*1.0).*(rr-0.1*1.0)./(0.01*0.01*1.0*1.0)); u=exp(-4*rr.*rr).*exp(
% -(zz-0.1*1.0).*(zz-0.1*1.0)./(0.01*0.01*1.0*1.0));
% 
% h=pcolor(r,z,u); set(h, 'EdgeColor', 'none');