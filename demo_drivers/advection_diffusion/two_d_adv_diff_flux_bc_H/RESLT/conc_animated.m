clearvars;close all;
RECORD=1;

fstruct = dir('soln*.dat');
data_length=length(fstruct)-1;

xxr=0:0.01:2;
yyr=0:0.1:1200;

f=figure;
set(gcf,'Position',[680 103 235 775]);

if RECORD
    videohold = VideoWriter('myFile');
    set(videohold,'FrameRate',1)
    open(videohold)
    af=struct([]);
end
kk=0;
for ii=0:data_length
kk=kk+1;
filename=strcat(['soln' num2str(ii) '.dat']);
newfile=strcat(['sol' num2str(ii) 'p.dat']);
readdata=fileread(filename);
readdata=regexprep(readdata,'ZONE I=5, J=5','');
filewrite=fopen(newfile,'w');
fwrite(filewrite,readdata);
fclose(filewrite);

clear readdata filename;

aa=load(newfile);
rangevec=1:length(aa);

[xx,yy]=meshgrid(xxr,yyr);
zz=griddata(aa(rangevec,1),aa(rangevec,2),aa(rangevec,3),xx,yy);

clf;
h=pcolor(xx,yy,zz);
set(h, 'EdgeColor', 'none');
colorbar;
set(gca, 'YDir','reverse')

aframe=getframe(f);
if kk==1 % jj is loop counter used for the animation
   af=aframe;
else
   af(kk)=aframe;
end

clear aa;
delete(newfile);

end

writeVideo(videohold,af) %saving
close(videohold) %closing
clearvars;