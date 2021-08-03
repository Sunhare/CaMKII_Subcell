clear all;
XN=10;
YN=10;
if (exist('ci.dat.xz')==2)
system('xz --decompress -k ci.dat.xz');
end

c=colormap(hsv(2400));
c(1601:2400,:)=[];
c=flipud(c);
colormap(c);
clen=(length(c)-1);
fid = fopen('ci.dat','r');
fseek(fid,0,'eof');
t=ftell(fid)/1/XN/YN;
fclose(fid);

cimm=load('cimm.txt');
cimmmax=cimm(:,1);
cimmmin=cimm(:,2);
cimax=0.2;
cimin=0.1;





fidci = fopen('ci.dat','r');

scrsz = get(0,'ScreenSize');
fig=figure('Position',[200 200 300 300],'Color','k');

mov = VideoWriter('tmp.avi','Uncompressed AVI')
open(mov);
[x,y]=meshgrid(1:XN,1:YN);


for frm=1:t

if (mod(frm,10)==0)frm/t*100
end

[ciz,count] = fread(fidci,[YN,XN],'uchar');
ciz=ciz*(cimmmax(frm)-cimmmin(frm))/255+cimmmin(frm);
surf(x,y,ciz,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
shading flat
axis([1,XN,1,YN])
axis 'auto z'
daspect([100 100 1])
caxis([cimin cimax])
view(2)
camlight
axis off




F = getframe(fig);
writeVideo(mov,F);
end



fclose(fidci);

close(mov);


if (exist('ci.dat.xz')==2)
delete('ci.dat');
end



comm=['/usr/local/bin/ffmpeg -i tmp.avi -c:v libx264 -vf format=yuv420p -preset placebo -y ca.mp4'];
system(comm);

if (exist('ca.mp4')==2)
delete('tmp.avi');
end

close all;
