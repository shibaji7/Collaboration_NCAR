clc;
clear;
load('/home/shibaji7/PHARLAP_SAMI_INTG/pharlap_4.1.3/sami_parsed.mat');
ntm = 1;
dene= dene(:,:,ntm);
deneold = dene;

prompt = 'Center latitude of eclipse (deg):';
d_lat = input(prompt);
% d_lat = 0;

attn = 0.8/30*abs(lat-d_lat)+0.2;
size(attn)
for i=1:nz
    for j=1:nf
        if (attn(i,j))>1
            attn(i,j)=1;
        end
    end 
end
dene = dene.*attn;

figure(1)
subplot(1,2,1)
pcolor(lat,alt,deneold)
shading interp
title('no eclipse')
xlabel('latitude')
ylabel('altitude')
axis([-inf inf 0 1000])
colorbar

subplot(1,2,2)
pcolor(lat,alt,dene)
shading interp
title('with eclipse')
xlabel('latitude')
ylabel('altitude')
axis([-inf inf 0 1000])
colorbar

figure(2)
pcolor(lat,alt,attn)
shading interp
title('density attenuation')
xlabel('latitude')
ylabel('altitude')
%axis([-inf inf 0 1000])
colorbar