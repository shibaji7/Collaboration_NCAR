lat_lon_fix1(:,1)=fix((lat_lon_sort(:,1)*1))/1;
lat_lon_fix1(:,2)=fix((lat_lon_sort(:,2)*1))/1;
lat_lon_fix1(:,3)=lat_lon_sort(:,3);
lat_lon_fix_unq1=unique(lat_lon_fix1,'rows');
latr_min= min(lat_lon_fix_unq1(:,1));
latr_max= max(lat_lon_fix_unq1(:,1));
lonr_min= min(lat_lon_fix_unq1(:,2));
lonr_max= max(lat_lon_fix_unq1(:,2));
i=1;
j=1;
z_mean = lat_lon_fix_unq1(1,3);
for latr_num = latr_min:1:latr_max
for lonr_num = lonr_min:1:lonr_max
ind = find(lat_lon_fix_unq1(:,1)==latr_num & lat_lon_fix_unq1(:,2)==lonr_num);

if length(ind) ~= 0
z_mean = mean(lat_lon_fix_unq1(ind,3));
z_std = std(lat_lon_fix_unq1(ind,3));
lat_lon_new(i,:) = [latr_num lonr_num z_mean];
lat_lon_new1(j,:) = [latr_num lonr_num z_mean];
i=i+1;
j=j+1;
else
lat_lon_new1(j,:) = [latr_num lonr_num z_mean];
j=j+1;
end
end
end