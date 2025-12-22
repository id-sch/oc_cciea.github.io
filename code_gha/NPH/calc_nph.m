% Use this to read the netcdf


%-----------------------------------------------
%--Load Packages
%-----------------------------------------------
pkg load image
pkg load netcdf
%-----------------------------------------------
%--Path to functions
%-----------------------------------------------
addpath ("./code_gha/NPH/")
%-----------------------------------------------


%--lon and lat of UI locations that are used to calculated the gradient
lat_ui = [33,36,39,42,45,48,51]';
lon_ui = [241,238,235,235,235,235,229]';
num_ui = i_size(lat_ui);

%--load the fnmoc data from erddap
% load slpM_1967_2018_5.5_65_150_250.mat

% read the netcdf files
% fn = '/home/isaac/data_files/Work/TS/data/python_erddap/slp_fnmoc_6hourly/pmsl_monthly.nc'
fn = './data_gha/NPH/pmsl_monthly_update.nc'

lat_slp = ncread(fn, 'lat_vec');
lon_slp = ncread(fn, 'lon_vec');
data = ncread(fn, 'pmsl');
time = ncread(fn, 'time');
time1 = datevec(time+datenum(1967,1,1));
yy1 = double(time1(:,1));
mm1 = double(time1(:,2));
dd1 = double(time1(:,3));
jd_slp = date2jd(yy1,mm1,dd1);

nt_all = size(yy1)(1);
yr_end = yy1(nt_all);
mon_end = mm1(nt_all);


% permute slp so that it is (lat x lon x time)
slp = permute(data, [2,1,3]);

#--get the year and month
[yy_slp,mm_slp] = jd2date(jd_slp); 


#--get time wanted
in_yr = find(jd_slp<=date2jd(yr_end,mon_end));

#--get size of domain
in_lat =find(lat_slp >=3 & lat_slp<=50);
in_lon = find(lon_slp >= 200 & lon_slp <= 250);

#--create final lon,lat,slp data
slp_final = slp(in_lat,in_lon,in_yr);
lon_slp = lon_slp(in_lon);
lat_slp = lat_slp(in_lat);
jd_slp = jd_slp(in_yr);
[yy_slp,mm_slp] = jd2date(jd_slp);


[ny,nx,nt] = size(slp_final);

%--need this to not get an error when contouring
[long,latg] = meshgrid(lon_slp,lat_slp);  
				
%--find the index (on the SLP grid) of the UI locations
ilat_ui = zeros(size(lat_ui));
ilon_ui = zeros(size(lon_ui));
for i=1:num_ui
   [mmm,ilat]  = min(abs(lat_ui(i)-lat_slp));
   [mmm,ilon]  = min(abs(lon_ui(i)-lon_slp));
   ilat_ui(i,1) = ilat;
   ilon_ui(i,1) = ilon;
end

%--Parameter: This is an important parameter, you need to set this to constrain the
%size of the area which the NPH is located, should not be to small, then
%it can conflict with very small areas of high pressure
d_area = 5000;  %--region has to be bigger than this to be
                %--considered the NPH

%--setup meshgrid of xp and yp, this will be used to calculate the
%center of mass and strenght like the old way, doing this as a
%comparioson
[xp,yp] = meshgrid(1:nx,1:ny);

%--Initialize nph location vectors and matrices
x_nph = zeros(nt,1);
y_nph = zeros(nt,1);
a_nph = zeros(nt,1);
s_nph = zeros(nt,1);
max_nph = zeros(nt,1);
num_nph = zeros(nt,1);
area_nph = zeros(nt,1);
x_max_nph = zeros(nt,1);
y_max_nph = zeros(nt,1);
grad_slp  = zeros(nt,num_ui);
dgrad_slp = zeros(nt,num_ui);
sgrad_slp = zeros(nt,num_ui);
dx_slp = zeros(nt,num_ui);

%--inititalize vectors of old location and strength
xbar = zeros(nt,1);
ybar = zeros(nt,1);
mbar = zeros(nt,1);
slp_lvl_vec = zeros(nt,1);

%--Parameter: Look at the SLP within this contour
slp_lvl_wnt = 1020;
slp_lvl_wnto = 1020;

nan_flag = 0;
for jj=1:nt
   disp([jj nt])
   slp1 = squeeze(slp_final(:,:,jj));

   % check to see if slp1 has data or is all NaN, see NOTE at top
   in_nan = isnan(slp1);
   [ny1, nx1] = size(in_nan);
   in_nan_re = reshape(in_nan, [ny1*nx1, 1]);
   [ndata1, ndata2] = size(find(in_nan_re==0));

   % if slp1 is all NaN, use the previous months data, see NOTE at top
   if(ndata1==0)
     slp1 = squeeze(slp_final(:,:,jj-1));
     indx0 = jj;
     nan_flag = 1;
   end

   slp_loop = 1;   
   while(slp_loop)   %--while there is a region of SLP as defined by slp_lvl_wnt
      inm1 = find(slp1 < slp_lvl_wnt);
      slp1(inm1) = 0;


      [mnrght,in_rght] = min(abs(lon_slp-235));
      [mntop,in_top]   = min(abs(lat_slp-45));
      slp1(:,in_rght:end) = 0;
      slp1(in_top:end,:)  = 0;

      islp1 = slp1 > 0;      

      [L,num] = bwlabel(islp1);  %--label each region that has slp_lvl_wnt
      strct_prop = regionprops(L);
      unque_L=unique(L); 
      num_rgns = i_size(unque_L)-1;

      %--find the NPH
      % if slp1 is all NaN, use the previous months data, see NOTE at top
      if(nan_flag)
          slp1 = squeeze(slp_final(:,:,jj-1));
      else
          slp1 = squeeze(slp_final(:,:,jj));
      end

      
      area_rgns = zeros(num_rgns,1);
      xpos_rgns = zeros(num_rgns,1);
      ypos_rgns = zeros(num_rgns,1);
      for i=1:num_rgns  %--get location and area of each region
	 L1=L==i;
	 slp1L = slp1.*L1;
	 slp2L = slp1L > 0;
	 strct_rgns = regionprops(slp2L); 
	 area_rgns(i) = strct_rgns.Area;
	 xpos_rgns(i) = strct_rgns.Centroid(1);
	 ypos_rgns(i) = strct_rgns.Centroid(2);
      end
      
      %--eliminate certain regions that are too far north and east
      %--might need to change (240,40)
      [mneast,in_east] = min(abs(lon_slp-240));
      [mnwest,in_west] = min(abs(lon_slp-200));
      [mnnorth,in_north] = min(abs(lat_slp-40));
      in_ysml = find(ypos_rgns < in_north);
      in_xsml = find(xpos_rgns(in_ysml) > in_west & xpos_rgns(in_ysml) < in_east);

      %--eliminate regions that are too small as defined by parameter d_area
      if (i_size(in_xsml) > 1)
	 in_asml = find(area_rgns(in_ysml(in_xsml))>d_area);
	 in_nph  = in_ysml(in_xsml(in_asml));
      else
	 in_nph  = in_ysml(in_xsml);
      end
 
      slp_lvl_wntf = slp_lvl_wnt;      
      if(in_nph)  %--exit the while loop if there is a region on nph
         slp_loop = 0;
         slp_lvl_wnt = slp_lvl_wnto; 
      else        %--continue the loop with a new slp_lvl_wnt if there
		  %--is no region of nph
         slp_loop = 1;
         slp_lvl_wnt = slp_lvl_wnt-0.25;
         slp_lvl_wntf = slp_lvl_wnt;
         fff(jj,:) = [jj,slp_lvl_wnt,slp_lvl_wntf,yy_slp(jj),mm_slp(jj)];
      end
   end
   slp_lvl_vec(jj,1) = slp_lvl_wntf;
   %--now that we have the region of the NPH, then we can get some
   %--information on it such as its (x,y) position and area
   x_nph(jj,1)    = xpos_rgns(in_nph);
   y_nph(jj,1)    = ypos_rgns(in_nph);
   a_nph(jj,1)    = area_rgns(in_nph);   

   %--now we want more information about the NPH, such as maximum SLP, 
   %--location of the maximum and gradients between the center and UI locations   
   L1=L==in_nph;  %--setup logical matrix that is the region of the NPH
   slp_nph = slp1.*L1;  %--just get slp in this region

   %--this finds the area of the SLP with in the contour, units are km^2
   area_all=fun_calc_area(slp_nph,long,latg);
   area_nph(jj,1) = area_all;

   %--this finds the mean of the SLP in the NPH region
   in_p = find(slp_nph > slp_lvl_wntf);
   sum_p = sum(slp_nph(in_p));
   
   mn_p  = sum_p./i_size(in_p);
   s_nph(jj,1)    = mn_p;
   num_nph(jj,1)  = i_size(in_p);
   %--this finds the (x,y) position of the maximum SLP in the NPH region
   max_nph1 = max(max(slp_nph));
   [jp,ip] = find(slp_nph == max_nph1);
   max_nph(jj) = max_nph1;
   x_max_nph(jj) = ip(1);
   y_max_nph(jj) = jp(1);


   %------------------calculate the center of mass and strength like the
   %old way
   sum_xm = 0;
   sum_ym = 0;
   sum_m  = 0;
   in_lrg  = in_p;
   num_lrg = i_size(in_lrg);
   for k=1:num_lrg
      x1  = xp(in_lrg(k));
      y1  = yp(in_lrg(k));
      m1  = slp_nph(in_lrg(k));
      sum_xm = sum_xm+slp_nph(in_lrg(k))*xp(in_lrg(k));
      sum_ym = sum_ym+slp_nph(in_lrg(k))*yp(in_lrg(k));     
      sum_m  = sum_m+slp_nph(in_lrg(k));
   end
   xbar(jj,1) = sum_xm/sum_m;
   ybar(jj,1) = sum_ym/sum_m;
   mbar(jj,1) = sum_m;

   %--this finds the gradients between the center of the NPH region and
   %--the UI locations -- might want to find the gradients between the
   %--maximum SLP and the UI locations -- then need to comment out (ip,jp)
   %--lines below
   ip = round(x_nph(jj));
   jp = round(y_nph(jj));
   for iii=1:num_ui   
      ip_ui = ilon_ui(iii);
      jp_ui = ilat_ui(iii);
      
      if (ip_ui == ip)
         ip=ip-1;
      end
      if (jp_ui == jp)
         jp=jp-1;
      end

      dx = (ip_ui-ip)/abs(ip_ui-ip);
      dy = (jp_ui-jp)/abs(jp_ui-jp);
      xm = ip:dx:ip_ui;
      ym = jp:dy:jp_ui;
      numx = i_size(xm);
      numy = i_size(ym);

      m = (jp_ui-jp)/(ip_ui-ip);
      b = jp_ui-m*ip_ui;

      if(numx<numy)
         yline = ym;
         xline = round((yline-b)/m);
      else
         xline = xm;
         yline = round(m*xline+b);
      end

      num_xline = i_size(xline);
      slp_line  = zeros(num_xline,1);
      for jjj=1:num_xline
         slp_line(jjj) = slp1(yline(jjj),xline(jjj));
      end

      %--grad as defined as (slp_nph-slp_ui)/dist_between
      dist_line = vdist(lat_slp(yline(1)),lon_slp(xline(1)),lat_slp(yline(end)),lon_slp(xline(end)))/1000;
      %grad_slp(jj,iii)  = mean(slp_line)/dist_line;
      grad_slp(jj,iii)  = (slp_line(1)-slp_line(end))/dist_line;
      dgrad_slp(jj,iii) = dist_line;
      sgrad_slp(jj,iii) = mean(slp_line);     

      %--grad as defined as 
      num_line = i_size(slp_line);
      mid_pnt  = num_line/2;
      bgn_pnt  = round(mid_pnt-num_line/5);
      end_pnt  = round(mid_pnt+num_line/5);
      dff_slp  = diff(slp_line);
      dx_slp(jj,iii)   = mean(dff_slp(bgn_pnt:end_pnt)); 
   end

   nan_flag = 0;
   
end

   
% convert xy indices to lon and lat
nx = i_size(lon_slp);
ny = i_size(lat_slp);
dx = 0.1;
dy = 0.1;
inx = 1:dx:nx;
iny = 1:dy:ny;
lati = interp1(1:ny, lat_slp, iny);
loni = interp1(1:nx, lon_slp, inx);

latf = zeros(nt, 1);
lonf = zeros(nt, 1);
for j=1:nt
	[min_lat, in_lat] = min(abs(iny-y_nph(j)));
        latf(j, 1) = lati(in_lat);
	[min_lon, in_lon] = min(abs(inx-x_nph(j)));
        lonf(j, 1) = loni(in_lon);
end

% save('-mat4-binary','mat_files/nph_1020_fmnoc_from_6hr_new.mat','x_nph','y_nph','a_nph','s_nph','max_nph','x_max_nph','y_max_nph','grad_slp','dgrad_slp','sgrad_slp','dx_slp','lon_slp','lat_slp','slp_lvl_wnt','d_area','jd_slp','xbar','ybar','mbar','num_nph','slp_lvl_vec','area_nph','fff')


% save area in csv
output1 = [yy_slp,mm_slp,area_nph,max_nph,x_nph,y_nph, lonf, latf];
csvwrite('year_mon_area_max_x_y_lon_lat.csv',output1)
