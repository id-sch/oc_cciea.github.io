%%% fun_calc_area.m --- 
%% 
%% Filename: fun_calc_area.m
%% Description: 
%% Author: isaac.schroeder
%% Maintainer: 
%% Created: Fri Jun  3 15:53:27 2011
%% Version: $Id$
%% Last-Updated: Fri Jun  3 15:55:36 2011
%%           By: isaac.schroeder
%%     Update #: 2
%% Keywords: 
%% Compatibility: 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Commentary: 
%% 
%% 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Change log:
%% 
%% RCS $Log$
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Code:




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function area_all=fun_calc_area(slp_nph,long,latg);

%load test_area.mat

ind = find(slp_nph>0);
slp_nph(ind) = 1;

lat = latg(:,1);
lon = long(1,:)';

nx = i_size(lon);
ny = i_size(lat);


area_mtrx = zeros(ny,nx);
for k=1:ny-1
   for j=1:nx-1
      pt1 = slp_nph(k,j);
      pt2 = slp_nph(k,j+1);
      pt3 = slp_nph(k+1,j);
      pt4 = slp_nph(k+1,j+1);
      sum_pts = pt1+pt2+pt3+pt4;
      if(sum_pts == 3)
          %         disp([k,j])
         vd_x = vdist(lat(k),lon(j),lat(k),lon(j+1))/1000;
         vd_y = vdist(lat(k),lon(j),lat(k+1),lon(j))/1000;
         area_rec = vd_x*vd_y/2;  
      else
          area_rec = 0;
      end
      if(sum_pts == 4)
          %         disp([k,j])
         vd_x = vdist(lat(k),lon(j),lat(k),lon(j+1))/1000;
         vd_y = vdist(lat(k),lon(j),lat(k+1),lon(j))/1000;
         area_rec = vd_x*vd_y;  
      else
         area_rec = 0;
      end
      area_mtrx(k,j) = area_rec;
   end
end


ind = find(area_mtrx > 0);
area_all = sum(area_mtrx(ind));