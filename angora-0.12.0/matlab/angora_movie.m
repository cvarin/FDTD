% Copyright (c) 2012, Ilker R Capoglu 
% This script reads an Angora movie file and plots the data.

fclose('all');
clear all;

framejump = 1;
% framelimits = [320 320];
draw_geometry = false;
showcolorbar = false;
flip_colormap = false;
% record AVI movie?
record_movie = false;

% full path to file
moviefilename = '/my_path/MovieFile_Ey_yz_0_0.amv';
% open the file
moviefile = fopen(moviefilename);
if (moviefile==-1)
    error(['Error: Cannot open movie file ' moviefilename ]);
end

% read the package version
version_major = fread(moviefile ,1,'int');
version_minor = fread(moviefile ,1,'int');
version_rev = fread(moviefile ,1,'int');

numofbytes = fread(moviefile ,1,'int');
dx = fread(moviefile ,1,'double');
dt = fread(moviefile ,1,'double');
initial_time_value = fread(moviefile ,1,'double');
maxvalue = fread(moviefile ,1,'double');
minvalue = fread(moviefile ,1,'double');    
length_1 = fread(moviefile ,1,'int');   %this is the first dimension of the recorded array (x if xy section is recorded, for example)
length_2 = fread(moviefile ,1,'int');   %this is the first dimension of the recorded array (y if xy section is recorded, for example)
length_time = fread(moviefile ,1,'int');
PML = fread(moviefile ,1,'int');
dim1_range = fread(moviefile,length_1,'double');    % first dimension has length "length_1"
dim2_range = fread(moviefile,length_2,'double');    % second dimension has length "length_2"

% is the number of bytes valid (1 or 8?)
if ((numofbytes~=1)&&(numofbytes~=8))
    error(['Number of bytes (' num2str(numofbytes) ') invalid']);
end
% is the field indexed (from 0 to 255) or absolute?
if (numofbytes>1)
    indexed_field = false;
else
    indexed_field = true;
end
% override the max. and min. values? (if recording is in double)
if (numofbytes==8)
    maxvalue = 310;
    minvalue = maxvalue-60; 
end

if (~exist('framelimits'))
    frameoffsets = 0:framejump:length_time;
else
    frameoffsets = framelimits(1):framejump:framelimits(2);
end

% read material properties
num_of_properties = 2;
property_array = zeros(length_2,length_1,num_of_properties);
for i=1:num_of_properties
    % In the C++ code, the second dimension is written first (due to the
    % row-major ordering of the blitz++ array library). Here, this
    % dimension is read into the "columns" (first dimension) of the matlab
    % array.
    property_array(:,:,i)= fread(moviefile ,[length_2,length_1],'double');
end
property_array = property_array(PML+1:length_2-PML,PML+1:length_1-PML,:);
rel_perm = property_array(:,:,1);   %relative permittivity
refr_index = sqrt(rel_perm);    %refractive index
cond = property_array(:,:,2);   %conductivity

% In the C++ code, the second dimension is written first, because of the
% row-major ordering of the blitz++ array library. The function "fread"
% reads this information column-wise into the matlab array.
% Therefore, the "height" of the image (size of each column, or number of
% rows) is length_2 (the size of the second dimension in the recorded field
% array).
plot_width = length_1-2*PML;
plot_height = length_2-2*PML;
width_range = dim1_range(PML+1:length_1-PML);
height_range = dim2_range(PML+1:length_2-PML);

ncolors = 256;

%some parameters related to the geometry
n0 = 1.38;
dn = 0.005;

geom = uint8(zeros(plot_height,plot_width,3));
if(draw_geometry)
    blendfactor = 1;
    dn_norm = (refr_index-n0)/n0;
    max_dn_norm = 1*dn;
    min_dn_norm = -1*dn;
    dn_norm_discr = round(ncolors*(-0.5+(dn_norm-min_dn_norm)/(max_dn_norm-min_dn_norm)));
    % regions where normalized fluctuation is positive and negative, assign different colors
    dn_norm_discr_pos = zeros(size(dn_norm_discr));
    dn_norm_discr_neg = zeros(size(dn_norm_discr));
    dn_norm_discr_pos(dn_norm_discr>0) = dn_norm_discr(dn_norm_discr>0);
	dn_norm_discr_neg(dn_norm_discr<0) = dn_norm_discr(dn_norm_discr<0);
    geom(:,:,1) = abs(dn_norm_discr_pos);    % red
	geom(:,:,2) = abs(dn_norm_discr_neg);    % green
    geom(:,:,3) = 255*cond/max(max(cond));   % blue
end

initial_offset = ftell(moviefile);
h=figure('DoubleBuffer','on');

rgb_array = uint8(zeros(plot_height,plot_width,3));
for n=1:length(frameoffsets)
    H=image(width_range*1e6,height_range*1e6,zeros(plot_height,plot_width),'CDataMapping','direct');

    fontsize = 16;
    set(gca,'LineWidth',1,'FontName','Times','FontSize',fontsize);% ...
%         ,'TickLength',[0.01 0.01],'Units','inches','Position',[0.75,0.9,10,5]);%,...
%         'XTick',0:2:12,'YTick',0:2:16);
    if (~isempty(strfind(moviefilename,'xz')))
        xlabel('\it x \rm (\mum)','FontSize',fontsize);
        ylabel('\it z \rm (\mum)','FontSize',fontsize);
    else if (~isempty(strfind(moviefilename,'yz')))
        xlabel('\it y \rm (\mum)','FontSize',fontsize);
        ylabel('\it z \rm (\mum)','FontSize',fontsize);
    else if (~isempty(strfind(moviefilename,'xy')))
        xlabel('\it x \rm (\mum)','FontSize',fontsize);
        ylabel('\it y \rm (\mum)','FontSize',fontsize);
    end,end,end            
    axis([1 plot_height 1 plot_width]);
    
    colormap(gray(ncolors));
    if (flip_colormap)
        colormap(flipud(colormap));
    end
    if (showcolorbar)
        cbar = colorbar('Location','EastOutside');
        cbarfontsize = 15;
        set(cbar,'FontName','Times','FontSize',cbarfontsize);
        title(cbar,'dB');
    end

    caxis([minvalue maxvalue]);
    if (showcolorbar)
        set(cbar,'YTick',[round(minvalue):10:round(maxvalue)])
    end
    
    if (n~=length(frameoffsets))
        colorbar off;
    end
    shading interp;
    axis image xy;  % "image" sets the aspect ratio to that of the xy axes
                    % "xy" interprets the dimensions CData matrix as increasing Cartesian coordinates
                    % CData(1,1) becomes the lower-left corner

    pixelsize = numofbytes;
    framesize = length_2*length_1*pixelsize;
    fseek(moviefile,initial_offset+min(frameoffsets(n),length_time-1)*framesize,'bof');
    % read the field array from the C++ output file
    % see the explanation of the array dimensions above
    % basically, the "height" of the image is always equal to the second
    % dimension of the recorded array (z if yz section is recorded, or y if xy section is recorded)
    if (length_time==0)
        indexarray = zeros(length_2,length_1);
    else
        if (numofbytes==1)
            indexarray = fread(moviefile ,[length_2,length_1],'uint8');
        else if (numofbytes==8)
                indexarray = fread(moviefile ,[length_2,length_1],'double');
                % convert to indexed array (between 0 and 255)
                indexarray = uint8(255/(maxvalue-minvalue)*(indexarray-minvalue)+.5);
            end
        end
    end
        
    indexarray = indexarray(PML+1:length_2-PML,PML+1:length_1-PML);

	% convert into RGB (all slices equal to indexarray, therefore in graytone)
    rgb_array(:,:,1) = indexarray;
    rgb_array(:,:,2) = indexarray;
    rgb_array(:,:,3) = indexarray;

    if (draw_geometry)
         rgb_array = rgb_array + geom;
    end
    
    if (length_time==0)
        title('Geometry');
    else
        title(['Time step : ',num2str(min(frameoffsets(n),length_time-1))]);
    end
    
    set(H,'CData',rgb_array);drawnow;
    
    if (record_movie) movieframes(n)=getframe(gcf); end
end
if (record_movie) movie2avi(movieframes,[moviefilename '.avi'],'fps',200/framejump); end

fclose(moviefile);
