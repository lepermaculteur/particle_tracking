function [x,y,t,area,eccentricity,focus_mes] = ParticleFinder(inputnames,bgname,framerange,bigstrel,threshold,smallstrel,arealim,thresh_focus,outputname)

% Usage: [x,y,t,ang] = ParticleFinder(inputnames,threshold,[framerange],[outputname],[bground_name],[arealim],[invert],[noisy])
% Given a movie of particle motions, ParticleFinder identifies the
% particles, returning their positions, times, and orientations in x, y,
% and t, respectively. The movie must be saved as a series of image files,
% an image stack in .tif or .gif format, or an uncompressed .avi file;
% specify the movie in "inputnames" (e.g., '0*.png' or 'stack.tif', or
% 'movie.avi'). To be identified as a particle, a part of the image must
% have brightness that differs from the background by at least "threshold".
% If invert==0, ParticleFinder seeks particles brighter than the
% background; if invert==1, ParticleFinder seeks particles darker than the
% background; and if invert==-1, ParticleFinder seeks any sort of contrast.
% The background is read from the file "bground_name"; see BackgroundImage.
% Frames outside the range specified by the two-element vector "framerange"
% are ignored. If arealim==1, ParticleFinder seeks single-pixel particles
% by comparing brightness to adjacent pixels (fast and good for small
% particles); otherwise ParticleFinder seeks particles having areas bounded
% by the two elements of the vector "arealim" (in square pixels; this
% method is better for tracking large particles). If "outputname" is not
% empty, particle positions are also saved as a binary file of that name.
% The file begins with a long int giving the frame count, then each frame
% begins with a long int giving its particle count, and continues with
% double-precision floats giving the x and y coordinates of each particle.
% If noisy~=0, the movie is repeated with particle locations overlaid. If
% noisy>1, each movie frame is also saved to disk as an image. See also
% BackgroundImage.m and PredictiveTracker.m. Requires
% read_uncompressed_avi.m for use with .avi movies. 
%
% Written 20 October 2011 by Doug Kelley, largely based on 
% PredictiveTracker.m.
% Renamed ParticleFinder and incorporated FindParticles function 27 October 
% 2011. 
% Fixed plotting bug in bug plotting 15 November 2011. 
% Added "ang" output (particle orientation) 18 November 2011. 
% Added "framerange" input 28 November 2011. 
% Updated to use weighted centroid 2 December 2011. 
% Included invert==-1 option 22 February 2012.
% Made compatible with tiff & gif stacks and squelched regionprops 
% divide-by-zero warning 7 March 2012. 
% Thanks go to Max Gould for spotting a bug in the main loop; fixed 26
% March 2012. 
% Changed 29 March 2012 to plot bars if minarea>1. 
% Updated 30 March 2012 to use *weighted* orientaion and pre-allocate
% arrays.
% Updated 29 May 2012 with "arealim" input instead of "minearea".

% Next: write angle to output file!

% -=- Set defaults -=-----------------------------------------------------

framerange_default = [1 inf]; % by default, all frames
bigstrel_default=5;
threshold_default=45;
smallstrel_default=2;
arealim_default=1; 
thresh_focus_default=15
pausetime=1/30; % seconds to pause between frames when plotting
if nargin<2
    error(['Usage: [x,y,t,ang] = ' mfilename ...
        '(inputnames,threshold,[framerange],[outputname],' ...
        '[bground_name],[arealim],[invert],[noisy])'])
end
if ~exist('framerange','var') || isempty(framerange)
    fprintf("No framerange input, framerange set to default :\n")
    framerange=framerange_default
elseif numel(framerange)==1
    framerange=framerange*[1 1];
end
if ~exist('bigstrel','var') || isempty(bigstrel)
    fprintf("No bigstrel input, bigstrel set to default :\n")
    bigstrel=bigstrel_default
end
if ~exist('threshold','var') || isempty(threshold)
    fprintf("No threshold input, threshold set to default :\n")
    threshold=threshold_default
end
if ~exist('smallstrel','var') || isempty(smallstrel)
    fprintf("No smallstrel input, smallstrel set to default :\n")
    smallstrel=smallstrel_default
end
if ~exist('arealim','var') || isempty(arealim)
    fprintf("No arealim input, arealim set to default :\n")
    arealim=arealim_default;
end
if ~exist('bgname','var') || isempty(bgname)
    fprintf('Il est où le bhackground freroooo ?!?\n')
    return
end
if ~exist('outputname','var') || isempty(outputname) || outputname==0
    writefile=false;
else
    writefile=true;
end

% -=- Decide whether avi, stack, or images; set up -=---------------------
[filepath,junk,ext]=fileparts(inputnames);
names=dir(inputnames);
if strcmpi(ext,'.avi')
    movtype='avi';
    if isempty(which('read_uncompressed_avi.m')) % check for req'd helper function
        error(['Sorry, reading .avi files requires ' ...
            'read_uncompressed_avi.m.'])
    end
    movinfo=aviinfo(fullfile(filepath,names.name));
    color_depth=movinfo.NumColormapEntries;
    ht=movinfo.Height;
    wd=movinfo.Width;
    tmin=max([framerange(1) 1]);
    tmax=min([framerange(2) movinfo.NumFrames]);
elseif numel(names)==1 && ( strcmpi(ext,'.tif') || ...
        strcmpi(ext,'.tiff') || strcmpi(ext,'.gif') ) % single file, looks like an image stack
    movtype='stack';
    movinfo=imfinfo(fullfile(filepath,names.name));
    color_depth=2^(movinfo(1).BitDepth);
    ht=movinfo(1).Height;
    wd=movinfo(1).Width;
    tmin=max([framerange(1) 1]);
    tmax=min([framerange(2) movinfo.NumFrames]);
else
    movtype='images';
    movinfo=imfinfo(fullfile(filepath,names(1).name));
    color_depth=2^(movinfo.BitDepth);
    ht=movinfo.Height;
    wd=movinfo.Width;
    for ii=1:numel(names);
        if strcmp(fullfile(filepath,names(ii).name),bgname)
            names(ii)=[]; % don't try to track the background file
            break
        end
    end
    tt=NaN(size(names));
    for ii=1:numel(names)
        [junk,myname]=fileparts(names(ii).name);
        tt(ii)=str2double(myname);
    end
    if any(isnan(tt)) % at least one filename was not a number, so just make ordinals
        tt=1:numel(tt);
    end
    tmin=max([framerange(1) min(tt)]);
    tmax=min([framerange(2) max(tt)]);
    names = names( (tt>=tmin) & (tt<=tmax) );
end
Nf=tmax-tmin+1 % frame count

% -=- Pre-compute logarithms for locating particle centers -=-------------
if arealim==1
    logs = 1:color_depth;
    logs = [log(0.0001) log(logs)];
end

%Read background
bk=imread(bgname);

% -=- Set up output file and variables -=---------------------------------
if writefile
    if exist(outputname,'file')
        yn=input(['File ' outputname ' exists. Overwrite (y/n)? '],'s');
        if ~strcmpi(yn(1),'y')
            disp('Not overwritten.')
            writefile=false;
        else
            disp(['Replacing file ' outputname '.'])
            fid=fopen(outputname,'w');
            fwrite(fid,Nf,'int32'); % header is integer frame count
        end
        pause(1)
    else
        fid=fopen(outputname,'w');
        fwrite(fid,Nf,'int32'); % header is integer frame count
    end
end
begins=ones(Nf+1,1);
mark=[];

for ii=1:Nf
    tt=tmin+ii-1; % current time
    switch movtype
        case('avi')
            [im,mark] = read_uncompressed_avi( ...
                fullfile(filepath,names.name),tt,mark);
            im = double(im.cdata);
        case('stack')
            im = double(imread(fullfile(filepath,names.name),tt));
        case('images')
            im = imread(fullfile(filepath,names(ii).name));
    end % switch movtype
    if ndims(im)==3
        if size(im,3)>3
            error('Sorry, only grayscale and RGB images are supported.')
        end
        im=round(mean(im,3)); % convert to grayscale if necessary
    end
    im_nobk=bk-im;
    ampl=255/max(im_nobk(:));
    nobk_ampl = im_nobk*ampl;
    if arealim==1
        pos = FindParticles(im,threshold,logs);
    else
         tic
         [pos,area_ii,eccentricity_ii,focus_mes_ii] = FindRegions(nobk_ampl,bigstrel,threshold,smallstrel,arealim,thresh_focus);
         toc
    end
    N=size(pos,1);
    if ii==1 % if first frame, pre-allocate arrays for speed
        x=NaN(N*Nf,1);
        y=NaN(N*Nf,1);
        t=NaN(N*Nf,1);
        if arealim~=1
            area=NaN(N*Nf,1);
            eccentricity=NaN(N*Nf,1);
            focus_mes=NaN(N*Nf,1);
        else
            area=[];
            eccentricity=[];
            focus_mes=[]; % no angles
        end
        memloc=1;
    end
    if N>0
        x(memloc:memloc+N-1)=pos(:,1);
        y(memloc:memloc+N-1)=pos(:,2);
        t(memloc:memloc+N-1)=tt;
        if max(arealim)>1
            area(memloc:memloc+N-1)=area_ii;
            eccentricity(memloc:memloc+N-1)=eccentricity_ii;
            focus_mes(memloc:memloc+N-1)=focus_mes_ii;
        end
        memloc=memloc+N;
    end
    begins(ii+1)=begins(ii)+N;
    switch movtype
        case('avi')
            disp(['Found ' num2str(N,'%.0f') ' particles in ' ...
                names.name ' frame ' num2str(tt) ' ('  ...
                num2str(ii) ' of ' num2str(Nf) ').'])
        case('stack')
            disp(['Found ' num2str(N,'%.0f') ' particles in ' ...
                names.name ' frame ' num2str(tt) ' (' ...
                num2str(ii) ' of ' num2str(Nf) ').'])
        case('images')
            disp(['Found ' num2str(N,'%.0f') ' particles in ' ...
                names(ii).name ' (' num2str(ii) ' of ' num2str(Nf) ').'])
    end % switch movtype
    if writefile
        fwrite(fid,N,'int32'); % start the frame w/ particle count
        fwrite(fid,pos','float32'); % then write particle locations x1 y1 x2 y2  ...
        % Next: write area_ii eccentricity_ii focus_mes_ii  as well!
    end
end % for ii=1:tmax
x(memloc:end)=[];
y(memloc:end)=[];
t(memloc:end)=[];
if arealim~=1
    area(memloc:end)=[];
    eccentricity(memloc:end)=[];
    focus_mes(memloc:end)=[];
end
if writefile
    fclose(fid);
end

end % function ParticleFinder


% -=- function FindRegions -=---------------------------------------------
function [pos,area_ii,eccentricity_ii,focus_mes_ii] = FindRegions(nobk_ampl,bigstrel,threshold,smallstrel,arealim,thresh_focus)
% Given an image "im", FindRegions finds regions that are brighter than
% "thresold" and have area larger than "arealim". Region centroids are 
% returned in the two-column array "pos" (with x-coordinates in the first
% column and y-coordinates in the second). Region orientations are 
% returned in radians, in the vector "ang".
    %%Remove anything bigger than a disk of radius 'bigstrel'
      dd = strel ("disk", bigstrel, 0);   
      imod = imopen(nobk_ampl,dd);
      im_op=nobk_ampl-imod; 
      ampl_contrast=255/max(max(im_op));
      im_op_ampl=im_op*ampl_contrast;  

    %%Applying a 'threshold' to remove the background noise
      im_thresh=im_op_ampl;
      im_thresh(im_thresh<threshold)=0;

    %%Remove anything smaller than a disk of radius 'smallstrel'
      dd = strel ("disk", smallstrel, 0);   
      im_finale = imopen(im_thresh,dd); 
      ampl_contrast=255/max(max(im_finale));
      im_finale_ampl=im_finale*ampl_contrast;  

    %%Find Particules in the processed image
    part=regionprops(im2bw(im_finale_ampl,0), 'Centroid', 'Area', 'Eccentricity')

    %Get the position of particles, their area and eccentricity
    pos=vertcat(part.Centroid);
    area_ii=vertcat(part.Area);
    eccentricity_ii=vertcat(part.Eccentricity);
    
    % remove regions  too small, or too big
    if numel(pos)>0
        good = area_ii>arealim(1) & area_ii<arealim(2) & eccentricity_ii<0.8;
          
        part=part(good);
        pos=pos(good,:);
        area_ii=area_ii(good);
        eccentricity_ii=eccentricity_ii(good);
    end
    
    nb_part=length(part);
    
    %Grab a small square around each particles and compute the focus estimation of each paritcles
     focus_measures = {'ACMO', 'BREN', 'CURV', 'GDER', ...
       'GLVA', 'GLLV', 'GLVN', 'GRAE', 'GRAT', 'GRAS', 'HELM', ...
       'HISE', 'HISR', 'LAPE', 'LAPM', 'LAPV', 'LAPD', 'SFIL', ...
       'SFRQ', 'TENG', 'TENV', 'VOLA'};
       
    choosen_method='HELM';
    
    rect_size=12;
    rect=zeros(nb_part,4);
    for pj=1:nb_part
      
      anchor=[pos(pj,1)-rect_size pos(pj,2)-rect_size];
      rect(pj,:) = [anchor 2*rect_size 2*rect_size];
    end
    rect_bound=rect;
    rect_bound(:,3)=rect_bound(:,1)+2*rect_size;
    rect_bound(:,4)=rect_bound(:,2)+2*rect_size;
    rect_bound(rect_bound<0.5)=1;
    rect_bound(rect_bound>=1023+1.5)=1023+1;
    rect_bound = round(rect_bound);
    
    focus_part_ii=zeros(nb_part,1);  
    for pj=1:nb_part
      sample=nobk_ampl(rect_bound(pj,2):rect_bound(pj,4), rect_bound(pj,1):rect_bound(pj,3));
      focus_part_ii(pj)=fmeasure(sample, choosen_method);
    end
    
    idx = focus_part_ii>thresh_focus;
    
    pos=pos(idx,:);
    area_ii=area_ii(idx);
    eccentricity_ii=eccentricity_ii(idx);
    focus_mes_ii=focus_part_ii(idx);
   
end % function FindRegions



% -=- function FindParticles -=-------------------------------------------
function pos = FindParticles(im, threshold, logs)
% Given an image "im", FindParticles finds small particles that are 
% brighter than their four nearest neighbors and also brighter than
% "threshold". Particles are located to sub-pixel accuracy by applying a 
% Gaussian fit in each spatial direction. The input "logs" depends on 
% the color depth and is re-used for speed. Particle locations are 
% returned in the two-column array "pos" (with x-coordinates in the first
% column and y-coordinates in the second). 

    s = size(im);

    % identify the local maxima that are above threshold  
    maxes = find(im >= threshold & ...
        im > circshift(im,[0 1]) & ...
        im > circshift(im,[0 -1]) & ...
        im > circshift(im,[1 0]) & ...
        im > circshift(im,[-1 0]));

    % now turn these into subscripts
    [x,y] = ind2sub(s, maxes);

    % throw out unreliable maxes in the outer ring
    good = find(x~=1 & y~=1 & x~=s(1) & y~=s(2));
    x = x(good);
    y = y(good);

    % find the horizontal positions

    % look up the logarithms of the relevant image intensities
    z1 = logs(im(sub2ind(s,x-1,y)) + 1)';
    z2 = logs(im(sub2ind(s,x,y)) + 1)';
    z3 = logs(im(sub2ind(s,x+1,y)) + 1)';

    % compute the centers
    xcenters = -0.5 * (z1.*(-2*x-1) + z2.*(4*x) + z3.*(-2*x+1)) ./ ...
        (z1 + z3 - 2*z2);

    % do the same for the vertical position
    z1 = logs(im(sub2ind(s,x,y-1)) + 1)';
    z3 = logs(im(sub2ind(s,x,y+1)) + 1)';
    ycenters = -0.5 * (z1.*(-2*y-1) + z2.*(4*y) + z3.*(-2*y+1)) ./ ...
        (z1 + z3 - 2*z2);

    % make sure we have no bad points
    good = find(isfinite(xcenters) & isfinite(ycenters));

    % fix up the funny coordinate system used by matlab
    pos = [ycenters(good), xcenters(good)];

end % function FindParticles

