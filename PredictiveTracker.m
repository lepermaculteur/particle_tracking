function [vtracks,ntracks,meanlength,rmslength] = PredictiveTracker(inputnames, bground_name, framerange, bigstrel,threshold,smallstrel,arealim,thresh_focus,max_disp,outputname)
% Usage: [vtracks,ntracks,meanlength,rmslength] = PredictiveTracker(inputnames,threshold,max_disp,[bground_name],[minarea],[invert],[noisy])
% Given a movie of particle motions, PredictiveTracker produces Lagrangian
% particle tracks using a predictive three-frame best-estimate algorithm.
% The movie must be saved as a series of image files, an image stack in
% .tif or .gif format, or an uncompressed .avi file; specify the movie in
% "inputnames" (e.g., '0*.png' or 'stack.tif', or 'movie.avi'). To be
% identified as a particle, a part of the image must have brightness that
% differs from the background by at least "threshold". If invert==0,
% PredictiveTracker seeks particles brighter than the background; if
% invert==1, PredictiveTracker seeks particles darker than the background;
% and if invert==-1, PredictiveTracker seeks any sort of contrast. The
% background is read from the file "bground_name"; see BackgroundImage. If
% minarea==1, PredictiveTracker seeks single-pixel particles by comparing
% brightness to adjacent pixels (fast and good for small particles);
% otherwise PredictiveTracker seeks particles having areas larger than
% "minarea" (in square pixels; this method is better for tracking large
% particles). Once identified, each particle is tracked using a kinematic
% prediction, and a track is broken when no particle lies within "max_disp"
% pixels of the predicted location. The results are returned in the
% structure "vtracks", whose fields "len", "X", "Y", "T", "U", and "V"
% contain the length, horizontal coordinates, vertical coordinates, times,
% horizontal velocities, and vertical velocities of each track,
% respectively. If minarea~=1, "vtracks" is returned with an additonal
% field, "Theta", giving the orientation of the major axis of the particle
% with respect to the x-axis, in radians. The total number of tracks is
% returned as "ntracks"; the mean and root-mean-square track lengths are
% returned in "meanlength" and "rmslength", respectively. If noisy~=0, the
% movie is repeated with overlaid velocity quivers and the tracks are
% plotted. If noisy==2, each movie frame is also saved to disk as an image.
% Requires ParticleFinder.m; also requires read_uncompressed_avi.m for use
% with .avi movies. This file can be downloaded from
% http://leviathan.eng.yale.edu/software.

% Written by Nicholas T. Ouellette September 2010. 
% Updated by Douglas H. Kelley 13 April 2011 to plot tracks. 
% Fixed bug in accounting active tracks 14 April 2011. 
% Allowed for inputnames and/or bground_name in other directory, 4 May 
% 2011. 
% Fixed minor bug providing default input parameters 27 May 2011. 
% Added movie with overlaid velocity quivers 17 August 2011. 
% Added "finder" option 18 August 2011. 
% Changed "finder" to "minarea", added "invert" option, and combined two
% plots into one 1 September 2011. 
% Enabled noisy>1 for saving a movie 2 September 2011. 
% Added "Theta" field (if minarea~=1) 7 September 2011. 
% Made compatible with uncompressed avi movies (using 
% read_uncompressed_avi.m) 19 October 2011. 
% Changed plotting method to use color and show only active tracks, 9 
% January 2012. 
% Improved compatibility with images having more than 255 grays, 10 
% January 2012. 
% Made compatible with tiff and gif stacks 13 February 2012. 
% Fixed bugs associated with frames in which no particles are found,
% 6 March 2012.
% Offloaded particle finding to ParticleFinder.m, 7 March 2012.

% -=- Set defaults -=-----------------------------------------------------
framerange_default = [1 inf]; % by default, all frames
bigstrel_default=5;
threshold_default=45;
smallstrel_default=2;
arealim_default=1; 
thresh_focus_default=15
pausetime=1/30; % seconds to pause between frames when plotting
filterwidth = 1; fitwidth = 0; % parameters for the differentiation kernel
minarea=arealim(2);
% -=- Parse inputs -=-----------------------------------------------------
if nargin<1
    error(['Usage: [vtracks,ntracks,meanlength,rmslength] = ' mfilename ...
        '(inputnames,threshold,max_disp,[bground_name],[minarea],[invert],[noisy])'])
end
if ~exist('bground_name','var') || isempty(bground_name)
    fprintf('Il est où le bhackground freroooo ?!?\n')
    return
end
if ~exist('framerange','var') || isempty(framerange)
    framerange=framerange_default;
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
if ~exist('output','var') || isempty(output)
    fprintf("No output input, output set to default :\n")
    outputname=0
end

% -=- Find particles in all frames -=-------------------------------------
[x,y,t,area,eccentricity,focus_mes] = ParticleFinder(inputnames,bground_name,framerange,bigstrel,threshold,smallstrel,arealim,thresh_focus,outputname);
[tt,begins]=unique(t);
ends=circshift(begins,[-1 0])-1;
if length(ends)>0
  ends(end)=length(t);
end
Nf=numel(tt);
if Nf < (2*fitwidth+1)
##    error(['Sorry, found too few files named ' inputnames '.'])
    ntracks=0;
    meanlength=0;
    rmslength=0;
    vtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'Area',[],'Eccentricity',[],'Focus',[]),ntracks,1);
    return
end

% -=- Set up struct array for tracks -=-----------------------------------
ind=begins(1):ends(1);
nparticles = numel(ind);
if minarea==1
    tracks = repmat(struct('len',[],'X',[],'Y',[],'T',[]),nparticles,1);
    for ii = 1:nparticles
        tracks(ii) = struct('len',1,'X',x(ind(ii)),'Y',y(ind(ii)),'T',1);
    end
else
    tracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'Area',[],'Eccentricity',[],'Focus',[]), ...
        nparticles,1);
    for ii = 1:nparticles
        tracks(ii) = struct('len',1,'X',x(ind(ii)),'Y',y(ind(ii)),'T',1, ...
            'Area',area(ind(ii)), 'Eccentricity',eccentricity(ind(ii)), 'Focus',focus_mes(ind(ii)));
    end
end

% -=- Keep track of which tracks are active -=----------------------------
active = 1:nparticles;
n_active = numel(active);
disp(['Processed frame 1 of ' num2str(Nf) '.'])
disp(['    Number of particles found: ' num2str(nparticles,'%.0f')])
disp(['    Number of active tracks: ' num2str(n_active,'%.0f')])
disp(['    Total number of tracks: ' num2str(numel(tracks),'%.0f')])

% -=- Loop over frames -=-------------------------------------------------
for t = 2:Nf
         
    ind=begins(t):ends(t);
    nfr1 = numel(ind);
    if nfr1==0
        warning('MATLAB:PredictiveTracker:noParticles', ...
            ['Found no particles in frame ' num2str(t) '.']);
    end % if nfr1==0
    
    format = size(x(ind));
    if format(2)>format(1)
      fr1=[x(ind)' y(ind)'];
    else
      fr1=[x(ind) y(ind)];
    end
      
    if minarea~=1
        area_t=area(ind);
        eccentricity_t=eccentricity(ind);
        focus_mes_t=focus_mes(ind);
    end

% -=- Match the tracks with kinematic predictions -=----------------------

    % for convenience, we'll grab the relevant positions from the tracks
    now = zeros(n_active,2);
    prior = zeros(n_active,2);
    for ii = 1:n_active
        tr = tracks(active(ii));
        now(ii,1) = tr.X(end);
        now(ii,2) = tr.Y(end);
        if tr.len > 1
            prior(ii,1) = tr.X(end-1);
            prior(ii,2) = tr.Y(end-1);
        else
            prior(ii,:) = now(ii,:);
        end
    end

    % estimate a velocity for each particle in fr0
    velocity = now - prior;
    % and use kinematics to estimate a future position
    estimate = now + velocity;

    % define cost and link arrays
    costs = zeros(n_active,1);
    links = zeros(n_active,1);

    if nfr1>0
        % loop over active tracks
        for ii = 1:n_active
            % now, compare this estimated positions with particles in fr1
            dist_fr1 = (estimate(ii,1)-fr1(:,1)).^2 + (estimate(ii,2)-fr1(:,2)).^2;
            % save its cost and best match
            costs(ii) = min(dist_fr1);
            if costs(ii) > max_disp^2
                continue;
            end
            bestmatch = find(dist_fr1 == costs(ii));
            % if there is more than one best match, we are confused; stop
            if numel(bestmatch) ~= 1
                continue;
            end
            % has another track already matched to this particle?
            ind = links == bestmatch;
            if sum(ind) ~= 0
                if costs(ind) > costs(ii)
                    % this match is better
                    links(ind) = 0;
                else
                    continue;
                end
            end
            links(ii) = bestmatch;
        end

        % now attach the matched particles to their tracks
        matched = zeros(nfr1,1);
        for ii = 1:n_active
            if links(ii) ~= 0 
                % this track found a match
                tracks(active(ii)).X(end+1) = fr1(links(ii),1);
                tracks(active(ii)).Y(end+1) = fr1(links(ii),2);
                tracks(active(ii)).len = tracks(active(ii)).len + 1;
                tracks(active(ii)).T(end+1) = t;
                if minarea~=1
                    tracks(active(ii)).Area(end+1) = area_t(links(ii));
                    tracks(active(ii)).Eccentricity(end+1) = eccentricity_t(links(ii));
                    tracks(active(ii)).Focus(end+1) = focus_mes_t(links(ii));                    
                end
                matched(links(ii)) = 1;
            end
        end
        active = active(links~=0);

        % and start new tracks with the particles in fr1 that found no match
        unmatched = find(matched == 0);
        if minarea==1
            newtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[]), ...
                numel(unmatched),1);
            for ii = 1:numel(unmatched)
                newtracks(ii) = struct('len',1,'X',fr1(unmatched(ii),1),...
                    'Y',fr1(unmatched(ii),2),'T',t);
            end
        else
            newtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[], ...
                'Area',[], 'Eccentricity',[], 'Focus',[]),numel(unmatched),1);
            for ii = 1:numel(unmatched)              
                newtracks(ii) = struct('len',1,'X',fr1(unmatched(ii),1),...
                    'Y',fr1(unmatched(ii),2),'T',t,'Area',area_t(unmatched(ii)),'Eccentricity',eccentricity_t(unmatched(ii)), 'Focus',focus_mes_t(unmatched(ii)));
            end
        end
    else % if nfr1>0
        active=[];
        newtracks=[];
        unmatched=[];
    end
    active = [active (numel(tracks)+1):(numel(tracks)+numel(newtracks))];
    tracks = [tracks ; newtracks];
    n_active = numel(active);

    disp(['Processed frame ' num2str(t) ' of ' num2str(Nf) '.'])
    disp(['    Number of particles found: ' num2str(nfr1,'%.0f')])
    disp(['    Number of active tracks: ' num2str(n_active,'%.0f')])
    disp(['    Number of new tracks started here: ' ...
        num2str(numel(unmatched),'%.0f')])
    disp(['    Number of tracks that found no match: ' ...
        num2str(sum(links==0),'%.0f')])
    disp(['    Total number of tracks: ' num2str(numel(tracks),'%.0f')])

end

% -=- Prune tracks that are too short -=----------------------------------
disp('Pruning...');
##tracks = tracks([tracks.len] >= (2*fitwidth+1));
tracks = tracks([tracks.len] >= 5);
ntracks = numel(tracks);
meanlength = mean([tracks.len]);
rmslength = sqrt(mean([tracks.len].^2));

% -=- Compute velocities -=-----------------------------------------------
disp('Differentiating...');

% define the convolution kernel
##Av = 1.0/(0.5*filterwidth^2 * ...
##    (sqrt(pi)*filterwidth*erf(fitwidth/filterwidth) - ...
##    2*fitwidth*exp(-fitwidth^2/filterwidth^2)));
##vkernel = -fitwidth:fitwidth;
##vkernel = Av.*vkernel.*exp(-vkernel.^2./filterwidth^2);

% loop over tracks
if minarea==1
    vtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'U',[],'V',[]),ntracks,1);
else
##    vtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'U',[], ...
##        'V',[],'Theta',[]),ntracks,1);
    vtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'Area',[],'Eccentricity',[],'Focus',[]),ntracks,1);
end
for ii = 1:ntracks
##    u = -conv(tracks(ii).X,vkernel,'valid');
##    v = -conv(tracks(ii).Y,vkernel,'valid');
    if minarea==1
        vtracks(ii) = struct('len',tracks(ii).len - 2*fitwidth, ...
            'X',tracks(ii).X(fitwidth+1:end-fitwidth), ...
            'Y',tracks(ii).Y(fitwidth+1:end-fitwidth), ...
            'T',tracks(ii).T(fitwidth+1:end-fitwidth), ...
            'U',u, ...
            'V',v);
    else
##        vtracks(ii) = struct('len',tracks(ii).len - 2*fitwidth, ...
##            'X',tracks(ii).X(fitwidth+1:end-fitwidth), ...
##            'Y',tracks(ii).Y(fitwidth+1:end-fitwidth), ...
##            'T',tracks(ii).T(fitwidth+1:end-fitwidth), ...
##            'U',u, ...
##            'V',v, ...
##            'Theta',tracks(ii).Theta(fitwidth+1:end-fitwidth));
        vtracks(ii) = struct('len',tracks(ii).len - 2*fitwidth, ...
            'X',tracks(ii).X(fitwidth+1:end-fitwidth), ...
            'Y',tracks(ii).Y(fitwidth+1:end-fitwidth), ...
            'T',tracks(ii).T(fitwidth+1:end-fitwidth), ...
            'Area',tracks(ii).Area(fitwidth+1:end-fitwidth), ...
            'Eccentricity',tracks(ii).Eccentricity(fitwidth+1:end-fitwidth), ...
            'Focus',tracks(ii).Focus(fitwidth+1:end-fitwidth));
    end
end
disp('Done.')

end

