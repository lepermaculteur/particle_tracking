function [tracks,ntracks,meanlength,rmslength] = PredictiveTracker(inputnames, bground_name, framerange, bigstrel,threshold,smallstrel,arealim,thresh_focus,max_disp,outputname)
% Usage: [tracks,ntracks,meanlength,rmslength] = PredictiveTracker(inputnames, bground_name, framerange, bigstrel,threshold,smallstrel,arealim,thresh_focus,max_disp,outputname)

%Adaptation and extension of original code written by Nicholas T. Ouellette, September 2010
% 
% Noticeable add-on (February 2022): 
% -Added a selection on the eccentricity of the detected regions to decide if they are real particles or not
% -Implementation of a "in-focus" threshold 
%   For each particle, a small square sample of 12x12pixels is taken from the background removed image  
%   (centered on the particle center) and sent to a "focus measurement algorithm" (see 'Helmli and Scherer's 
%   mean method').If the value returned is below a certain threshold, the particle is considered as 'out of focus' and 
%   is rejected. The "focus measurement algorithm" has been selected among about twenty other algorithms presented in
%   "Pertuz et al. / Pattern Recognition (2013)" after being tested on calibration images.
% -Correction of  minor indexation errors in the initialization and finalization tracking step 

% -=- Set defaults -=-----------------------------------------------------
framerange_default = [1 inf]; % by default, all frames
bigstrel_default=5;
threshold_default=45;
smallstrel_default=2;
arealim_default=[0 7]; 
thresh_focus_default=15
max_disp_default=2;

% -=- Parse inputs -=-----------------------------------------------------
if nargin<1
    error(['Usage: [tracks,ntracks,meanlength,rmslength] = ' mfilename ...
        '(inputnames, bground_name, framerange, bigstrel,threshold,'...
        'smallstrel,arealim,thresh_focus,max_disp,outputname)'])
end
if ~exist('inputnames','var') || isempty(inputnames)
    fprintf('You should provide a inputnames of your images\n')
    return
end
if ~exist('bground_name','var') || isempty(bground_name)
    fprintf('You should provide a background\n')
    return
end
if ~exist('framerange','var') || isempty(framerange)
    fprintf("No framerange input, framerange set to default :\n")
    framerange=framerange_default
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
if ~exist('thresh_focus','var') || isempty(thresh_focus)
    fprintf("No thresh_focus input, thresh_focus set to default :\n")
    thresh_focus=thresh_focus_default;
end
if ~exist('max_disp','var') || isempty(max_disp)
    fprintf("No max_disp input, max_disp set to default :\n")
    max_disp=max_disp_default;
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
if Nf < 5
##    error(['Sorry, found too few files named ' inputnames '.'])
    ntracks=0;
    meanlength=0;
    rmslength=0;
    tracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'Area',[],'Eccentricity',[],'Focus',[]),ntracks,1);
    return
end

% -=- Set up struct array for tracks -=-----------------------------------
ind=begins(1):ends(1);
nparticles = numel(ind);
tracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'Area',[],'Eccentricity',[],'Focus',[]), ...
    nparticles,1);
for ii = 1:nparticles
    tracks(ii) = struct('len',1,'X',x(ind(ii)),'Y',y(ind(ii)),'T',1, ...
        'Area',area(ind(ii)), 'Eccentricity',eccentricity(ind(ii)), 'Focus',focus_mes(ind(ii)));
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
     
    area_t=area(ind);
    eccentricity_t=eccentricity(ind);
    focus_mes_t=focus_mes(ind);


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

                tracks(active(ii)).Area(end+1) = area_t(links(ii));
                tracks(active(ii)).Eccentricity(end+1) = eccentricity_t(links(ii));
                tracks(active(ii)).Focus(end+1) = focus_mes_t(links(ii));                    

                matched(links(ii)) = 1;
            end
        end
        active = active(links~=0);

        % and start new tracks with the particles in fr1 that found no match
        unmatched = find(matched == 0);
        
        newtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[], ...
            'Area',[], 'Eccentricity',[], 'Focus',[]),numel(unmatched),1);
        for ii = 1:numel(unmatched)              
            newtracks(ii) = struct('len',1,'X',fr1(unmatched(ii),1),...
                'Y',fr1(unmatched(ii),2),'T',t,'Area',area_t(unmatched(ii)),'Eccentricity',eccentricity_t(unmatched(ii)), 'Focus',focus_mes_t(unmatched(ii)));
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
tracks = tracks([tracks.len] >= 5);
ntracks = numel(tracks);
meanlength = mean([tracks.len]);
rmslength = sqrt(mean([tracks.len].^2));


% loop over tracks
tracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'Area',[],'Eccentricity',[],'Focus',[]),ntracks,1);

for ii = 1:ntracks
  tracks(ii) = struct('len',tracks(ii).len - 2*0, ...
      'X',tracks(ii).X, ...
      'Y',tracks(ii).Y, ...
      'T',tracks(ii).T, ...
      'Area',tracks(ii).Area, ...
      'Eccentricity',tracks(ii).Eccentricity, ...
      'Focus',tracks(ii).Focus);
end
disp('Done.')

end

