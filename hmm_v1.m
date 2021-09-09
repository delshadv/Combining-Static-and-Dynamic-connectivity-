
%% Pre-Processing of MEG data and defining confounds

% Combining Static and Dynamic FC using kernel-based late combination to classify MCI vs Control
% (BioFIND dataset)
%
% This script contains MEG preprocessing steps for reproducing results of the
% paper ""

% Henson R.N 2020, Vaghari D 2020, et al.,

%% Define Paths ands variables

% Assumed you are currently in the directory including BioFIND data,
% OSL and MKL directories as described in readme.md

%restoredefaultpath
bwd = pwd;
wd  = fullfile(bwd,'MKL');
addpath(wd)

% Setup OSL
addpath(fullfile(bwd,'osl','osl-core'))
osl_startup
osl_check_installation

% BIDS and Processed directories
bidspth = fullfile(bwd,'BioFIND','MCIControls'); %BIDS Path
BIDS   = spm_BIDS(bidspth); % (If it does not work with OSL's SPM, so copy last version of spm_BIDS)
subs   = spm_BIDS(BIDS,'subjects', 'task', 'Rest');
nsub   = numel(subs);
subdir = cellfun(@(s) ['sub-' s], subs, 'UniformOutput',false);
procpth = fullfile(bidspth,'derivatives','meg_derivatives'); % If want maxfiltered files

% Define participant variables
participants = spm_load(fullfile(wd,'participants-imputed.tsv'));
group_num    = grp2idx(participants.group);
mri_num      = grp2idx(participants.sImaging);

% Remove noisy-MRIs and non-MRI subjects
mri_num([23 197]) = 2;
y = group_num(mri_num==1);  % Group labels

%% Create Processed directory
% Please do all analysis in a separate directory from BIDS
% Here, we call it "Processed"

processed_pth = fullfile(bwd,'Processed_new');

if ~exist(processed_pth,'dir')
    
    mkdir('Processed_new');
    cd ('Processed_new')
    for s=1:nsub
        mkdir(sprintf('sub-Sub%04d',s))
    end
end

cd (processed_pth)

%% PreProcess- Part 1 (Convert, Downsample, Filter)

parfor sub = 1:nsub
    
    % Read event & json file to extract desirable length of MEG Recordings
    tmp = spm_jsonread(fullfile(procpth,sprintf('sub-Sub%04d',sub),'ses-meg1','meg',[sprintf('sub-Sub%04d',sub) '_ses-meg1_task-Rest_proc-sss_meg.json']));
    event_file = spm_load(fullfile(bidspth,sprintf('sub-Sub%04d',sub),'ses-meg1','meg',[sprintf('sub-Sub%04d',sub) '_ses-meg1_task-Rest_events.tsv']));
    onset = (event_file.onset*tmp.SamplingFrequency)+1;
    
    % offset = onset + event_file.duration*tmp.SamplingFrequency;
    offset = (onset + 120 *tmp.SamplingFrequency)-1; % we put 120 seconds due to min length of raw data
    
    % Converting
    S = [];
    S.outfile = fullfile(processed_pth,sprintf('sub-Sub%04d',sub),'spmeeg');
    S.dataset = fullfile(procpth,sprintf('sub-Sub%04d',sub),'ses-meg1','meg',[sprintf('sub-Sub%04d',sub) '_ses-meg1_task-Rest_proc-sss_meg.fif']);
    S.mode = 'epoched';
    S.channels = {'EOG', 'ECG', 'MEGMAG', 'MEGPLANAR'}; % EEG was removed
    S.checkboundary = 0;
    S.trl = [onset offset 0];
    try
        S.conditionlabels = event_file.stim_type;
    catch
        S.conditionlabels = event_file.trial_type;
    end
    D = spm_eeg_convert(S);
    
    % Set channel types and bad channels
    S = [];
    S.D    = D;
    S.task = 'bidschantype';
    S.save = 1;
    S.filename = fullfile(bidspth,sprintf('sub-Sub%04d',sub),'ses-meg1','meg',[sprintf('sub-Sub%04d',sub) '_ses-meg1_task-Rest_channels.tsv']);
    D = spm_eeg_prep(S);
    D = chantype(D,indchantype(D,'megmag'),'MEGMAG');
    D = chantype(D,indchantype(D,'megplanar'),'MEGPLANAR');
    D.save
    
    % Downsampling the data
    S = [];
    S.D = D;
    S.method = 'resample';
    S.fsample_new = 500;
    D = spm_eeg_downsample(S);
    delete(S.D)
    
    % High-pass filter
    S = [];
    S.D = D;
    S.type = 'butterworth';
    S.band = 'high';
    S.freq = 0.5; % Cutoff frequency
    S.dir = 'twopass';
    S.order = 5;
    S.prefix = 'f';
    D = spm_eeg_filter(S);
    delete(S.D)
    
    
    % Low-pass filter
    S = [];
    S.D = D;
    S.type = 'butterworth';
    S.band = 'low';
    S.freq = 98; % Cutoff frequency
    S.dir = 'twopass';
    S.order = 5;
    S.prefix = 'f';
    D = spm_eeg_filter(S);
    delete(S.D)
    
end


%% PreProcess- Part 2 - OSL Artifacts detection

parfor sub = 1:nsub
    
    infile = fullfile(processed_pth,(sprintf('sub-Sub%04d',sub)),'ffdspmeeg');
    D = spm_eeg_load(infile);
    
    % OSL artifact detection
    D = osl_detect_artefacts(D,'modalities',unique(D.chantype(D.indchantype('MEGANY'))),'badchannels',false);
    D.save;
    
end

%% Rhino Co-Reg
anat = {'MRI'};
ana=1;
Nanat = length(anat);

UseHeadshape = 1; %load('fake_mesh') % latter just to speed up
UseRhino = 1;
RemoveNose = 0;

parfor sub=1:nsub
    
    
    infile = fullfile(processed_pth,(sprintf('sub-Sub%04d',sub)),'ffdspmeeg');
    D = spm_eeg_load(infile);
    if isfield(D,'inv')
        D = rmfield(D,'inv');
    end
    S0 = []; S0.D = D;
    
    S = [];
    S.useheadshape  = 1;
    S.use_rhino     = 1;
    S.fid.label.nasion  = 'Nasion';
    S.fid.label.lpa     = 'LPA';
    S.fid.label.rpa     = 'RPA';
    S.fid.coordsys      = 'Native';
    S.forward_meg = 'Single Shell';
    S.do_plots = 0;
    
    
    if RemoveNose
        fid = D.fiducials;
        %size(fid.pnt,1), figure,hold on; plot3(fid.pnt(:,1),fid.pnt(:,2),fid.pnt(:,3),'r.') % can comment out this line when running for all subjects - just to show you
        fid.pnt(find(fid.pnt(:,2)>0 & fid.pnt(:,3)<0),:)=[];
        %size(fid.pnt,1), plot3(fid.pnt(:,1),fid.pnt(:,2),fid.pnt(:,3),'bo'); rotate3d % can comment out this line when running for all subjects - just to show you
        D=fiducials(D,fid);
        D.save
    end
    megfid = D.fiducials;
    
    S.D = D;
    
    %        T1file = fullfile('/imaging/rh01/Projects/Defacing',anat{ana},sprintf('sub-Sub%04d_ses-meg1_T1w.nii',sub))
    T1file = fullfile('MRI',sprintf('sub-Sub%04d_ses-meg1_T1w.nii',sub));
    
    if exist(T1file,'file')
        
        fids = spm_jsonread(sprintf('BioFIND/MCIControls/sub-Sub%04d/ses-meg1/anat/sub-Sub%04d_ses-meg1_T1w.json',sub,sub));
        
        mrifid = []; mrifid.fid.label = {'Nasion';'LPA';'RPA'};
        V = spm_vol(T1file);
        
        nasion = V.mat*[fids.AnatomicalLandmarkCoordinates.Nasion; 1]; nasion = nasion(1:3)';
        lpa    = V.mat*[fids.AnatomicalLandmarkCoordinates.LPA; 1];    lpa    = lpa(1:3)';
        rpa    = V.mat*[fids.AnatomicalLandmarkCoordinates.RPA; 1];    rpa    = rpa(1:3)';
        mrifid.fid.pnt = [nasion; lpa; rpa];
        
        
        
        if UseHeadshape
            
            
            
            if UseRhino
                S.mri = T1file;
                S.fid.coords.nasion = nasion;
                S.fid.coords.lpa    = lpa;
                S.fid.coords.rpa    = rpa;
                if strcmp(anat{ana},'MRI') % ensure SPM mesh normalisation parameters not affected by defacing
                    S.refMRI = T1file;
                else
                    try, S = rmfield(S,'refMRI'); end
                end
                D = osl_headmodel(S);
                D.save;
                
                %D = rhino(S); %D.save not removed in this version
            else
                D.inv{1}.mesh = spm_eeg_inv_mesh(T1file,2);
                D.inv{1}.mesh.fid.fid.label{1} = 'Nasion'; D.inv{1}.mesh.fid.fid.label{2} = 'LPA'; D.inv{1}.mesh.fid.fid.label{3} = 'RPA'; % assumes spm_eeg_inv_mesh always returns in this order
                %if strcmp(anat{ana},'MRI') % if use same MRI, all SPM results should be same for all defacing
                %M = spm_eeg_inv_mesh(S.refMRI,2);
                
                S.template = 0;
                S.useheadshape = 1;
                S.sourcefid = megfid;
                S.targetfid = D.inv{1}.mesh.fid;
                M1 = spm_eeg_inv_datareg(S);
                %else
                %    D.inv{1}.datareg.fid_mri = mrifid;
                %    D.inv{1}.datareg.fid_meg = megfid;
                %end
                D.inv{1}.datareg(1).sensors = D.sensors('MEG');
                D.inv{1}.datareg(1).fid_eeg = S.sourcefid;
                D.inv{1}.datareg(1).fid_mri = ft_transform_headshape(inv(M1), S.targetfid);
                D.inv{1}.datareg(1).toMNI = D.inv{1}.mesh.Affine*M1;
                D.inv{1}.datareg(1).fromMNI = inv(D.inv{1}.datareg(1).toMNI);
                D.inv{1}.datareg(1).modality = 'MEG';
            end
            
            
            
        else
            
            try
                D.inv{1}.mesh = fake_mesh; % To speed up; doesn't matter for Fid Coreg below if incorrect mesh (unless want spm_eeg_inv_checkdatareg to be correct; just to get spm functions to work!)
            catch
                error('need a fake mesh')
                %D.inv{1}.mesh = spm_eeg_inv_mesh(T1file,2);
                %fake_mesh = D.inv{1}.mesh;
            end
            
            D.inv{1}.mesh.fid = rmfield(D.inv{1}.mesh.fid,'fid');
            for f = 1:3
                D.inv{1}.mesh.fid.fid.label{f} = mrifid.fid.label{f};
                D.inv{1}.mesh.fid.fid.pnt(f,:) = mrifid.fid.pnt(f,:);
            end
            %if strcmp(anat{ana},'MRI') % if use same MRI, all SPM results should be same for all defacing
            %M = spm_eeg_inv_mesh(S.refMRI,2);
            
            S.template = 0;
            S.useheadshape = 0;
            S.sourcefid = megfid;
            S.targetfid = D.inv{1}.mesh.fid;
            M1 = spm_eeg_inv_datareg(S);
            %else
            %    D.inv{1}.datareg.fid_mri = mrifid;
            %    D.inv{1}.datareg.fid_meg = megfid;
            %end
            D.inv{1}.datareg(1).sensors = D.sensors('MEG');
            D.inv{1}.datareg(1).fid_eeg = S.sourcefid;
            D.inv{1}.datareg(1).fid_mri = ft_transform_headshape(inv(M1), S.targetfid);
            D.inv{1}.datareg(1).toMNI = D.inv{1}.mesh.Affine*M1;
            D.inv{1}.datareg(1).fromMNI = inv(D.inv{1}.datareg(1).toMNI);
            D.inv{1}.datareg(1).modality = 'MEG';
            %D.save;  %spm_eeg_inv_checkdatareg(D)
            
            
            
        end
        
        %rhino_display(D)
        delete(sprintf('rhino*sub%04d*',sub))
        delete(sprintf('sub%04d*.nii',sub))
        delete(sprintf('sub%04d*.gii',sub))
        delete(sprintf('y_sub%04d*.nii',sub))
    end
    
end


%% Beamforming
modalities  = {{'MEG','MEGPLANAR'}};%, {'MEG'}, {'MEGPLANAR'}}; % Might be possible to store multiple montages per modality
p = parcellation(fullfile('Toolboxes/osl/parcellations/fMRI_parcellation_ds8mm.nii.gz'));

%mni_coords  = p.template_coordinates;
mni_coords = osl_mnimask2mnicoords(fullfile('Toolboxes/osl/std_masks/MNI152_T1_8mm_brain.nii.gz'));

parfor sub = 1:nsub
    
    T1file = fullfile('MRI',sprintf('sub-Sub%04d_ses-meg1_T1w.nii',sub));
    
    if exist(T1file,'file')
        
        infile = fullfile(processed_pth,(sprintf('sub-Sub%04d',sub)),'ffdspmeeg');
        D = spm_eeg_load(infile);
        D = D.montage('switch',0);
        D = osl_filter(D,[2 48]);
        
        SS = struct;
        SS.timespan          = [0 Inf];
        SS.pca_order         =  50;
        SS.type              = 'Scalar';
        SS.inverse_method    = 'beamform';
        SS.prefix            = '_';
        
        for m = 1:length(modalities)
            
            S = []; S.D = D;
            
            if length(modalities{m})>1
                SS.fuse = 'all';
                SS.modalities = modalities{m};
                SS.dirname = [infile '_bpPCA_248' strcat(modalities{m}{:}) '_BFn'];
                S.outfile = [infile '_bpPCA_248' strcat(modalities{m}{:})];
            else
                SS.fuse = 'no';
                SS.modalities = modalities{m}{1};
                SS.dirname = [infile '_' modalities{m}{1} '_BFn'];
                S.outfile = [infile '_' modalities{m}{1}];
            end
            
            DS = spm_eeg_copy(S);
            
          
            DS = osl_inverse_model(DS, mni_coords, SS); % coordis
            
            % Select montage
            DS = DS.montage('switch',1); % What difference between normalised and unnormalised weights?
            DS.save;
            
            DS = ROInets.get_node_tcs(DS,p.to_matrix(p.binarize),'pca');
           
            DS.save;
            
            
        end
    end
    
end

%% PreProcess- Part 3 - Despiking and calculate covariance matrix of ROIs

%freqbands = {[2 4],[4 8],[8 12],[12 30],[30 48],[2 48],[52 86]}
freqbands = {[2 48]};

modal = {'MEGMEGPLANAR','MEG', 'MEGPLANAR'};
nsub=324;

covariance = cell(1,numel(freqbands));
variance = cell(1,numel(freqbands));

for ii = 1:length(freqbands)
    
    Cov = cell(1,nsub); Var = Cov;
    
    parfor sub = 1:nsub
        
        T1file = fullfile('MRI',sprintf('sub-Sub%04d_ses-meg1_T1w.nii',sub));
        
        if exist(T1file,'file')
            
            infile = fullfile(processed_pth,(sprintf('sub-Sub%04d',sub)),'_ffdspmeeg_bpPCA2_98_MEGMEGPLANAR');
            D = spm_eeg_load(infile);
            D = D.montage('switch',3);
            
            % Remove bad badtrials
            chans = D.indchantype('LFP','GOOD'); % MEG :MAG or MEGPLANAR : GRD
            g = D(:,:);
            g = g(chans,good_samples(D,chans));

            g = ROInets.remove_source_leakage(g,'symmetric');
            
            % Filter to desired freq band
            y1 = ft_preproc_bandpassfilter(g, D.fsample, freqbands{ii}, 4, 'but');
            
            % Despiking
            y = filloutliers(y1,'clip','median','ThresholdFactor',3);
            
            % Extract envelops
            
            Hen_lc_sep = hilbenv(squeeze(y),1:D.nsamples,1,1);
            Hen_lc_sep = abs(Hen_lc_sep)
            
            
            %y1 = resample(Hen_lc_sep',20,D.fsample);
            
            % Calculate Covariance Matrix
            cm = cov(Hen_lc_sep);
            
            Cov{sub} = cm(find(triu(cm,1)))';
            Var{sub} = diag(cm)';
            
        end
    end
    
    Cov = cat(1,Cov{:});
    Var = cat(1,Var{:});
    covariance{ii} = Cov;
    variance{ii} = Var;
    
end

%% HMM inference

ii = 0;
data_new = cell(307,1);
T = cell(307,1);

for sub=1:324
    
    T1file = fullfile('MRI',sprintf('sub-Sub%04d_ses-meg1_T1w.nii',sub));
    
    if exist(T1file,'file')
        ii = ii + 1 ;
        infile = fullfile(processed_pth,(sprintf('sub-Sub%04d',sub)),'_ffdspmeeg_bpPCA2_98_MEGMEGPLANAR');
        
        D = spm_eeg_load(infile);
        D = D.montage('switch',3);
        chans = D.indchantype('LFP','GOOD'); % MEG :MAG or MEGPLANAR : GRD
        
        g = squeeze(D(:,:));
        
        g = g(chans,good_samples(D,chans));
        y1 = ft_preproc_bandpassfilter(g, D.fsample, [2 48], 4, 'but');
        
        y1 = ROInets.remove_source_leakage(y1,'symmetric');
        % Despiking
        y = filloutliers(y1,'clip','median','ThresholdFactor',3);
        
        % Extract envelops
        
        Hen_lc_sep = hilbenv(squeeze(y),1:D.nsamples,1,1);
        Hen_lc_sep = abs(Hen_lc_sep);
        
        % Smooth and normalise
        Hen_lc_sep = movmean(Hen_lc_sep,25,2,'omitnan'); % 100ms smoothing window
        for ll = 1:size(Hen_lc_sep,1) % num of parcels
            Hen_lc_sep(ll,:) = ( Hen_lc_sep(ll,:) - nanmean(Hen_lc_sep(ll,:)) ) ./ nanstd(Hen_lc_sep(ll,:));
        end
        
        
        data_new{ii} = Hen_lc_sep';
        T{ii} = size(Hen_lc_sep,2);
        
    end
end

data = transpose(data_new);
T = transpose(T);


% Prepare options structure
options = struct();
options.verbose = 1;

% These options specify the data and preprocessing that hmmmar might perform. Further options are discussed here
options.onpower = 0;
options.standardise = 0;
options.Fs = 500;

% Here we specify the HMM parameters
options.K = 8;  	  % The number of states to infer
options.order = 0; 	  % The lag used, this is only relevant when using MAR observations
options.zeromean = 0; 	  % We do want to model the mean, so zeromean is set off
options.covtype = 'full'; % We want to model the full covariance matrix
options.useParallel = 1;

options.BIGNinitbatch = 15;
options.BIGNbatch = 15;
options.BIGtol = 1e-7;
options.BIGcyc = 500;
options.BIGundertol_tostop = 5;
options.BIGdelay = 5;
options.BIGforgetrate = 0.7;
options.BIGbase_weights = 0.9;

% The following loop performs the main HMM inference. We start by
% estimating a 6 state HMM as used in the manuscript.
%states_to_infer = [8];

% Optionally, we can explore a wider range of values for K by looping through
% several values. This can be done by uncommenting the line below.
% Warning: This is likely to be extremely time-consuming to infer!

%states_to_infer = 2:2:12; % uncomment this line to explore different numbers of states
% The HMM inference is repeated a number of times and the results based on
% the iteration with the lowest free energy. Note that this can be
% extremely time-consuming for large datasets. For a quick exploration of
% results, nrepeats can be set to a smaller value or even 1. The full inference
% is run over 10 repeats.
%nrepeats = 1;

%for kk = states_to_infer
%   best_freeenergy = nan;
%   options.K = kk;

%    for irep = 1:nrepeats
% Run the HMM, note we only store a subset of the outputs
% more details can be found here: https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#estimation
[hmm, Gamma, ~, vpath, ~, ~, ~, ~, ~] = hmmmar (data',T',options);
%
%         if isnan(best_freeenergy) || fehist(end) < best_freeenergy
%             hmm = hmm_iter;
%             Gamma = Gamma_iter;
%             vpath = vpath_iter;
%         end
%     end
% Save the HMM outputs
%     hmm_outfile = fullfile( config.analysisdir, 'envelope_hmm', sprintf('envelope_HMM_K%d',options.K));
%     save( hmm_outfile ,'hmm','Gamma','vpath','T')
% end
%% HMM Dual for subject specific Gamma , HMM and states

clear hmm_s Gamma_s

parfor sub=1:nsub
    
    
    T1file = fullfile('MRI',sprintf('sub-Sub%04d_ses-meg1_T1w.nii',sub));
    
    if exist(T1file,'file')
        
        infile = fullfile(processed_pth,(sprintf('sub-Sub%04d',sub)),'_ffdspmeeg_bpPCA2_98_MEGMEGPLANAR');
        
        D = spm_eeg_load(infile);
        D = D.montage('switch',3);
        
        chans = D.indchantype('LFP','GOOD'); % MEG :MAG or MEGPLANAR : GRD
        
        g = squeeze(D(:,:));
        
        g = g(chans,good_samples(D,chans));
        y1 = ft_preproc_bandpassfilter(g, D.fsample, [2 48], 4, 'but');
        
        y1 = ROInets.remove_source_leakage(y1,'symmetric');
        % Despiking
        y = filloutliers(y1,'clip','median','ThresholdFactor',3);
        
        % Extract envelops
        
        Hen_lc_sep = hilbenv(squeeze(y),1:D.nsamples,1,1);
        Hen_lc_sep = abs(Hen_lc_sep);
        
        % Smooth and normalise
        Hen_lc_sep = movmean(Hen_lc_sep,25,2,'omitnan'); % 100ms smoothing window
        for ll = 1:size(Hen_lc_sep,1) % num of parcels
            Hen_lc_sep(ll,:) = ( Hen_lc_sep(ll,:) - nanmean(Hen_lc_sep(ll,:)) ) ./ nanstd(Hen_lc_sep(ll,:));
        end
        
        %D = squeeze(D(:,:,1));
        [hmm_s{sub},Gamma_s{sub},viterbi_s{sub}] = hmmdual(Hen_lc_sep',size(Hen_lc_sep,2),hmm);
    end
    
end

HMM_f = hmm_s(~cellfun('isempty',hmm_s));
Gamma_f = Gamma_s(~cellfun('isempty',Gamma_s));
%viterbi_s(~cellfun('isempty',viterbi_s));

%% Statistical Features

% Some useful information about the dynamics
for sub = 1:307
    
    
    FO{sub} = getFractionalOccupancy( Gamma_f{sub}, T{sub},options, 2);
    % Interval Time is the time between subsequent visits to a state
    IT = getStateIntervalTimes( Gamma_f{sub}, T{sub}, options,[]);
    ITmerged{sub} = cellfun(@mean,IT);clear IT
    % Life Times (or Dwell Times) is the duration of visits to a state
    LT = getStateLifeTimes( Gamma_f{sub}, T{sub}, []);
    LTmerged{sub} = cellfun(@mean,LT); clear LT
    
    SwitchingRate{sub} =  getSwitchingRate(Gamma_f{sub},T{sub},options); % rate of switching between stats
    
end

FO_mat = cat(1,FO{:});
IT_mat = cat(1,ITmerged{:});
LT_mat = cat(1,LTmerged{:});
SR_mat = cat(1,SwitchingRate{:});


SF = [FO_mat IT_mat LT_mat SR_mat];


%% Static FC using cov of states

clear Cov

for sub=1:307
    
    cm = [];
    cm = cov(Gamma_f{sub});
    %    Cov{sub} = (cm(:))';
    Cov{sub} = cm(find(triu(cm,1)))';
    %     Var{sub} = diag(cm)';
    
end

Cov = cat(1,Cov{:});

%% Static FC using cov of ROI envelops

clear Cov_itr2
parfor sub=1:324
    
    cm = [];
    % Extract envelops
    
    T1file = fullfile('MRI',sprintf('sub-Sub%04d_ses-meg1_T1w.nii',sub));
    
    if exist(T1file,'file')
        infile = fullfile(processed_pth,(sprintf('sub-Sub%04d',sub)),'_ffdspmeeg_bpPCA2_98_MEGMEGPLANAR');
        
        D = spm_eeg_load(infile);
        D = D.montage('switch',3);
        
        chans = D.indchantype('LFP','GOOD'); % MEG :MAG or MEGPLANAR : GRD
        
        g = squeeze(D(:,:));
        
        g = g(chans,good_samples(D,chans));
        y1 = ft_preproc_bandpassfilter(g, D.fsample, [2 48], 4, 'but');
        
        y1 = ROInets.remove_source_leakage(y1,'symmetric');
        % Despiking
        y = filloutliers(y1,'clip','median','ThresholdFactor',3);
        
        % Extract envelops
        
        Hen_lc_sep = hilbenv(squeeze(y),1:D.nsamples,1,1);
        Hen_lc_sep = abs(Hen_lc_sep);
        
        % Smooth and normalise
        Hen_lc_sep = movmean(Hen_lc_sep,25,2,'omitnan'); % 100ms smoothing window
        for ll = 1:size(Hen_lc_sep,1) % num of parcels
            Hen_lc_sep(ll,:) = ( Hen_lc_sep(ll,:) - nanmean(Hen_lc_sep(ll,:)) ) ./ nanstd(Hen_lc_sep(ll,:));
        end
        
        
        %y1 = resample(Hen_lc_sep',20,D.fsample);
        
        % Calculate Covariance Matrix
        cm = cov(Hen_lc_sep');
        
        Cov2{sub} = cm(find(triu(cm,1)))';
        %Var_itr2{sub} = diag(cm)';
        
    end
    
    
end

Cov2 = cat(1,Cov2{:});

%% Machine Learning 

% Needs MKL repository to be installed

clear V

V = {{Cov},{SF},{SF,Cov2}};
rng('default') % For reproducibility
[acc1,acc2] = mkl_ens(V,labels,'Hyper1',0.1,'Hyper2',1,...
    'Nfold',5,'Nrun',1000,'PCA_cut',0,'feat_norm',1,'ens',1);

mean(mean(acc2,3))
% ~64 without ortho and with downsampling
% ~61 with ortho and ds
% ~62  with ortho but no ds
% ~64.5  with neither ortho nor ds

titles = {'MEG(COV)','MEG(HMM-temporal feats)','MEG(COV),MEG(HMM-temporal feats)'};
pos_titles = {'MEG(HMM-temporal feats),MEG(COV)>MEG(COV)'};
%  define contrasts
c = [-1 0 1 ];

f1 = plot_results(titles,acc2,pos_titles,c); % main figure 1
sgtitle('Late Combination')
sgt.FontSize = 20;
