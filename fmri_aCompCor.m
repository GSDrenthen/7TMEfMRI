function fmri_aCompCor(func_files,rp_file,seg_file,brain_mask_file,ALVIN_mask)
% code is adapted from the CNIR-fmri toolbox (https://github.com/KKI-CNIR/CNIR-fmri_preproc_toolbox)

nm=load(rp_file);
nm=detrend(nm,'constant');
nui.names={'mp_x','mp_y','mp_z','mp_yaw','mp_pitch','mp_roll'};
nui.tc=nm;
nui.tc=[nui.tc, [0,0,0,0,0,0; diff(nui.tc)]];
n=1;
for i=(length(nui.names)+1):(length(nui.names)+6)
	nui.names{i}=['d',nui.names{n}];
    n=n+1;
end
nui.tc=[nui.tc,nm.^2];
n=1;
for i=(length(nui.names)+1):(length(nui.names)+6)
	nui.names{i}=['s',nui.names{n}];
    n=n+1; 
end

Vm = spm_vol(brain_mask_file);
brain_mask = spm_read_vols(Vm);
brain_mask=logical(brain_mask);

Vm = spm_vol(seg_file);
seg_mask = spm_read_vols(Vm);
wm_mask = (seg_mask==2) + (seg_mask==41);
wm_mask = double(imerode(logical(wm_mask),ones(3,3,3)));
csf_mask = (seg_mask==4) + (seg_mask==43);
csf_mask = double((logical(csf_mask)));

alvin_mask = spm_read_vols(spm_vol(ALVIN_mask));
%multiple by ALVIN mask of the ventricles
csf_mask = (csf_mask .* alvin_mask) > 0;
csf_mask = double(imerode(logical(csf_mask),ones(1,1,1))); %Erode

V = spm_vol(func_files);
Y = spm_read_vols(V);
Y = brain_mask.*Y;
wm_mask = (wm_mask .* brain_mask) > 0;
csf_mask = csf_mask .* brain_mask;

tc_wm=zeros(size(Y,4),sum(wm_mask(:)));
wm_mask=logical(wm_mask);

tc_csf=zeros(size(Y,4),sum(csf_mask(:)));
csf_mask=logical(csf_mask);

mY=zeros(size(Y,4),1);

for i=1:size(Y,4)
    tY=Y(:,:,:,i);
    tY(isnan(tY))=0;
    mY(i) = mean(tY(brain_mask))'; 
    tc_wm(i,:)=tY(wm_mask);
    tc_csf(i,:)=tY(csf_mask);
end
mY=detrend(mY,'constant');%demean
clear tY;

tc_wm = bsxfun(@minus, tc_wm, mean(tc_wm, 1));
        
npca_wm=mean(tc_wm,2);
[u,s] = svd(cov(tc_wm'));
eigen_values = diag(s);
ncomp = 5;
npca_wm = u(:,1:ncomp);
nui.wm_percent = sum(eigen_values(1:ncomp))/sum(eigen_values);
clear u s ncomp;

%concatenate white matter timeseries to variable nui
nui.tc=[nui.tc,npca_wm];
for i=1:size(npca_wm,2)
    nui.names{1+length(nui.names)}=['wm',num2str(i)];
end

npca_csf=mean(tc_csf,2);
[u,s] = svd(cov(tc_csf'));
eigen_values = diag(s);
ncomp = 5;
npca_csf = u(:,1:ncomp);
nui.csf_percent = sum(eigen_values(1:ncomp))/sum(eigen_values);
clear u s ncomp

%concatenate csf timeseries to variable nui
nui.tc=[nui.tc,npca_csf];
for i=1:size(npca_csf,2)
    nui.names{1+length(nui.names)}=['csf',num2str(i)];
end

nui.tc(:,i)=nui.tc(:,i)./max(squeeze(nui.tc(:,i)));

nt = size(nui.tc,1);
ntc = [nui.tc,ones(nt,1), (0:1/nt:(1-1/nt))'];% nuisances, constant, linear
precal_inv = ntc*((ntc'*ntc)\ntc');

W = single(reshape(spm_read_vols(V), V(1).dim(1)*V(1).dim(2)*V(1).dim(3), size(V,1)));
mY_time = reshape(mean(W,2), V(1).dim(1), V(1).dim(2), V(1).dim(3));

W=transpose(W);
Y = W-precal_inv*W;
Y(:, logical(~brain_mask)) = 0;
Y = permute(reshape(Y, size(V,1), V(1).dim(1), V(1).dim(2), V(1).dim(3)), [2, 3, 4, 1]);

for i_time = 1:size(V,1)
    Vo=V(i_time);
    [pathname, filename, ext] = fileparts(V(i_time).fname);
    %%if GSR included
    npref = 'aCompCor-';
    Vo.fname = fullfile(pathname,[npref,filename,ext]);
    
    Vo.private.dat.fname = Vo.fname;
    Y(:,:,:,i_time) = Y(:,:,:,i_time) + mY_time;
    spm_write_vol(Vo, Y(:, :, :, i_time));
end
