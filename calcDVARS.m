function DVARS = calcDVARS(fMRI, gm_mask)

X = size(fMRI,1); Y = size(fMRI,2); Z = size(fMRI,3); T = size(fMRI,4);
Ydata = reshape(fMRI,[X*Y*Z,T]); 
Mdata = reshape(gm_mask,[X*Y*Z,1]);

Ydata = Ydata(Mdata>0,:);

I1   = size(Ydata,1); 

DY    = diff(Ydata,1,2);
DVARS = sqrt(sum(DY.^2)./I1);
std(DVARS);