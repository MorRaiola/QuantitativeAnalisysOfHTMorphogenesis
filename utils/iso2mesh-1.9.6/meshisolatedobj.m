function [node,elem]=meshisolatedobj(img,opt,maxvol)
%
% Usage:
%    [node,elem]=meshisolatedobj(img,opt,maxvol)
%    Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%
img=logical(img);
[bin,regionnum]=bwlabeln(img,6);
node=[];
elem=[];

for regionid=1:regionnum
    fprintf(1,'meshing region #%d ...\n',regionid);

    [no,el]=v2m(bin==regionid,0.5,opt,maxvol);

    % merge the resulting mesh with other regions
    el(:,5)=regionid;
    el(:,1:4)=el(:,1:4)+size(node,1);
    if(regionid>1)
        elem=[elem;el];
        node=[node;no];
    else
        node=no;
        elem=el;
    end
end