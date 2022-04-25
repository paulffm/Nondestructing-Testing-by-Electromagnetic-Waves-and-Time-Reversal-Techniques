function [exedg,ghostedg]= getGhostEdges(nx,ny,nz)
% getGhostEdges(msh) returns the canonical indices of all ghost edges in
% the considered mesh
%   
% Input:
%   msh             msh struct as created by cartMesh
%
% Output:
%   idxGhostEdges   column vector that contains all ghost edges' indices
%
% See also cartMesh

Mx = 1;
My = nx;
Mz = nx*ny;
np=nx*ny*nz;

% calculates indices for boundary xmax
indy = repmat(1:ny,1,nz);
indz = reshape(repmat(1:nz,ny,1),1,ny*nz);
n_xmax = 1+(nx-1)*Mx+(indy-1)*My+(indz-1)*Mz;

% calculates indices for boundary ymax
indx = repmat(1:nx,1,nz);
indz = reshape(repmat(1:nz,nx,1),1,nx*nz);
n_ymax = 1+(indx-1)*Mx+(ny-1)*My+(indz-1)*Mz;

% calculates indices for boundary zmax
indx = repmat(1:nx,1,ny);
indy = reshape(repmat(1:ny,nx,1),1,nx*ny);
n_zmax = 1+(indx-1)*Mx+(indy-1)*My+(nz-1)*Mz;

% return all indices in one column vector
ghostedg = [n_xmax,np+n_ymax,2*np+n_zmax]';
alledg=1:1:3*np;
exedg=setdiff(alledg,ghostedg);
