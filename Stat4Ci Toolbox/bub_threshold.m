function [Zpeak,Zextent]=bub_threshold(area,fwhm,bubres,pC,Ppeak,Pcluster,Pextent);

% Finds thresholds for peaks in scale-space bubbles Z images 
% using the saddle-point method of 
%   Chamandy, N., Worsley, K.J., Taylor, J.E. & Gosselin, F. (2007). Tilted
%   Euler characteristic densities for Central Limit random fields, with an 
%   application to 'bubbles'. Annals of Statistics, accepted.
% and the scale-space Lipschitz-Killing curvatures in 
%   Taylor, J.E. & Worsley, K.J. (2007). Random fields of multivariate test 
%   statistics, with applications to shape analysis. Annals of Statistics, 
%   in press.
%
% Usage: [Zpeak, Zextent]=
%   bub_threshold([area[,fwhm[,bubres[,pC[,Ppeak[,Pcluster[,Pextent]]]]]])
%
% area = either a single number, in which case the search region is
% assumed to be a 2D square, or a vector of D components giving the
% intrinsic volumes, e.g. for D=2: [ 1  perimeter/2  area]. 
% Default is [1 0 0] i.e. just a single point.
%
% fwhm = FWHM of smoothing kernel, in pixels. 
% If the area (above) is in resels, then set fwhm=1 (default). 
% If fwhm is a vector of two components [fwhm1 fwhm2] with fwhm1<=fwhm2, 
% then scale-space thresholds are calculated for FWHM in [fwhm1, fwhm2]. 
% 
% bubres = total number of bubbles per resel, i.e. #bubles per image x 
% #images / (area/fwhm^D) . Default is Inf, which gives Gaussian results.
% 
% pC = proportion of correct responses in the contrast of the images. 
% If pC is a vector it is the contrast values for each bubble image, e.g.
% 1/#correct repeated #correct times, then -1/#incorrect repeated
% #incorrect times, which is the same as the scalar
% pC=#correct/(#correct+#incorrect). If pC==1 then there is no subtraction,
% just an averaging of bubbles. Note that if pC is a vector, multplying all
% elements by a factor or doubling its length gives the same results, but
% the elements should sum to 0. Default is 0.5, i.e. guessing.
% 
% Ppeak: If the first component is <1 then Ppeak is treated as a vector of 
% P-values for peaks and Zpeak returns the corresponding thresholds. If the
% first component is >=1 then Ppeak is treated as peaks and Zpeak returns
% the P-value. Default is 0.05. 
% 
% Pcluster: If <1, then the threshold for clusters is chosen so that 
% P(Z>threshold)=Pcluster where Z~N(0,1). If Pcluster>=1 then Pcluster is
% the threshold. Default is 0.001. 
%     
% Pextent: If the first component is <1 then Pextent is treated as a vector
% of P-values for cluster extents and Zextent returns the corresponding
% thresholds, in the same units as fwhm. If the first component is >=1 then
% Pextent is treated as extents (in the same units as fwhm) and Zextent
% returns the P-value. Default is 0.05.  
% 
% Example: For a 240 x 380 image with bubble fwhm=24.9766, 750 images with
% 54 bubbles per image, pC=0.6827 correct responses, gives 
% 
% area = [1 240+380 240*380] 
% resels = 240*380/fwhm^2 = 146.1935
% bubres = 54*750/resels = 277.0300 
% 
%     bub_threshold(area,fwhm,bubres,pC)
% 
% Zpeak = 3.8121      P=0.05 threshold for peaks
% Pcluster = 0.001    P-value at a point used to define Zcluster
% Zcluster = 3.0319   threshold used to define clusters
% Zextent = 331.3720  P=0.05 threshold of cluster extent in pixels
% Pbound = 0.1396     probability that a cluster of size Zextent hits
%                     the boundary. This should be small, say <0.2, 
%                     otherwise cluster inference will be too conservative. 
%
% Example: For a scale space search for fwhm from 5 to 100 pixels:
% 
%     bub_threshold(area,[5 100],bubres,pC)
%
% Zpeak = 4.6227      (a higher threshold). 
% Pcluster = 0.001   
% Zcluster = 3.0319
% Zextent = 120.6872  (pixels in scale space)
% Pbound = 2.1315     (issues a warning:)
% Warning: Probability of a cluster intersecting the boundary  
% >0.2 so cluster inference is likely to be too conservative.

if nargin<1; area=[1 0 0]; end
if nargin<2; fwhm=1; end
if nargin<3; bubres=Inf; end
if nargin<4; pC=0.5; end
if nargin<5; Ppeak=0.05; end
if nargin<6; Pcluster=0.001; end
if nargin<7; Pextent=0.05; end

if length(area)==1
    area=[1 2*sqrt(area) area];
else
    area=area(1:max(find(area)));
end
if length(fwhm)==1
    fwhm=[fwhm fwhm];
end

D=length(area)-1;
fw=fwhm(1)/sqrt(4*log(2));
lkc=area./fw.^(0:D);

% For scale space:
scale=fwhm(2)/fwhm(1);
if scale<=1
    Ds=D;
    lkcs=lkc;
else
    Ds=D+1;
    lkcs=[lkc 0];
    scale=fwhm(2)/fwhm(1);
    Kappa=D/2;
    dlkc(1)=log(scale)*lkc(1);
    dlkc(2:(D+1))=(1-scale.^(-(1:D)))./(1:D).*lkc(2:(D+1));
    for d=1:Ds
        j=0:floor((D-d+1)/2);
        k=d+2*j-1;
        lkcs(d+1)=(1+scale^(-d))/2*lkcs(d+1) + sum( ...
            (-1).^j./(1-2*j).*exp(gammaln(k+1)-gammaln(j+1)-gammaln(d) ...
            +(1/2-j)*log(Kappa)-j*log(4*pi)).*dlkc(k+1));
    end
end

% EC densities
z=10:-0.01:-10;
nz=length(z);

theta=z;
if bubres==Inf
   I=z.^2/2;
   c2=ones(1,nz);
else
   % cumulants:
   nk=20;
   j=1:nk;
   if length(pC)==1
       if pC==1
           mu=ones(1,nk);
       else if pC==0
               mu=(-1).^j;
           else
               mu=pC.^(1-j)+(-1).^j.*(1-pC).^(1-j);
           end
       end
   else
      mu=zeros(1,nk);
      for jj=1:nk
         mu(jj)=mean(pC.^jj);
      end
   end
   kappa=mu./mu(2).^(j/2).*2.^(D*j/4)./j.^(D/2)./  ...
       ((pi/4/log(2))^(D/2)*bubres).^(j/2-1);
   % saddle point approximation:
   niter=5;
   for iter=1:niter
      c1=zeros(1,nz);
      c2=ones(1,nz);
      thetaj=ones(1,nz);
      for j=1:nk
         thetaj=thetaj.*theta/j; 
         if (j+1)<=nk; c1=c1+kappa(j+1).*thetaj; end;
         if (j+2)<=nk; c2=c2+kappa(j+2).*thetaj; end;
      end
      theta=theta-(c1-z)./c2;
   end
   c0=zeros(1,nz);
   cdot=ones(1,nz);
   thetaj=ones(1,nz);
   for j=1:nk
      thetaj=thetaj.*theta/j; 
      if j>=2; c0=c0+kappa(j).*thetaj; end;
   end
   I=theta.*c1-c0;
end

zt=sqrt(c2).*theta;
ec=zeros(Ds+1,nz,4);
% using \t\th and \lambda:
ec(1,:,1)=exp(zt.^2/2-I).*erfc(zt/sqrt(2))/2;
% using \t\th and \hat\lambda:
ec(1,:,2)=ec(1,:,1);
% using u^2/2 and \lambda:
ec(1,:,3)=exp(z.^2/2-I).*erfc(z/sqrt(2))/2;
% using u^2/2 and \hat\lambda:
ec(1,:,4)=ec(1,:,3);

if Ds>0
    %Hermite polynomials
    Hez=ones(Ds,nz);
    Hezt=ones(Ds,nz);
    if Ds>1
        Hez(2,:)=z;
        Hezt(2,:)=zt;
        for d=1:Ds-2
            Hez(d+2,:)=z.*Hez(d+1,:)-d*Hez(d,:);
            Hezt(d+2,:)=zt.*Hezt(d+1,:)-d*Hezt(d,:);
        end
    end
    f=((2*pi).^(-(2:(Ds+1))/2))'*exp(-I);
    c=exp(-(1:Ds)'*log(c2)/2);
    ec(2:(Ds+1),:,1)=f.*c.*Hezt;
    ec(2:(Ds+1),:,2)=f.*Hezt;
    ec(2:(Ds+1),:,3)=f.*c.*Hez;
    ec(2:(Ds+1),:,4)=f.*Hez;
end

peak=zeros(4,length(Ppeak));
for i=1:4
    ps=lkcs*ec(:,:,i);
    if Ppeak(1)<1
        peak(i,:)=minterp1(ps,z,Ppeak);
    else
        peak(i,:)=interp1(z,ps,Ppeak);
    end
end
Zpeak=min(peak)

% cluster
if Ds>0
    rho0=min(ec(1,:,:),[],3);
    if Pcluster<1
        Pcluster
        Zcluster=minterp1(rho0,z,Pcluster)
    else
        Zcluster=Pcluster
        Pcluster=interp1(z,rho0,Zcluster)
    end
    if Pcluster>0.01
        ['Warning: Pcluster>0.01 so cluster threshold is too low';
         'and so cluster inference is likely to be too liberal. ']
    end

    if Pextent(1)<1
        Zxs=zeros(4,length(Pextent));
        for i=1:4
            rho0=interp1(z,ec(1,:,i),Zcluster);
            rhoDs=interp1(z,ec(Ds+1,:,i),Zcluster);
            Zxs(i,:)=rho0/rhoDs/gamma(Ds/2+1)*...
                (-log(-log(1-Pextent)/lkcs(Ds+1)/rhoDs)).^(Ds/2);
        end
        Zx=min(Zxs);  
        Zextent=Zx*fw^Ds
    else
        Zx=Pextent/fw^Ds;
        Pxs=zeros(4,length(Zextent));
        for i=1:4
            Pxs(i,:)=1-exp(-lkcs(Ds+1)*rhoD*...
                exp(-(Zx*rhoDs*gamma(Ds/2+1)/rho0).^(2/Ds)));
        end
        Pextent=min(Pxs)
        Zextent=Pextent;
    end

    radius=(Zx*gamma(Ds/2+1)).^(1/Ds)/sqrt(pi);
    Pbound=lkcs(Ds)*2*radius/lkcs(Ds+1)
    if any(Pbound>0.2)
        ['Warning: Probability of a cluster intersecting the boundary'
         '>0.2 so cluster inference is likely to be too conservative.']
    end  
end

return

function iy=minterp1(x,y,ix);
% interpolates only the monotonically increasing values of x at ix
n=length(x);
mx=x(1);
my=y(1);
xx=x(1);
for i=2:n
   if x(i)>xx
      xx=x(i);
      mx=[mx xx];
      my=[my y(i)];
   end
end
iy=interp1(mx,my,ix,'linear',0); 
return
