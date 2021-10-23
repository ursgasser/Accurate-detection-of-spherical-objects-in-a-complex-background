; This code is based on parts of J. Crocker's feature3d.pro.
; See  J. C. Crocker and D. G. Grier, “Methods of digital video
; microscopy for colloidal studies,” J. Colloid Interface Sci. 179,
; 298–310 (1996).

; This version was first written by D. Reinke and U. Gasser.
; 30.1.2013 changes to save memory (U. Gasser)

function getlocalmax,pic,dia,sep,tops,hurdle=hurdle
;pic: the image-stack containing the features
;dia: size of diameter-mask used for excluding an area once a local
;     maximum is found.
;sep: size of compare-mask used for ...
;     1) searching for local maxima and 
;     2) for excluding voxels that cannot be a local maximum 
;        (in the case that no local maximum was found)
;     sep must be smaller or equal than dia !!
;     sep is increased to the next larger odd integer
;tops: (output) the values of the voxels containing local maxima
;hurdle: gives the percentage of dark voxels that are not considered
;     as potential local maxima (percentile=hurdle)
;     set percentile to 0.0 if you want every voxel to be a potential maximum
;     set it to 0.8 or so if you have lots of tiny spikes, to run faster

  if n_elements(hurdle) eq 0 then percentile=0.8 else percentile=hurdle 
  print,'getlocalmax: percentile=',percentile

  nx = long(n_elements(pic(*,0,0)))
  ny = long(n_elements(pic(0,*,0)))
  nz = long(n_elements(pic(0,0,*)))	
  ;set up masks
  extentd = fix(dia) + 1
  extentd = extentd + (extentd+1) mod 2
  extents = fix(sep) + 1
  extents = extents + (extents+1) mod 2
  print,'getlocalmax: diameter mask:',extentd
  print,'getlocalmax: separation mask:',extents
  drs=mydlrsqd3d(extentd,yratio=dia[1]/dia[0],zratio=dia[2]/dia[0])
  srs=mydlrsqd3d(extents,yratio=sep[1]/sep[0],zratio=sep[2]/sep[0])
  smask = srs lt (extents[0]/2)^2
  ddmask = drs gt (extentd[0]/2)^2
  ;Put padding around the image to prevent mask out-of-bounds
  image_ext = fltarr( nx+extentd(0), ny+extentd(1), nz+extentd(2) )
  for i = 0,nz-1 do $	; do as a loop to reduce memory piggage
    image_ext[extentd(0)/2:(extentd(0)/2)+nx-1, $
      extentd(1)/2:(extentd(1)/2)+ny-1, extentd(2)/2+i] = float( pic[*,*,i] )
  nx = nx + extentd[0]
  ny = ny + extentd[1]
  nz=  nz + extentd[2]
  nxy=nx*ny
  indexa=sort(image_ext)
  indexc=indexa[round(percentile*n_elements(indexa)):n_elements(indexa)-1]
  indexa = 0
  index=[reverse(indexc),0]; so it knows how to stop!
  indexc = 0
  allmin=image_ext[index[n_elements(index)-1]]
  idx=0L
  rr=index[idx]
  m=image_ext[rr]
  r = [-1L,-1L,-1L]  ;will contain the coords. of local maxima
  hash = bytarr(nx,ny,nz)+1B
  erwidx=n_elements(index)-1L
  z=rr/nxy & y=(rr-z*nxy)/nx & x=(rr-z*nxy-y*nx)
  repeat begin &$
    region = image_ext[x-extents[0]/2:x+extents[0]/2, $
      y-extents[1]/2:y+extents[1]/2, z-extents[2]/2:z+extents[2]/2]*smask &$
    wmax=where(region gt m,ngood) &$
    if ngood eq 0 then begin &$ ;it is a local maximum
      hhh=hash[x-extentd[0]/2:x+extentd[0]/2,y-extentd[1]/2:y+extentd[1]/2, $
        z-extentd[2]/2:z+extentd[2]/2] &$
      hash[x-extentd[0]/2:x+extentd[0]/2,y-extentd[1]/2:y+extentd[1]/2, $
        z-extentd[2]/2:z+extentd[2]/2]=hhh*ddmask &$
      r = [[r],[x,y,z]] &$
    endif else begin &$  ;it's not a local maximum
      ;print,'strange' &$
      w = where( region lt m and smask eq 1, nw) &$
      if nw gt 0 then begin &$
        ssmask = smask &$
        ssmask[w] = 0 &$
        hash[x-extents[0]/2:x+extents[0]/2,y-extents[1]/2:y+extents[1]/2, $
          z-extents[2]/2:z+extents[2]/2] = ssmask &$
      endif &$
    endelse &$
    repeat begin &$
      idx = idx+1
      rr=index[idx]
      z=rr/nxy & y=(rr-z*nxy)/nx & x=(rr-z*nxy-y*nx) &$
    endrep until(hash[x,y,z] eq 1 or idx ge erwidx) &$
    if (idx lt erwidx) then begin &$
      m = image_ext(x,y,z) &$
    endif else begin &$
      m = allmin &$
    endelse &$
  endrep until (m le allmin )
  r=r[*,1:n_elements(r[0,*])-1]
  print,'getlocalmax: ',n_elements(r[0,*]),'particles found'
  returnv=[r[0,*]-extentd[0]/2,r[1,*]-extentd[1]/2,r[2,*]-extentd[2]/2]
  tops=image_ext[r[0,*],r[1,*],r[2,*]]
  return, returnv
end
