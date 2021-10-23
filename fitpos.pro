; This code is based on parts of J. Crocker's feature3d.pro.
; See  J. C. Crocker and D. G. Grier, “Methods of digital video
; microscopy for colloidal studies,” J. Colloid Interface Sci. 179,
; 298–310 (1996).

; A first version of this code was first written by D. Reinke and U. Gasser.

function fitpos,image,pt,diameter,masscut = masscut, threshold = threshold, $
  info=info,silent=silent,nvoxel=nvoxel

if not keyword_set(masscut) then masscut = 0
if keyword_set(threshold) then $
	if threshold le 0 or threshold ge 0.9 then $
	  message,'Threshold value must be between 0.0 and 0.9!'

; make extents be the smallest odd integers bigger than diameter
if n_elements(diameter) eq 1 then diameter = fltarr(3)+diameter
extents = fix(diameter) + 1
extents = extents + (extents+1) mod 2
sz = size( image )
nx = sz(1)
ny = sz(2)
nz = sz(3)
if keyword_set(info) then print,'fitsphere(diameters)',extents

;	Put a border around the image to prevent mask out-of-bounds
a = fltarr( nx+extents(0), ny+extents(1), nz+extents(2) )
for i = 0,nz-1 do $	; do as a loop to reduce memory piggage
a(extents(0)/2:(extents(0)/2)+nx-1,extents(1)/2:(extents(1)/2)+ny-1, $
  extents(2)/2+i) = float( image(*,*,i) )
nx = nx + extents[0]
ny = ny + extents[1]
nz=  nz + extents[2]
x=pt[0,*]+extents[0]/2
y=pt[1,*]+extents[1]/2
z=pt[2,*]+extents[2]/2

; 	Set up some stuff....
nmax=n_elements(x)

xl = x - fix(extents(0)/2) 
xh = xl + extents(0) -1
yl = y - fix(extents(1)/2) 
yh = yl + extents(1) -1
zl = z - fix(extents(2)/2) 
zh = zl + extents(2) -1
m  = fltarr(nmax)
pd = fltarr(nmax)
thresh = fltarr(nmax)
nthresh = fltarr(nmax)
nvoxel = fltarr(nmax)

;	Set up some masks
rsq   = mydlrsqd3d(extents,yratio=diameter(1)/diameter(0),$
	zratio=diameter(2)/diameter(0))
mask  = rsq lt ((diameter(0)/2.))^2 +1.
nask  = total(mask)
rmask = (rsq+1./6.) * mask

imask = make_array( extents(0), extents(1), extents(2), /float, /index ) $
	mod ( extents(0) ) + 1.
xmask = mask * imask
imask = make_array( extents(1), extents(0), extents(2), /float, /index ) $
	mod ( extents(1) ) + 1.
ymask = mask * transpose(imask,[1,0,2])
imask = make_array( extents(2), extents(1), extents(0), /float, /index ) $
	mod ( extents(2) ) + 1.
zmask = mask * transpose(imask,[2,1,0])

;	Get the 'tops', i.e. peak heights
tops = a(x,y,z)
tops = tops(*)
if keyword_set(threshold) then thresh = tops*threshold

; 	if 'threshold' then get the max shell values
if keyword_set(threshold) then $
  for i=0L,nmax-1L do nthresh(i) = $
    total(((a(xl(i):xh(i),yl(i):yh(i), zl(i):zh(i))) * mask) gt thresh(i))/ $
    total(((a(xl(i):xh(i),yl(i):yh(i), zl(i):zh(i))) * mask) gt 0)

if keyword_set(nvoxel) then for i=0L,nmax-1L do nvoxel(i) = $
  total(((a(xl(i):xh(i),yl(i):yh(i), zl(i):zh(i))) * mask) gt 0)



;	Estimate the mass	
for i=0L,nmax-1L do m(i) = total( (a(xl(i):xh(i),yl(i):yh(i),$
	zl(i):zh(i))-thresh(i) > 0) * mask )

; do a masscut, and prevent divide by zeroes in the centroid calc.
w = where( m gt masscut, nmax )
if nmax eq 0 then begin
	message,'No features found!',/inf
	return,-1
endif
xl = xl(w)
xh = xh(w)
yl = yl(w)
yh = yh(w)
zl = zl(w)
zh = zh(w)
x = x(w)
y = y(w)
z = z(w)
m = m(w)
tops = tops(w)
thresh = thresh(w)
nthresh = nthresh(w)

if not keyword_set(silent) then $
  message, strcompress( nmax ) + ' features found.',/inf

;	Setup some result arrays
xc = fltarr(nmax)
yc = fltarr(nmax)
zc = fltarr(nmax)
rg = fltarr(nmax)

;	Calculate the radius of gyration^2
for i=0L,nmax-1L do rg(i) = total( (a(xl(i):xh(i),yl(i):yh(i),$
	zl(i):zh(i)) - thresh(i) > 0) * rmask )/m(i)

;	Calculate peak centroids
for i=0L,nmax-1L do begin
	xc(i) = total( (a(xl(i):xh(i),yl(i):yh(i),$
		zl(i):zh(i)) - thresh(i) >0) * xmask )
	yc(i) = total( (a(xl(i):xh(i),yl(i):yh(i),$
		zl(i):zh(i)) - thresh(i) >0) * ymask )
	zc(i) = total( (a(xl(i):xh(i),yl(i):yh(i),$
		zl(i):zh(i)) - thresh(i) >0) * zmask )
endfor

;	Correct for the 'offset' of the centroid masks
xc = xc / m - ((float(extents(0))+1.)/2.)
yc = yc / m - ((float(extents(1))+1.)/2.)
zc = zc / m - ((float(extents(2))+1.)/2.)

;	Update the positions and correct for the width of the 'border'
x = x + xc - extents(0)/2
y = y + yc - extents(1)/2 
z = z + zc - extents(2)/2 

if keyword_set(threshold) then $
return,[transpose(x),transpose(y),transpose(z),$
   transpose(m),transpose(rg),transpose(tops),transpose(nthresh)] else $
return,[transpose(x),transpose(y),transpose(z),$
	transpose(m),transpose(rg),transpose(tops)]

end
