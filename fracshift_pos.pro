; get feature positions with fracshift
; The feature positions are fitted iteratively by shifting the image by
; a fraction of a pixel to the actual position.

; U. Gasser
; 19.9.2019
; following the OpticsExpress paper of Gao and Kilfoil 2009

function fracshift_pos, image, posi, diameter, nrep, $
  verbose=verbose, allpos=allpos
; image: the image containing all features
; posi: approximate positions of all features
; diameter: expected diameter of the features (in pixel units)
; nrep: number of frac-shift iterations
; verbose: set this keyword for some output
; allpos: give a variable to return the position after each iteration
;         of fracshift.
; output: [x, y, z, brightness, Rg^2, peak intensity,
;          dr (displacement due to fracshift)]

  if n_elements(diameter) eq 1 then diameter = fltarr(3)+diameter
  ; make extents be the smallest odd integers larger than diameter
  extents = round(diameter)
  extents = extents + (extents+1) mod 2
  if keyword_set(verbose) then print, 'fitsphere (diameters):', extents

  rsq = mydlrsqd3d(extents, yratio=diameter[1]/diameter[0],$
    zratio=diameter[2]/diameter[0])
  mask  = rsq lt (diameter[0]/2.)^2 +1.
  rmask = (rsq+1./6.)*mask

  ; x-, y-, and z-masks to calculate corrections of feature positions
  imask = make_array( extents, /float, /index ) mod ( extents[0] ) + 1.
  xmask = mask * imask
  imask = make_array( extents[1], extents[0], extents[2], /float, /index ) $
    mod ( extents[1] ) + 1.
  ymask = mask * transpose(imask,[1,0,2])
  imask = make_array( extents[2], extents[1], extents[0], /float, /index ) $
    mod ( extents[2] ) + 1.
  zmask = mask * transpose(imask,[2,1,0])

  ; add padding around the image
  dim = size(image, /dimensions)
  nx = dim[0] & ny = dim[1] & nz = dim[2]
  nxp = nx+extents[0]-1 & nyp = ny+extents[1]-1 & nzp = nz+extents[2]-1
  img = fltarr( nx+extents[0], ny+extents[1], nz+extents[2] )
  ave = fltarr(dim[2])
  for i = 0, dim[2]-1 do ave[i] = total(image[*,*,i]) / product(dim[0:1])
  img[*,*,0:extents[2]/2-1] = ave[0]
  img[*,*,dim[2]+extents[2]/2:*] = ave[dim[2]-1]
  for i = 0, dim[2]-1 do img[*,*,i+extents[2]/2] = ave[i]
  for i = 0L,nz-1 do $  ; do as a loop to save memory
    img[extents[0]/2:(extents[0]/2)+nx-1, extents[1]/2:(extents[1]/2)+ny-1, $
    extents[2]/2+i] = float( image[*,*,i] )

  ; positions in padded image
  pos = posi[0:2,*]
  for i = 0,2 do pos[i,*] += extents[i]/2
  npos = n_elements(pos[0,*])

  ; get peak heights for all features
  tops = img[pos[0,*], pos[1,*], pos[2,*]]

  corr = fltarr(3,npos)  ; the corrections of all positions
  bri = fltarr(npos)  ; the total brightness of all features
  rg = fltarr(npos)  ; radius of gyration^2 for all features
  frac_img = fltarr( extents )  ; array for frac-images

  if arg_present(allpos) then $
    ; will contain all positions for all iterations
    allpos = fltarr(3,npos,nrep)

  ; loop for repeated centroid fitting
  for ic = 0L, nrep-1 do begin &$
    if arg_present(allpos) then allpos[*,*,ic] = pos &$
    frac_pos = pos[0:2,*] mod 1. &$
    ; lower and upper limits to get shifted images for all features
    pl = long(pos[0:2,*] - (extents/2)#replicate(1.,npos)) &$
    pu = pl + extents#replicate(1,npos) - 1 &$
    ; loop over all features
    for ind = 0L, npos-1 do begin &$
      ll = pl[*,ind] &$  ; lower bounds for fracshift image
      ul = pu[*,ind] &$  ; upper bounds for fracshift image
      if min(ll) ge 0 and ul[0] lt nxp and ul[1] lt nyp and $
        ul[2] lt nzp then begin &$
        frac = [[1.-frac_pos[*,ind]], [frac_pos[*,ind]]] &$
        frac_img[*] = 0. &$  ; start with empty fracshift image
        ; loops ix, iy, iz for pixel at position and the pixel next to it
        for ix = 0,1 do for iy = 0,1 do for iz = 0,1 do begin &$
          frac_img += frac[0,ix] * frac[1,iy] * frac[2,iz] * $
            img[ll[0]+ix:ul[0]+ix, ll[1]+iy:ul[1]+iy, ll[2]+iz:ul[2]+iz] &$
        endfor &$
        ; update the total brightness of the feature
        bri[ind] = total( frac_img * mask ) &$
        ; calculate correction of pos
        corr[0,ind] = total(frac_img * xmask) &$
        corr[1,ind] = total(frac_img * ymask) &$
        corr[2,ind] = total(frac_img * zmask) &$
        ; calculate Rg^2 in last iteration
        if ic eq nrep-1 then rg[ind] = total(frac_img * rmask) / bri[ind] &$
      endif else $
        corr[*,ind] = (extents+1.)/2.*bri[ind] &$  ; leave this pos as it is
    endfor &$
    ; normalize corrections with brightness and subtract offset
    corr /= ([1.,1,1]#bri) &$
    corr -= ((extents+1.)/2.) # replicate(1.,npos) &$
    pos += corr &$             ; correct the positions
     if keyword_set(verbose) then $
       print, ic, total(corr[0,*])/npos, total(corr[1,*])/npos, $
       total(corr[2,*])/npos
  endfor

  ; the corrected positions
  h = (extents/2)#replicate(1.,npos)
  pos -= h
  if arg_present(allpos) then $
    for i = 0,nrep-1 do allpos[*,*,i] -= h

  ; effect of position fitting
  diff = pos[0:2,*] - posi[0:2,*]
  d2 = total(diff^2, 1)

  return, [pos, transpose(bri), transpose(rg), tops, transpose(d2)]
end
