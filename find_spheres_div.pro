; find spheres in 3d image using intensity gradients and the 'flow' of
; the gradient towards the center of the particle combined with the
; Crocker-Grier method. The gradient method finds the edge of the
; sphere while the Crocker-Grier method finds bright features.

; 17.11.2019: This version includes the suppression of 'bridges between
; particles using the Hessian matrix (matrix of second derivatives)
; and the eigenvalues and determinant of the Hessian matrix in each
; pixel.

; Urs Gasser
; 5.7.2019
; 17.11.2019 added Hessian matrix and suppression of 'bridges'
; 1.4.2020 calculate Hessian matrix using a faster method, see
;   https://en.wikipedia.org/wiki/Eigenvalue_algorithm#cite_note-2 for
;   eigenvalues of a 3x3 matrix.
; 7.5.2020 i) pad width reduced from object_size to object_size/2
;   ii) save memory: bg, bp as float instead of double
; 28.11.2020 improved kernel for gradient image with just one FFT.

function find_spheres_div, img, rad, object_size, pixel_size, $
  exclude_mask, compare_mask, fit_mask, threshold, $
  noise_size=noise_size, percentile=percentile, bp_div=bp_div, $
  kernels=kernels, fkbg=fkbg, suppress_bridges=suppress_bridges, $
  img_smooth=img_smooth, verbose=verbose
; img: 3d array with image
; rad: radius of spheres given in the units used in pixel_size
;      (e.g. microns or nm)
; object_size: object size parameters for bpass3d (Crocker-Grier method)
; pixel_size: x, y, and z size of a voxel
; exclude_mask: exclude mask for getlocalmax.pro
; compare_mask: compare mask for getlocalmax.pro
; fit_mask: fit mask for fitpos.pro
; threshold: threshold for fitpos.pro
; noise_size: size of noise to remove pixel noise, default is [1.,1.,1.]
; percentile: percentile parameter for getlocalmax.pro, default is 0.8
; bp_div: give a variable to return the bpass image and the gradient
;         flow image.
; kernels: -1: calculate kernels; variable given: take kernels from
;          the variable to save time
; fkbg: -1: calculate the kbg kernel in bpass3d and its Fourier transform.
;       The Fourier transform of kbg can be given in fkbg to save time.
; suppress_bridges: give a float for the strength of bride suppression
;       (default: 0.)
; img_smooth: give a variable to return the smoothed image (with
;       suppressed 'bridges' between particles)

  if not keyword_set(verbose) then verbose = 0
  if not keyword_set(noise_size) then noise_size = [1.,1.,1.]
  if n_elements(percentile) eq 0 then percentile = 0.8

  flag_make_ker = 1
  if keyword_set(kernels) then begin
    if typename(kernels) ne 'INT' then begin
      flag_make_ker = 0
      fga = kernels.fga
      fkg = kernels.fkg
    endif
  endif

  if not keyword_set(suppress_bridges) then supp_br = 0. $
  else supp_br = float(suppress_bridges)

  flag_calc_kbg = 1
  if keyword_set(fkbg) then begin
    if typename(fkbg) ne 'INT' then begin
      flag_calc_kbg = 0
    endif
  endif

  dims = size(img, /dimensions)

  ; add padding to the image
  padw = round( object_size/2 > 2.*noise_size )  ; width of padding
  ave = fltarr(dims[2])
  for i = 0, dims[2]-1 do ave[i] = total(img[*,*,i]) / product(dims[0:1])
  img_pad = fltarr(dims[0]+2*padw[0], dims[1]+2*padw[1], dims[2]+2*padw[2])
  img_pad[*,*,0:padw[2]-1] = ave[0]
  img_pad[*,*,dims[2]+padw[2]:*] = ave[dims[2]-1]
  for i = 0, dims[2]-1 do img_pad[*,*,i+padw[2]] = ave[i]
  img_pad[padw[0]:padw[0]+dims[0]-1, padw[1]:padw[1]+dims[1]-1, $
    padw[2]:padw[2]+dims[2]-1] = float(img)

  ; 3d arrays of x-, y-, and z-coordinates
  xx = fltarr(dims+2*padw)
  yy = fltarr(dims+2*padw)
  zz = fltarr(dims+2*padw)
  hx = (findgen(dims[0]+2*padw[0]) - (dims[0]/2+padw[0])) # $
    replicate(1., dims[1]+2*padw[1])
  hy = replicate(1., dims[0]+2*padw[0]) # $
    (findgen(dims[1]+2*padw[1]) - (dims[1]/2+padw[1]))
  for iz = 0,dims[2]+2*padw[2]-1 do begin
    xx[*,*,iz] = hx
    yy[*,*,iz] = hy
    zz[*,*,iz] = iz - (dims[2]/2+padw[2])
  endfor

  if flag_make_ker then begin
    if verbose then $
      print, 'find_spheres_div: calculating convolution kernels'
    ; Gaussian for smoothing the image
    ga = exp(-0.25*(xx^2/noise_size[0]^2 + yy^2/noise_size[1]^2 + $
      zz^2/noise_size[2]^2))
    ; sqrt(2) * noise_size is the standard deviation used in the Gaussian.
    ; This is also done in bpass3d_spherical.pro.
    ga /= total(ga)
    ga = shift(ga, -dims/2-padw)
    ; Fourier transform of Gaussian
    fga = fft(ga)
    ga = 0  ; free memory

    ; kernel for flow into a point
    ; phi for each pixel
    ph = atan(yy*pixel_size[1], xx*pixel_size[0])
    ; theta for each pixel
    th = -atan(zz*pixel_size[2], $
      sqrt((xx*pixel_size[0])^2 + (yy*pixel_size[1])^2)) + !pi/2

    ; unit vector towards center at each voxel
    ks = fltarr([dims+2*padw, 3])
    sth = sin(th)
    ks[*,*,*,0] = -sth * cos(ph)
    ks[*,*,*,1] = -sth * sin(ph)
    ks[*,*,*,2] = -cos(th)
    sth = 0 & ph = 0 & th = 0  ; free memory
    ; calculate distance^2 for each pixel
    rr = (xx*pixel_size[0])^2 + (yy*pixel_size[1])^2 + (zz*pixel_size[2])^2
    w = where(rr gt (1.1*rad)^2)
    rr = 0  ; free memory
    ks[w] = 0.
    ks[w+long(product(dims+2*padw))] = 0.
    ks[w+2*long(product(dims+2*padw))] = 0.
    h = dims/2+padw
    ks[h[0], h[1], h[2], [0,1,2]] = 0.  ; zero vector at origin
    for i = 0,2 do ks[*,*,*,i] = shift(ks[*,*,*,i], -dims/2-padw)
    ; Fourier transform of kernel ks
    fks = dcomplexarr([dims+2*padw,3])
    for i = 0,2 do fks[*,*,*,i] = fft(ks[*,*,*,i])
    ks = 0  ; free memory

    fkg = complex(0.,shift(xx, -dims/2-padw)) * fks[*,*,*,0] + $
      complex(0.,shift(yy, -dims/2-padw)) * fks[*,*,*,1] + $
      complex(0.,shift(zz, -dims/2-padw)) * fks[*,*,*,2]
    fks = 0  ; free memory

    if arg_present(kernels) then kernels = {fga: fga, fkg: fkg}
  endif

  ; FT of smoothed image
  fimg = fft(img_pad)
  fimg_s = fimg * fga

  if arg_present(img_smooth) then begin
    img_s = float(fft(fimg_s, /inverse))
    img_smooth = img_s[padw[0]:dims[0]+padw[0]-1, padw[1]:dims[1]+padw[1]-1, $
      padw[2]:dims[2]+padw[2]-1]
  endif

  if supp_br gt 0. then begin
    if verbose then print,'find_spheres_div: with Hessian matrix'
    ; calculate Hessian matrix of smoothed image
    img_he = dblarr([dims+2*padw, 6])
    img_he[*,*,*,0] = fft( -fimg_s* shift(xx^2, -dims/2-padw), /inverse)
    img_he[*,*,*,1] = fft( -fimg_s* shift(yy^2, -dims/2-padw), /inverse)
    img_he[*,*,*,2] = fft( -fimg_s* shift(zz^2, -dims/2-padw), /inverse)
    img_he[*,*,*,3] = fft( -fimg_s* shift(yy, -dims/2-padw)* $
      shift(zz, -dims/2-padw), /inverse)
    img_he[*,*,*,4] = fft( -fimg_s* shift(zz, -dims/2-padw)* $
      shift(xx, -dims/2-padw), /inverse)
    img_he[*,*,*,5] = fft( -fimg_s* shift(xx, -dims/2-padw)* $
      shift(yy, -dims/2-padw), /inverse)

    ; get eigenvalues and eigenvectors of Hessian matrix
    p1 = total(img_he[*,*,*,3:5]^2, 4)  ; A_xy^2 + A_xz^2 + A_yz^2
    w0 = where(p1 eq 0, nw0, complement=w0compl)
    eval0 = dblarr(dims+2*padw)
    eval1 = dblarr(dims+2*padw)
    eval2 = dblarr(dims+2*padw)
    if nw0 gt 0 then begin
      ; cases with diagoal matrix
      eval0[w0] = (img_he[*,*,*,0])[w0]
      eval1[w0] = (img_he[*,*,*,1])[w0]
      eval2[w0] = (img_he[*,*,*,2])[w0]
    endif
    if n_elements(w0compl) gt 0 then begin
      ; Hessian matrix is not diagonal
      ; q = trace/3
      q = total(img_he[*,*,*,0:2], 4) / 3.
      q = q[w0compl]
      h0 = (img_he[*,*,*,0])[w0compl]
      h1 = (img_he[*,*,*,1])[w0compl]
      h2 = (img_he[*,*,*,2])[w0compl]
      h3 = (img_he[*,*,*,3])[w0compl]
      h4 = (img_he[*,*,*,4])[w0compl]
      h5 = (img_he[*,*,*,5])[w0compl]
      p = (h0 - q)^2 + (h1 - q)^2 + (h2 - q)^2 + 2.*p1
      p = sqrt(p/6.)
      ; transform h0...5 to B = (1/p) * (A - q*I)
      h0 -= q  &  h0 /= p
      h1 -= q  &  h1 /= p
      h2 -= q  &  h2 /= p
      h3 /= p
      h4 /= p
      h5 /= p
      ; det(B) / 2
      r = ( h0 * (h1*h2 - h3^2) - h5 * (h5*h2 - h4*h3) + $
        h4 * (h5*h3 - h4*h1) ) / 2.
      ; in exact arithmetic for a symmetric matrix -1 <= r <= 1
      ; but computation error can leave it slightly outside this range.
      wm = where(r le -1., nwm)
      wp = where(r ge 1., nwp)
      phi = acos(r) / 3.
      if nwm gt 0 then phi[wm] = !pi/3.
      if nwp gt 0 then phi[wp] = 0.
      ; eigenvalues with eig3 <= eig2 <= eig1
      eval0[w0compl] = q + 2.*p * cos(phi)
      eval2[w0compl] = q + 2.*p * cos(phi + (2.*!pi/3.))
      ; trace = eig1 + eig2 + eig3
      eval1[w0compl] = 3.*q - eval0[w0compl] - eval2[w0compl]
      h0 = 0 & h1 = 0 & h2 = 0 & h3 = 0 & h4 = 0 & h5 = 0 ; free memory
      q = 0 & p = 0 & r = 0
    endif
    img_he = 0  ; free memory
    p1 = 0
    ; determinant
    determ = eval0 * eval1 * eval2

    ; factor to weaken pixels that are not in a 'maximum area'
    w = where(eval0 lt 0 and eval1 lt 0 and eval2 lt 0, nw, compl=cw)
    eval0 = 0  ; free memory
    eval1 = 0
    eval2 = 0
    factor_img = dblarr(size(img_pad, /dim))
    factor_img[cw] = abs(determ[cw])
    factor_img *= -supp_br / max(factor_img)
    factor_img = float(1. + factor_img)
    ; correct smoothed image
    img_s = float(fft(fimg_s, /inverse))
    img_s *= factor_img
    factor_img = 0  ; free memory
    fimg_s = fft(img_s)
  endif else begin
    if verbose then print, 'find_spheres_div: without Hessian matrix'
    if not arg_present(img_smooth) then img_s = float(fft(fimg_s, /inverse))
  endelse

  img_div = -double(fft(fimg_s * fkg, /inverse)) > 0
  img_div *= 255./max(img_div)


  ; bpass image with spherical convolution kernel (Crocker-Grier method)
  if flag_calc_kbg then begin
    if verbose then $
      print, 'find_spheres_div: calculating kernel for bpass background'
    ; kernel for background
    kbg = fltarr(dims+2*padw)
    hy = yy * (float(2*padw[0]+1) / float(2*padw[1]+1))
    hz = zz * (float(2*padw[0]+1) / float(2*padw[2]+1))
    rr2 = xx^2 + hy^2 + hz^2
    w = where(rr2 lt (0.5*object_size[0]+1)^2, nw)
    kbg[w] = 1./nw
    kbg = shift(kbg, -dims/2-padw)
    ; calculate the background image
    fkbg = fft(kbg)
    kbg = 0  ; free memory
  endif

  ; background image brought to intensity comparable to img_s
  bg = float(fft(fimg_s*fkbg, /inverse)) * float(product(dims+2*padw))

  ; the enhanced image (bpass image)
  bp = (img_s-bg) > 0.
  bp *= 255./max(bp)

  ; combine images and find the feature coordinates
  img_combined = img_div*bp
  ; remove padding
  img_combined_nopad = img_combined[padw[0]:dims[0]+padw[0]-1, $
    padw[1]:dims[1]+padw[1]-1, padw[2]:dims[2]+padw[2]-1]
  img_combined_nopad *= 255./max(img_combined)
  pos = getlocalmax(img_combined_nopad, exclude_mask, compare_mask, tops, $
    hurdle=percentile)
  ;optimize positions
  coo = fitpos(img_combined_nopad, pos, fit_mask, threshold=threshold)

  if arg_present(bp_div) then begin
    bp_div = {bp: bp[padw[0]:dims[0]+padw[0]-1, padw[1]:dims[1]+padw[1]-1, $
      padw[2]:dims[2]+padw[2]-1], div: img_div[padw[0]:dims[0]+padw[0]-1, $
      padw[1]:dims[1]+padw[1]-1, padw[2]:dims[2]+padw[2]-1], $
      bg: bg[padw[0]:dims[0]+padw[0]-1, padw[1]:dims[1]+padw[1]-1, $
      padw[2]:dims[2]+padw[2]-1], i2: img_combined_nopad}
  endif
  return, coo
end
