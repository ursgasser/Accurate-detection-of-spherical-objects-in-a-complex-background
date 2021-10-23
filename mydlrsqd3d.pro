;
;	produce a 3d, anisotropic parabolic mask
;	anisotropic masks are 'referenced' to the x-axis scale.
;	ratios less than one squash the thing relative to 'x'
;	using float ratios allows precise 'diameter' settings 
;	in an odd-sized mask.
;
; create mask for mydlocalmax.pro
; mydlrsqd3d.pro is the same function as lrsqd3d.pro in feature3d.pro

function mydlrsqd3d,extent,yratio=yratio,zratio=zratio

if not keyword_set(yratio) then yratio = 1.
if not keyword_set(zratio) then zratio = 1.
if n_elements(extent) eq 1 then ext = intarr(3)+extent else ext = extent
x = ext(0)
y = ext(1)
z = ext(2)

r2 = fltarr(x,y,z,/nozero)
xc = float(x-1) / 2.
yc = float(y-1) / 2.
zc = float(z-1) / 2.

yi = fltarr(x) +1
xi = fltarr(y) +1

xa = (findgen(x) - xc)
xa = xa^2
ya = (findgen(y) - yc)/yratio
ya = ya^2
za = (findgen(z) - zc)/zratio
za = za^2

for k=0,z-1 do begin
	r2(*,*,k) = (xi ## xa) + (ya ## yi) + za(k) 
endfor

return,r2
end
