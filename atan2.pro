; written October 2009 by Marissa Vogt

FUNCTION atan2, y, x

  z = double(atan(y/x))
  if n_elements(y) eq 1 then begin
    if y gt 0d and x lt 0d then z = z+!dpi
    if y lt 0d and x lt 0d then z = z+!dpi
  endif else begin
    if total(where(y gt 0d and x lt 0d)) ne -1 then z(where(y gt 0d and x lt 0d)) = z(where(y gt 0d and x lt 0d))+!dpi
    if total(where(y lt 0d and x lt 0d)) ne -1 then z(where(y lt 0d and x lt 0d)) = z(where(y lt 0d and x lt 0d))+!dpi
  endelse
  return, double(z)

END