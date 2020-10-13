function atmos_venus_vega2, z, g

; This subroutine calculates the ambient atmospheric
; temperature and pressure for the average Venus Atmosphere
; Pressure and Temperatures are given for every 1 km of altitude
; Pressures are given in Pa/10^5
; Vaues for T/P are taken from Linkin et al. (1987) "Thermal
; structure of the atmosphere of Venus from the results of 
; measurements taken by the landing vehicle VEGA-2", Cosmic
; Research, v25, no5, 501-512.

; Define atmospheric parameters

T0 = 734.5
P0 = 89.025*100000
zkm = z/1000
zi = fix(zkm)
deltazi = zkm-float(zi)
i = zi

Parray = [89.025,83.664,78.243,73.64,68.758,64.429,60.311,56.478,$
  52.486,48.979,45.459,42.405,39.505,36.74,34.165,31.538,29.303,$
  27.227,25.003,23.232,21.606,19.876,18.272,16.716,15.401,14.153,$
  12.923,11.916,10.788,9.864,8.973,8.249,7.55,6.755,6.119,5.611,$
  5.009,4.554,4.031,3.611,3.301,2.94,2.56,2.321,2.052,1.83,1.607,$
  1.459,1.278,1.136,0.972,0.863,0.745,0.66,0.554,0.474,0.417,0.342,$
  0.292,0.245,0.216,0.175,0.15,0.126]

Tarray = [734.5,726.5,716.8,708.6,698.8,689.9,681.7,673.7,665.2,658.6,$
  650.6,643.6,635.9,628.8,622.5,615.5,610.4,603.1,594.1,585.8,578.5,$
  570.3,561.7,553.1,545.3,536.8,528.0,520.3,510.8,502.4,493.6,486.3,$
  478.2,468.8,460.1,453.6,445.0,438.0,428.4,420.2,414.2,406.7,401.8,$
  397.3,391.6,387.4,380.7,372.8,363.6,355.6,344.2,335.9,325.3,316.2,$
  304.4,294.2,286.1,278.2,271.1,264.2,259.9,255.1,255.9,262.8]

; Calculate the atmospheric temperature and pressure gradients:

if (zkm gt 63) then begin
	if (z lt 63005) then print, ' Plume extends beyond atmos data.'
    P=Parray(63)
    Ta=Tarray(63)
endif else begin
    P1=Parray(i)
    T1=Tarray(i)
    if (deltazi eq 0) then begin
      P = P1*100000
      Ta = T1
    endif else begin
      P2 = Parray(i+1)
      T2 = Tarray(i+1)
      z1 = i
      z2 = (i+1)
      mP = (z2-z1)/(alog(P2)-alog(P1))
      lnP = (zkm-z1)/mP + alog(P1)
      P = (exp(lnP))*100000
      mT = (z2-z1)/(T2-T1)
      Ta = (zkm-z1)/mT + T1
    endelse
endelse

   return, [Ta,P]

end