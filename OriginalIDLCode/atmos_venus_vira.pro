function atmos_venus_vira, z, g

; This subroutine calculates the ambient atmospheric
; temperature and pressure for the average Venus Atmosphere.
; Pressure and Temperatures are given for every 1 km of altitude.
; Pressures are given in Pa/10^5.
; Vaues for T/P are taken from Seiff et al. (1985) "Models of
; the structure of the atmosphere of Venus from the surface to
; 100 km altitude", Advances in Space Research, v5 no11, 3-58.
; Values from 0 - 32 km are taken from Table 1-1, values from
; 33 - 100 mk are taken from Table 1-2 for 45 degrees latitude.

; Define atmospheric parameters

T0 = 735.3
P0 = 92.1*100000
zkm = z/1000
zi = fix(zkm)
deltazi = zkm-float(zi)
i = zi

Parray = [92.1,86.45,81.09,76.01,71.2,66.65,62.35,58.28,54.44,$
  50.81,47.39,44.16,41.12,38.26,35.57,33.04,30.66,28.43,26.33,$
  24.36,22.52,20.79,19.17,17.66,16.25,14.93,13.7,12.56,11.49,10.5,$
  9.581,8.729,7.94,7.211,6.537,5.916,5.345,4.82,4.338,3.897,3.495,$
  3.127,2.793,2.491,2.216,1.969,1.745,1.544,1.364,1.202,1.057,0.9258,$
  0.8087,0.7036,0.6095,0.5255,0.4505,0.3839,0.3249,0.2733,0.2289,$
  0.1908,0.1591,0.1321,0.1096,0.0906,0.07489,0.06169,0.05082,0.04174,$
  0.03428,0.02807,0.02298,0.01876,0.01531,0.01245,0.01012,0.00819,$
  0.006625,0.005324,0.004278,0.00341,0.002718,0.002149,0.001699,$
  0.001332,0.001044,0.000812,0.0006317,0.0004881,0.0003772,0.0002902,$
  0.0002233,0.0001715,0.0001317,0.0001011,0.00007768,0.00005974,$
  0.00004595,0.00003541,0.00002729]

Tarray = [735.3,727.7,720.2,712.4,704.6,696.8,688.8,681.1,673.6,$
  665.8,658.2,650.6,643.2,635.5,628.1,620.8,613.3,605.2,597.1,589.3,$
  580.7,572.4,564.3,556.0,547.5,539.2,530.7,522.3,513.8,505.6,496.9,$
  488.3,479.9,471.7,463.4,455.0,446.8,438.6,430.8,423.3,415.5,408.1,$
  401.6,395.0,388.3,381.8,376.1,370.2,364.6,357.5,349.7,341.1,332.5,$
  321.9,312.3,301.4,290.2,278.3,267.4,258.7,253.3,249.8,246.2,243.5,$
  240.7,238.3,235.8,233.9,231.9,230.1,228.2,226.4,224.6,222.8,221.0,$
  218.6,216.2,213.3,210.4,206.5,202.5,199.0,195.5,192.1,188.6,185.2,$
  181.8,179.2,176.6,174.5,172.3,171.4,170.4,170.0,169.5,169.8,170.0,$
  170.5,170.9,171.6,172.2]

; Calculate the atmospheric temperature and pressure gradients:

if (zkm gt 100) then begin
	if (z lt 100005) then print, ' Plume extends beyond atmos data.'
    P=Parray(100)
    Ta=Tarray(100)
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