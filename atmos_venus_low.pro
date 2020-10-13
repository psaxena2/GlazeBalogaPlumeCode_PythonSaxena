function atmos_venus_low, z, g

; This subroutine calculates the ambient atmospheric
; temperature and pressure for the Venus Atmosphere at low latitudes
; Pressure and Temperatures are given for every 2 km of altitude
; 
; Define atmospheric parameters

T0 = 733.0
P0 = 92.1*100000
zkm = z/1000
zi = fix(zkm)
deltazi = zkm-2.*float(zi/2)
i = zi/2

Parray = [92.100,81.0900,71.2000,62.3500,54.400,47.300,41.00,35.400,$
	30.500,26.200,22.400,19.00,16.100,13.5500,11.3500,9.450,$
	7.8200,6.430,5.250,4.250,3.42,2.730,2.170,1.71,1.336,1.034,0.791,$
	0.595, 0.440,0.318,0.227,0.159,0.110,0.076,0.0518,0.0352,0.0238,0.0159,$
	0.01054,0.00691,0.00447,0.00284,0.00176,0.00107,0.000637,0.000374,$
	0.000217,0.000125,0.0000719,0.0000415,0.0000242]

Tarray = [733.0,717.0,701.0,685.0,669.0,653.0,637.0,623.0,609.0,595.0,579.0,$
	562.,545.,529.,512.,495.,478.,461.,446.,430.,416.,403.,391.,379.,366.,$
	349.,332.,312.,291.,274.,263.,252.,244.,236.,237.,233.,228.,223.,217.,$
	211.,203.,194.,184.,176.,171.,167.,164.,162.,162.,163.,170.]

; Calculate the atmospheric temperature and pressure gradients:


if (zkm gt 100) then begin
	if (z lt 100010) then print, ' Plume extends beyond atmos data.'
    P=Parray(51)
    Ta=Tarray(51)
endif else begin
    P1=Parray(i)
    T1=Tarray(i)
    if (deltazi eq 0) then begin
      P = P1*100000
      Ta = T1
    endif else begin
      P2 = Parray(i+1)
      T2 = Tarray(i+1)
      z1 = i*2
      z2 = (i+1)*2
      mP = (z2-z1)/(alog(P2)-alog(P1))
      lnP = (zkm-z1)/mP + alog(P1)
      P = (exp(lnP))*100000
      mT = (z2-z1)/(T2-T1)
      Ta = (zkm-z1)/mT + T1
    endelse
endelse
      
   return, [Ta,P]

end