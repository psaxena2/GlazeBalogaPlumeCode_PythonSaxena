function atmos_venus_high, z, g

; This subroutine calculates the ambient atmospheric
; temperature and pressure for the Venus Atmosphere at low latitudes
; Pressure and Temperatures are given for every 2 km of altitude

; Define atmospheric parameters

T0 = 737.0
P0 = 92.1*100000
zkm = z/1000
zi = fix(zkm)
deltazi = zkm-2.*float(zi/2)
i = zi/2

Parray = [92.100,81.200,71.300,62.500,54.600,47.500,41.200,35.600,$
	30.700,26.300,22.500,19.10,16.100,13.5600,11.3400,9.430,$
	7.800,6.40,5.220,4.230,3.40,2.710,2.150,1.69,1.316,1.014,0.772,$
	0.577,0.422,0.302,0.210,0.143,0.095,0.0626,0.0422,0.0289,0.0197,0.0134,$
	0.00906,0.00604,0.00396,0.00256,0.00162,0.00101,0.000613,0.000365,$
	0.000213,0.000123,0.0000715,0.0000421,0.0000253]

Tarray = [737.,722.,706.,691.,675.,659.,644.,627.,611.,594.,577.,$
	561.,544.,527.,509.,493.,477.,461.,445.,429.,415.,400.,388.,375.,359.,$
	344.,324.,303.,282.,262.,243.,228.,217.,219.,240.,237.,236.,234.,227.,$
	218.,210.,202.,194.,185.,177.,170.,165.,163.,165.,173.,177.]

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