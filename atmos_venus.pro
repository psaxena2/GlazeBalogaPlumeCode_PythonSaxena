function atmos_venus, z, g

; This subroutine calculates the ambient atmospheric
; temperature and pressure for the average Venus Atmosphere
; Pressure and Temperatures are given for every 1 km of altitude
; Pressures are given in Pa/10^5
; Glaze [1999] states that vaues for T/P are taken from Seiff 
; (1983) "Thermal structure of the atmosphere of Venus", 
; Chapter 11, in Venus edited by D. Hunten et al., 215 - 279, 
; U of Ariz Press. HOWEVER, there are no data tables given in Seiff (1983).
; Data tables in Seiff et al. (1985) - VIRA paper - are not the same
; as those in here.

; Define atmospheric parameters

T0 = 743.0
P0 = 98.1193*100000
zkm = z/1000
zi = fix(zkm)
deltazi = zkm-float(zi)
i = zi

Parray = [98.1193,92.100,86.4500,81.0900,76.0100,71.2000,66.6500,62.3500,$
	58.2800,54.4400,50.8100,47.3900,44.1600,41.1200,38.2600,35.5700,$
	33.0400,30.6600,28.4300,26.3300,24.3600,22.5200,20.7900,19.1700,$
	17.6600,16.2500,14.9300,13.7000,12.5800,11.4900,10.5000,9.5810,$
	8.7290,7.9400,7.2110,6.5370,5.9170,5.3460,4.8220,4.3420,3.9030,$
	3.5010,3.1350,2.8020,2.4990,2.2260,1.9790,1.7560]

Tarray = [743.0,735.3,727.7,720.2,712.4,704.6,696.8,688.8,681.1,673.6,$
	665.8,658.2,650.6,643.2,635.5,628.1,620.6,613.3,605.2,597.1,589.3,$
	580.7,572.4,564.3,556.0,547.5,539.2,530.7,522.3,513.8,505.6,496.9,$
	488.3,479.9,471.7,463.4,455.5,448.0,439.9,432.5,425.1,417.6,410.0,$
	403.5,397.1,391.2,385.4,379.7]

dlnParray = [.06331,.06331,.06401,.06469,.06537,.06604,.06669,.06750,$
	.06816,.06901,.06968,.07059,.07132,.07209,.07290,.07378,.07476,$
	.07551,.07674,.07777,.07854,.07993,.08113,.08204,.08321,.08472,$
	.08598,.08529,.09063,.09010,.09159,.09313,.09474,.09631,.09813,$
	.09965,.10148,.10316,.10485,.10659,.10870,.11042,.11230,.11444,$
	.11568,.11761,.11955,.12092]

dTarray = [7.7,7.6,7.5,7.8,7.8,7.8,8.0,7.7,7.5,7.8,7.6,7.6,7.4,7.7,7.4,$
	7.5,7.3,8.1,8.1,7.8,8.6,8.3,8.1,8.3,8.5,8.3,8.5,8.4,8.5,8.2,8.7,8.6,$
	8.4,8.2,8.3,7.9,7.5,8.1,7.4,7.4,7.5,7.6,6.5,6.4,5.9,5.8,5.7,6.6]


; Calculate the atmospheric temperature and pressure gradients:


if (zkm gt 47) then begin
	if (z lt 47005) then print, ' Plume extends beyond atmos data.'
    Pi=Parray(47)
    Ti=Tarray(47)
    dlnPi=-dlnParray(47)
    dTi=-dTarray(47)
endif else begin
    Pi=Parray(i)
    Ti=Tarray(i)
    dlnPi=-dlnParray(i)
    dTi=-dTarray(i)
endelse

      Ta=Ti+dTi*deltazi
      P=(exp(alog(Pi)+dlnPi*deltazi))*(100000)

   return, [Ta,P]

end