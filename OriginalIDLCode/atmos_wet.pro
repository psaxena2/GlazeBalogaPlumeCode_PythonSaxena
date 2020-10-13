function atmos_wet, z, g

; This subroutine calculates the latitudinally averaged
; ambient atmospheric temperature and pressure for Earth
; at 10 degrees latitude in July 1984. Values are
; estimated from Figure 2 [Glaze et al., 1997] and derived
; in the file "Mixing Ratio Data.xls". The Glaze et al. figure
; states that Temperature and associated mixing ratio
; data are unpublished data from Dave Crisp. 

; Define atmospheric parameters

T0 = 293.4006
P0 = 1*101300
zkm = z/1000
H = [2.1, 15.7, 20.]
mu = [3.06, 6.01, 1.98]

; Calculate the atmospheric temperature and pressure gradients:

CASE 1 of

	(zkm lt H(0))                  :  Ta = T0-mu(0)*zkm
	(zkm ge H(0)) and (zkm le H(1)):	Ta = T0-mu(0)*H(0)-mu(1)*(zkm-H(0))
	(zkm gt H(1)) and (zkm le H(2)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(zkm-H(1))
	(zkm gt H(2)):	stop

ENDCASE

P = P0*exp(-zkm/7.)
return, [Ta,P]
end