function atmos_dry, z, g

; This subroutine calculates the latitudinally averaged
; ambient atmospheric temperature and pressure for Earth
; at 60 degrees latitude in January 1984. Values are
; estimated from Figure 2 [Glaze et al., 1997] and derived
; in the file "Mixing Ratio Data.xls". The Glaze et al. figure
; states that Temperature and associated mixing ratio
; data are unpublished data from Dave Crisp. 

; Define atmospheric parameters

T0 = 258.7912
P0 = 1*101300
zkm = z/1000
H = [2.3, 9.2, 11.4, 13.4, 16.8, 20.]
mu = [2.14, 4.76, 2.55, 0.55, 0.,0.68]

; Calculate the atmospheric temperature and pressure gradients:

CASE 1 of

	(zkm lt H(0))                  :  Ta = T0-mu(0)*zkm
	(zkm ge H(0)) and (zkm le H(1)):	Ta = T0-mu(0)*H(0)-mu(1)*(zkm-H(0))
	(zkm gt H(1)) and (zkm le H(2)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(zkm-H(1))
	(zkm gt H(2)) and (zkm le H(3)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(H(2)-H(1))-mu(3)*(zkm-H(2))
	(zkm gt H(3)) and (zkm le H(4)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(H(2)-H(1))-mu(3)*(H(3)-H(2))-mu(4)*(zkm-H(3))
	(zkm gt H(4)) and (zkm le H(5)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(H(2)-H(1))-mu(3)*(H(3)-H(2))-mu(4)*(H(4)-H(3))-mu(5)*(zkm-H(4))
	(zkm gt H(5)):	stop

ENDCASE

P = P0*exp(-zkm/7.)
return, [Ta,P]
end