function atmos_jan_60, z, g

; This subroutine calculates the ambient atmospheric
; temperature and pressure for Earth at 60 degrees
; latitude in January. Values are taken from LSG
; PhD thesis fortran code appendix.

; Define atmospheric parameters

T0 = 255.0
P0 = 1*101300
zkm = z/1000
H = [1., 4., 9., 15., 25., 34., 50.]
mu = [-3.0, 2., 7.2, 0.0, 0.6, -1.1111, -2.5]

; Calculate the atmospheric temperature and pressure gradients:

CASE 1 of

	(zkm lt H(0)):						Ta = T0-mu(0)*zkm
	(zkm ge H(0)) and (zkm le H(1)):	Ta = T0-mu(0)*H(0)-mu(1)*(zkm-H(0))
	(zkm gt H(1)) and (zkm le H(2)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(zkm-H(1))
	(zkm gt H(2)) and (zkm le H(3)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(H(2)-H(1))-mu(3)*(zkm-H(2))
	(zkm gt H(3)) and (zkm le H(4)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(H(2)-H(1))-mu(3)*(H(3)-H(2))-mu(4)*(zkm-H(3))
	(zkm gt H(4)) and (zkm le H(5)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(H(2)-H(1))-mu(3)*(H(3)-H(2))-mu(4)*(H(4)-H(3))-mu(5)*(zkm-H(4))
	(zkm gt H(5)) and (zkm le H(6)):    Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(H(2)-H(1))-mu(3)*(H(3)-H(2))-mu(4)*(H(4)-H(3))-mu(5)*(H(5)-H(4))-mu(6)*(zkm-H(5))
	(zkm gt H(6)):	stop

ENDCASE

P = P0*exp(-zkm/7.)
return, [Ta,P]
end