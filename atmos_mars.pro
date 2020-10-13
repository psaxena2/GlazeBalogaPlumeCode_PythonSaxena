function atmos_mars, z, g

; This subroutine calculates the ambient atmospheric
; temperature and pressure for Mars that match those
; given in Hort and Weitz.

; Define atmospheric parameters

T0 = 225
P0 = .01*101300
zkm = z/1000
H = [5.,30.,44.,133.,140.,150]
mu = [4.,1.,2.36,0.169,-5.71,-2.4]

; the following H and mu parameters are for the estimate
; of Pathfinder temperatures:

;H = [5., 50., 62., 133., 140., 150.]
;mu = [4., 1., 2.92, 0., -6.71, -2.4]

; Calculate the atmospheric temperature and pressure gradients:

P = P0*exp(-zkm/8.0)

CASE 1 of

	(zkm lt H(0)):						Ta = T0-mu(0)*zkm
	(zkm ge H(0)) and (zkm le H(1)):	Ta = T0-mu(0)*H(0)-mu(1)*(zkm-H(0))
	(zkm gt H(1)) and (zkm le H(2)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(zkm-H(1))
	(zkm gt H(2)) and (zkm le H(3)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(H(2)-H(1))-mu(3)*(zkm-H(2))
	(zkm gt H(3)) and (zkm le H(4)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(H(2)-H(1))-mu(3)*(H(3)-H(2))-mu(4)*(zkm-H(3))
	(zkm gt H(4)) and (zkm le H(5)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(H(2)-H(1))-mu(3)*(H(3)-H(2))-mu(4)*(H(4)-H(3))-mu(5)*(zkm-H(4))
	(zkm gt H(5)):						stop

ENDCASE

return, [Ta,P]
end