function atmos_decaan, z, g

; This subroutine calculates the ambient atmospheric
; temperature and pressure for Earth during the Maastrichtian, 
; relevant to he Decaan Traps flood basalt eruptions

; Define atmospheric parameters

T0 = 305. ; Kelvin
P0 = 1*101300 ; Pascal
zkm = z/1000
H = [9.,15.5,21., 23.5, 27.5]
mu = [6.22,6.12, 3.63, 1.47, -1.64]

;exp1 = (g*1000)/(Ra*mu(0))
;exp2 = (g*1000)/(Ra*mu(1))

; Calculate the atmospheric temperature and pressure gradients:

CASE 1 of

	(zkm lt H(0)):						Ta = T0-mu(0)*zkm
	(zkm ge H(0)) and (zkm le H(1)):	Ta = T0-mu(0)*H(0)-mu(1)*(zkm-H(0))
	(zkm gt H(1)) and (zkm le H(2)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(zkm-H(1))
	(zkm gt H(2)) and (zkm le H(3)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(H(2)-H(1))-mu(3)*(zkm-H(2))
	(zkm gt H(3)) and (zkm le H(4)):	Ta = T0-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(H(2)-H(1))-mu(3)*(H(3)-H(2))-mu(4)*(zkm-H(3))
	(zkm gt H(4)):	stop

ENDCASE

P = P0*exp(-zkm/7.6)
return, [Ta,P]
end