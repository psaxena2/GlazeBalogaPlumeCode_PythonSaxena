function mixing_ratio, z, humidity

; This subroutine calculates the latitudinally averaged
; water vapor mixing ratio (in kg/kg) on Earth. Two choices
; are for a "dry" atmosphere at 60 degrees latitude in January
; 1984 or a "wet" atmosphere at 10 degrees lattiude in July.
; Values are estimated from Figure 2 [Glaze et al., 1997] and derived
; in the file "Mixing Ratio Data.xls". The Glaze et al. figure
; states that Temperature and associated mixing ratio
; data are unpublished data from Dave Crisp. 

; Define atmospheric parameters

zkm = z/1000

; Calculate the atmospheric mixing ratio gradient:

CASE humidity of
  'null': w_a = 0.0
  
  'wet':  BEGIN
        w_o = 0.0129
        H = [4.6, 7.2, 9.4, 12.,20]
        mu = [0.001964, 0.000882, 0.000383, 0.000189, 0.000041]
        CASE 1 of
        (zkm lt H(0))                  :  w_a = w_o-mu(0)*zkm
        (zkm ge H(0)) and (zkm le H(1)):  w_a = w_o-mu(0)*H(0)-mu(1)*(zkm-H(0))
        (zkm gt H(1)) and (zkm le H(2)):  w_a = w_o-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(zkm-H(1))
        (zkm gt H(2)) and (zkm le H(3)):  w_a = w_o-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(H(2)-H(1))-mu(3)*(zkm-H(2))
        (zkm gt H(3)) and (zkm le H(4)):  w_a = w_o-mu(0)*H(0)-mu(1)*(H(1)-H(0))-mu(2)*(H(2)-H(1))-mu(3)*(H(3)-H(2))-mu(4)*(zkm-H(3))
        (zkm gt H(4)):  stop
        ENDCASE
    END
  
  'dry':  BEGIN
        w_o = 0.0016
        H = [4.6,20]
        mu = [0.000250, 0.000064]
        CASE 1 of
        (zkm lt H(0))                  :  w_a = w_o-mu(0)*zkm
        (zkm ge H(0)) and (zkm le H(1)):  w_a = w_o-mu(0)*H(0)-mu(1)*(zkm-H(0))
        (zkm gt H(1)):  stop
        ENDCASE
   END
ENDCASE

return, w_a
end