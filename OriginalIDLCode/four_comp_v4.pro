function calc, z, y

common block,Rd,Rv,Cd,Cv,Cl,Cs,rhol,rhos,MsF,delz,g,alpha,atmos,humidity,tau,omega,latht

; Calculate the atmospheric temperature, pressure, and mixing ratio (if appropriate) at z:

if atmos eq 'mars' then TP = atmos_mars(z,g)
if atmos eq 'paleomars' then TP = atmos_paleomars(z,g)
if atmos eq 'venus-average' then TP = atmos_venus_vira(z,g)
if atmos eq 'venus-low' then TP = atmos_venus_low(z,g)
if atmos eq 'venus-high' then TP = atmos_venus_high(z,g)
if atmos eq 'venus-vega' then TP = atmos_venus_vega2(z,g)
if atmos eq 'earth' then TP = atmos_earth(z,g)
if atmos eq 'hekla' then TP = atmos_Jan_60(z,g)
if atmos eq 'wet' then TP = atmos_wet(z,g)
if atmos eq 'dry' then TP = atmos_dry(z,g)
if atmos eq 'decaan' then TP = atmos_decaan(z,g)
if atmos eq 'io' then TP = atmos_io(z,g)


Ta = TP(0)
P = TP(1)
wa = mixing_ratio(z,humidity)

; Calculate the ambient air density, rhoaB, at z:

rhoadphiad = (P/(Rv*Ta))/(wa + Rd/Rv)
rhoavphiav = rhoadphiad*wa
rhoaB = rhoavphiav + rhoadphiad
CaB = (Cd + wa*Cv)/(1 + wa)

MdF = y(0)
nuF = y(1)
TF  = y(2)
MlF = y(3)
MvF = y(4)

m = MsF + MvF + MdF + MlF
CB = (MsF*Cs + MvF*Cv + MdF*Cd + MlF*Cl)/m
varm = (P*m*CB)/(Rv*TF)
rhoB = m*varm/(MvF+MdF*Rd/Rv)

info = {MdF:MdF,MvF:MvF,MlF:MlF,nuF:nuF,TF:TF,m:m,rhoaB:rhoaB,varm:varm,rhoB:rhoB,CB:CB,$
	Ta:Ta,P:P,rhoadphiad:rhoadphiad,rhoavphiav:rhoavphiav,CaB:CaB}
return, info
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function outpro, z, y, NBH

; This procedure prints out values for the plume parameters

info = calc(z,y)

r = sqrt(info.m^2/(info.rhoB*info.nuF))
u = info.nuF/info.m
theta = info.TF/(info.m*info.CB)
rhoB = info.rhoB
rhoaB = info.rhoaB

printf,1, z, r, u, theta, rhoB, rhoaB
printf,2,z,rhoaB

if abs(rhoaB - rhoB) lt 0.0001 then NBH = z

return, NBH
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function revcalc, z, y, yout

common block,Rd,Rv,Cd,Cv,Cl,Cs,rhol,rhos,MsF,delz,g,alpha,atmos,humidity,tau,omega,latht

info_out = calc(z+delz, yout)
info = calc(z,y)

uout = info_out.nuF/info_out.m
u = info.nuF/info.m

CBout = (info_out.MdF*Cd + info_out.MvF*Cv + MsF*Cs + info_out.MlF*Cl)/info_out.m
CB = (info.MdF*Cd + info.MvF*Cv + MsF*Cs + info.MlF*Cl)/info.m

thout = info_out.TF/(info_out.m*CBout)
th = info.TF/(info.m*CB)

rout = sqrt((info_out.m^2)/(info_out.rhoB*info_out.nuF))
r = sqrt((info.m^2)/(info.rhoB*info.nuF))

info = {uout:uout,u:u,rout:rout,r:r,thout:thout,th:th}
return, info
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function toptest, z, y, yout, u_cutoff

; This function calculates the values of u, uout, r, and
; rout at each time step in order to determine if the top of
; the plume has been reached.

top = 0

result = revcalc (z, y, yout)
if (result.uout le u_cutoff) or $
	 (result.rout lt 0) or $
	 (abs((result.rout-result.r)/result.r) gt 0.4) then top = 1

return, top
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function plumederivs, z, y

; This function calculates the values of the derivatives of MdF, MvF,
; MlF, nuF, and TF, for each step in z.

common block,Rd,Rv,Cd,Cv,Cl,Cs,rhol,rhos,MsF,delz,g,alpha,atmos,humidity,tau,omega,latht

dydz = dblarr(5)

info = calc(z,y)

c1 = (info.MvF+info.MdF*(Rd/Rv))/info.varm

dMdFdz = 2*alpha*info.rhoadphiad*sqrt(info.nuF*c1/info.m)
dnuFdz = g*info.m*(info.rhoaB*c1-info.m)/info.nuF
dMlFdz = omega*info.MvF*info.m/info.nuF
dMvFdz = dMdFdz*info.rhoavphiav/info.rhoadphiad-dMlFdz
dTFdz = dMdFdz*info.rhoaB*info.CaB*info.Ta/info.rhoadphiad + latht*dMlFdz - info.rhoaB*g*(c1-(info.MlF/rhol + MsF/rhos))

dydz(0) = dMdFdz
dydz(1) = dnuFdz
dydz(2) = dTFdz
dydz(3) = dMlFdz
dydz(4) = dMvFdz

return, dydz
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function cond_rate, z,y

; This function determines the value of omega (condensation
; rate) based on whether or not there is sufficient
; water vapor available and if the temperature has dropped
; below the critical temperature.

common block,Rd,Rv,Cd,Cv,Cl,Cs,rhol,rhos,MsF,delz,g,alpha,atmos,humidity,tau,omega,latht

info=calc(z,y)
theta = info.TF/(info.m*info.CB)
MvF = info.MvF

; we use the Clausius-Clapeyron relationship to determine the critical
; temperature for condensation. In this equation, the reference
; pressure is 101300 Pa (standard pressure) and the boiling point
; of water at standard pressure, 373 K.

theta_crit = (latht*373.)/(latht-Rv*373.*alog(info.P/101300.))
if (info.MvF gt 0) and (theta lt theta_crit) then begin
  if tau ne 0 then omega = 1./tau
  if tau eq 0 then omega = 0 
endif else begin
  omega = 0.
endelse
   
return, omega
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function define, ifile

; This function opens up and reads the file containing
; all the initial conditions. It then defines all constants
; for the appropriate planet, and calculates the initial
; values for all variables.

common block,Rd,Rv,Cd,Cv,Cl,Cs,rhol,rhos,MsF,delz,g,alpha,atmos,humidity,tau,omega,latht

input = strarr(8)
close,1
openr,1,ifile
readf,1,input
close,1

tmp = strarr(2,8)
for i = 0,7 do tmp(*,i) = strtrim(str_sep(input(i),';'),2)

r0 = float(tmp(1,0))
u0 = float(tmp(1,1))
theta0 = float(tmp(1,2))
z0 = float(tmp(1,3))
n0 = float(tmp(1,4))
atmos = tmp(1,5)
tau   = tmp(1,6)
ofile = tmp(1,7)

print,'r0     = ',r0
print,'u0     = ',u0
print,'theta0 = ',theta0
print,'z0     = ',z0
print,'n0     = ',n0
print,'planet = ',atmos
print,'tau    = ',tau

openw,2,'rhoaB.txt'
openw,1,ofile
printf,1,'r0     = ',r0
printf,1,'u0     = ',u0
printf,1,'theta0 = ',theta0
printf,1,'z0     = ',z0
printf,1,'n0     = ',n0
printf,1,'planet = ',atmos
printf,1,'tau    = ',tau

printf,1,'       z(m)          r(m)            u(m/s)         theta(K)       rhoB(kg/m^3)    rhoaB(kg/m^3)'

rhol = 1000.
rhos = 1000.
alpha = 0.09

Cv = 2000.  ; water vapor
Cl = 4200.  ; liquid water
Cs = 920.	; pumice
Rv = 461.	; water vapor

latht = 2.257e6 ; latent heat of condensation for water

if r0 ge 50 then delz = 5.0 else delz = r0/10

Case atmos of

	'earth':	Begin
				g = 9.8
				Cd = 998.
				Rd = 287.
				TP = atmos_earth(z0,g)
				humidity = 'null'
				End

	'hekla':	Begin
				g = 9.8
				Cd = 998.
				Rd = 287.
				TP = atmos_Jan_60(z0,g)
				humidity = 'null'
        End

  'wet':  Begin
        g = 9.8
        Cd = 998.
        Rd = 287.
        TP = atmos_wet(z0,g)
        humidity = 'wet'
        End

  'dry':  Begin
        g = 9.8
        Cd = 998.
        Rd = 287.
        TP = atmos_dry(z0,g)
        humidity = 'dry'
        End

  'decaan':  Begin
        g = 9.8
        Cd = 998.  ; same as for US Standard
        Rd = 286.6 ; increased CO2 by four times 
        TP = atmos_decaan(z0,g)
        humidity = 'null'
        End

	'venus-average':	Begin
				;Cv = 1652. ; 70% H2O; 15% CO2; 8% N2; 7% SO2
        ;Rv = 385.  ; 70% H2O; 15% CO2; 8% N2; 7% SO2
        g = 8.4
				Cd = 835.
				Rd = 191.
				humidity = 'null'
        TP = atmos_venus_vira(z0,g)
				End

  'venus-low':  Begin
        ;Cv = 1652. ; 70% H2O; 15% CO2; 8% N2; 7% SO2
        ;Rv = 385.  ; 70% H2O; 15% CO2; 8% N2; 7% SO2
        g = 8.4
        Cd = 835.
        Rd = 191.
        TP = atmos_venus_low(z0,g)
        humidity = 'null'
        End

  'venus-high':  Begin
;        Cv = 1652. ; 70% H2O; 15% CO2; 8% N2; 7% SO2
;        Rv = 385.  ; 70% H2O; 15% CO2; 8% N2; 7% SO2
;        Cv = 1237. ; 42% H2O; 15% CO2; 0% N2; 43% SO2
;        Rv = 278.  ; 42% H2O; 15% CO2; 0% N2; 43% SO2
        g = 8.4
        Cd = 835.
        Rd = 191.
        TP = atmos_venus_high(z0,g)
        humidity = 'null'
        End

  'venus-vega':  Begin
        ;Cv = 1652. ; 70% H2O; 15% CO2; 8% N2; 7% SO2
        ;Rv = 385.  ; 70% H2O; 15% CO2; 8% N2; 7% SO2
        g = 8.4
        Cd = 835.
        Rd = 191.
        TP = atmos_venus_vega2(z0,g)
        humidity = 'null'
        End

	'mars':		Begin
				;Cv = 1400.	; same as atmos
				;Rv = 350.	; same as atmos
				g = 3.7
				Cd = 835.
				Rd = 191.
				TP = atmos_mars(z0,g)
				humidity = 'null'
        End

	'paleomars':Begin
				g = 3.7
				Cd = 835.
				Rd = 191.
				TP = atmos_paleomars(z0,g)
				humidity = 'null'
        End

	'io':		Begin
				Cv = 1317.	; 50/50 SO2 & water vapor
				Rv = 296.	; 50/50 SO2 & water vapor
				g = 1.8
				Cd = 634.
				Rd = 130.
				TP = atmos_io(z0,g)
				humidity = 'null'
        End

EndCase

P0 = double(TP(1))
Ta0 = TP(0)
z = z0

rhov0 = P0/(Rv*double(theta0))
rhod0 = P0/(Rd*double(theta0))

w_a = mixing_ratio(z,humidity)
rhoaB0 = (P0/(Rv*double(Ta0)))*(1+w_a)/(w_a+Rd/Rv)

phis0 = double(((1-n0)*rhov0)/(n0*rhos+(1-n0)*rhov0))
phid0 = double(0)
phil0 = double(0)
phiv0 = double(1-phis0-phid0-phil0)

rhoB0 = rhov0*phiv0 + rhod0*phid0 + rhol*phil0 + rhos*phis0

MsF = double(rhos*phis0*u0*r0^2)
MlF = double(rhol*phil0*u0*r0^2)
MvF = double(rhov0*phiv0*u0*r0^2)
MdF = double(rhod0*phid0*u0*r0^2)

m0 = double(MsF + MvF + MdF + MlF)
CB0 = (MsF*Cs + MvF*Cv + MdF*Cd + MlF*Cl)/m0

nuF = double(m0*u0)
TF = double(m0*CB0*theta0)

y = dblarr(5)

y(0) = MdF
y(1) = nuF
y(2) = TF
y(3) = MlF
y(4) = MvF

omega = cond_rate(z,y)

defs = {z:z, y:y}
return, defs
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro four_comp_v4, ifile

; This is the core of the code for solving the system of
; five differential equations describing the thermal
; dynamics of buoyant eruption columns.  In this case,
; the plumes are comprised of dry air, magmatic gas, liquid,
; and solids.
; 
; The ambient atmosphere may have a different composition
; to the magmatic gas.  The program takes as input the
; initial values for gas content, plume radius, bulk
; column velocity, and column temperature, and solves five
; differnetial equations, for plume density, radius, velocity
; and temperature as functions of z, using the Runge-Kutta
; method. This IDL code is based on the plumes.pro idl code 
; (known to be working at the beginning of this development
; on 6/17/2011) developed for 2 component plumes (ash and gas).
; I have used my old fortran code (PhD Thesis Appendix) as
; a guide for the subroutines required for a 4 component plume.
; 
; Progrm four_comp confirmed to be working:
; with no changes made from plumes.pro other than file name: 6/17/2011
; with changes to dry and vapor variable names: 6/17/2011
; with changes to density, vol fraction and flux variable names: 6/17/2011
; with change to 5 equations, omega = 0: 6/17/2011
; with calculation and printout of NBH (v2.1): 6/17/2011
; with input of tau instead of omega (v2.2): 6/17/2011
; with tau and omega in common block, start cond_rate function (v2.3): 6/17/2011
; with new cond_rate function appears to initiate condensation
; correctly and predict reasonable results (v2.4): 6/17/2011
; removed "Ra" from call to all atmos procedures (v3): 6/21/2011
; added mixing ratio calc to define procedure (v3.1): 6/21/2011
; changed definition of rhoaB in calc procedure (v3.2): 6/21/2011
; changed definitions of dervis in plumederivs procedure (v3.3): 6/21/2011
; v3.3 works for the 'null' case when the normal earth atmospheres
; are used. have not tried this with a humid atmosphere.
; v4.0 was tested for both wet and dry atmospheres and for no, to moderate,
; to rapid condensation. Results from Glaze et al. [1997] were reproduced
; to within 100 m predicted max H.
;
; The variables and constants used in the program are:
; g		  gravity [m s^(-2)] = 9.8 e; 3.7 m; 8.4 v
; rhoaB0	ambient density 0 km above mean radius, Ta = T0, P = P0)
; alpha	entrainment constant = 0.09 (Woods, 1988; Morton
;		    et al., 1956)
; Cd	  specific heat for ambient [J K^(-1) kg^(-1)] = 998 e, 835 m,
;		    835 v (m & v are a guess based on a mainly CO2 atmosphere)
; Cv	  specific heat of the magmatic gas = Cd, for
;		    similar composition, = 2000 J K^(-1) kg^(-1) for
;		    water vapor
;	Cl    specific heat of liquid water = 4200 J K^(-1) kg^(-1)
;	      (Handbook of Chemistry and Physics)
; Cs	  specific heat of ash particles = 920 J K^(-1) kg^(-1)
;		    (Riehle, 1973)
; rhoB	bulk density of the plume, defined by the weighted
;       average of the plume constituents
; rhol  density of liquid water = 1000 kg m^(-3)
; rhos	density of pumice = 1000 kg m^(-3) (Wilson, 1976)
; Rd	  gas constant for ambient [J K^(-1) kg^(-1)], 287 e, 191 m, 191 v
;		    (m & v are a guess based on a mainly CO2 atmosphere)
; Rv	  gas constant for magmatic gas = Rd for similar
;		    composition, = 461 J K^(-1) kg^(-1) for water vapor
; T0	  temperature of atmosphere 0 km above mean radius [K]
; P0	  pressure of atmospehre 0 km above mean radius [Pa]
; delz	step size increment = 5 m
; phi	  volume fraction of solids in the control volume
; latht latent heat of condensation for water vapor = 2.257e6 J kg^(-1)
; omega condensation rate s^(-1)

common block,Rd,Rv,Cd,Cv,Cl,Cs,rhol,rhos,MsF,delz,g,alpha,atmos,humidity,tau,omega,latht

defs = define(ifile)

top = 0
if delz le 1.0 then u_cutoff = 0.1 else u_cutoff = 10
z = defs.z
y = defs.y
NBH = 0.0
NBH = outpro(z,y,NBH)

repeat begin

	omega = cond_rate(z,y)
	dydz = plumederivs(z, y)
	yout = rk4(y, dydz, z, delz, 'plumederivs',/double)
	top = toptest(z, y, yout, u_cutoff)
	z = z + delz
	y = yout
	NBH = outpro(z,y,NBH)

endrep until (top eq 1)

print,'Neutral Buoyancy Height (km):  ',NBH/1000
print,'Plume Height (km):  ',z/1000

close,1
close,2

end