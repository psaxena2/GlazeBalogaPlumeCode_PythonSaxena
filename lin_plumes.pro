function calc, z, y

common block,Ra,Rm,Ca,Cm,Cs,rhos,MsF,MmF,delz,g,alpha,atmos

; Calculate the atmospheric temperature and pressure at z:

if atmos eq 'mars' then TP = atmos_mars(z,g)
if atmos eq 'paleomars' then TP = atmos_paleomars(z,g)
if atmos eq 'venus-average' then TP = atmos_venus_vira(z,g)
if atmos eq 'venus-low' then TP = atmos_venus_low(z,g)
if atmos eq 'venus-high' then TP = atmos_venus_high(z,g)
if atmos eq 'venus-vega' then TP = atmos_venus_vega2(z,g)
if atmos eq 'earth' then TP = atmos_earth(z,g)
if atmos eq 'hekla' then TP = atmos_Jan_60(z,g)
if atmos eq 'decaan' then TP = atmos_decaan(z,g)
if atmos eq 'miocene' then TP = atmos_miocene(z,g)
if atmos eq 'io' then TP = atmos_io(z,g)

Ta = TP(0)
P = TP(1)

; Calculate the ambient air density, rhoa, at z:

rhoa = P/(Ra*Ta)

MeF = y(0)
nuF = y(1)
TF  = y(2)

m = MsF + MmF + MeF
CB = (MsF*Cs + MmF*Cm + MeF*Ca)/m
varm = (P*m*CB)/(Rm*TF)
ub = (MmF + MeF*(Ra/Rm))/varm
rhoB = m/ub

info = {MeF:MeF,nuF:nuF,TF:TF,m:m,rhoa:rhoa,varm:varm,rhoB:rhoB,CB:CB,$
	Ta:Ta,P:P}
return, info
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro outpro, z, y

; This procedure prints out values for the plume parameters

info = calc(z,y)

b = info.m^2/(info.rhoB*info.nuF)
u = info.nuF/info.m
theta = info.TF/(info.m*info.CB)
rhoB = info.rhoB
rhoa = info.rhoa

printf,1, z, b, u, theta, rhoB; rhoa was removed 5/9/05 to compare to new model
printf,2,z,rhoa
;printf,1, z, b, u, theta, rhoB, rhoa ; this is the working line of code
;printf,1, z, b, u, y(0), rhoa

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function revcalc, z, y, yout

common block,Ra,Rm,Ca,Cm,Cs,rhos,MsF,MmF,delz,g,alpha,atmos

info_out = calc(z+delz, yout)
info = calc(z,y)

uout = info_out.nuF/info_out.m
u = info.nuF/info.m

CBout = (info_out.MeF*Ca + MmF*Cm + MsF*Cs)/info_out.m
CB = (info.MeF*Ca + MmF*Cm + MsF*Cs)/info.m

thout = info_out.TF/(info_out.m*CBout)
th = info.TF/(info.m*CB)

bout = (info_out.m^2)/(info_out.rhob*info_out.nuF)
b = (info.m^2)/(info.rhob*info.nuF)

info = {uout:uout,u:u,bout:bout,b:b,thout:thout,th:th}
return, info
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function toptest, z, y, yout, u_cutoff

; This function calculates the values of u, uout, b, and
; bout at each time step in order to determine if the top of
; the plume has been reached.

top = 0

result = revcalc (z, y, yout)
if (result.uout le u_cutoff) or $
	 (result.bout lt 0) or $
	 (abs((result.bout-result.b)/result.b) gt 0.4) then top = 1

return, top
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function plumederivs, z, y

; This function calculates the values of the derivatives of b, u,
; and theta, for each step in z.

common block,Ra,Rm,Ca,Cm,Cs,rhos,MsF,MmF,delz,g,alpha,atmos

dydz = dblarr(3)

info = calc(z,y)

c1 = (MmF+info.MeF*(Ra/Rm))/info.varm
c2 = (MsF/rhos)-c1

dMeFdz = info.rhoa*alpha*(info.nuF/info.m)
dnuFdz = g*info.m*(info.rhoa*c1-info.m)/info.nuF
dTFdz = ((info.rhoa*alpha*info.nuF*Ca*info.Ta)/info.m)+info.rhoa*g*c2

dydz(0) = dMeFdz
dydz(1) = dnuFdz
dydz(2) = dTFdz

return, dydz
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function define, ifile

; This function opens up and reads the file containing
; all the initial conditions. It then defines all constants
; for the appropriate planet, and calculates the initial
; values for all variables.

common block,Ra,Rm,Ca,Cm,Cs,rhos,MsF,MmF,delz,g,alpha,atmos

input = strarr(7)
close,1
openr,1,ifile
readf,1,input
close,1

tmp = strarr(2,7)
for i = 0,6 do tmp(*,i) = strtrim(str_sep(input(i),';'),2)

b0 = float(tmp(1,0))
u0 = float(tmp(1,1))
theta0 = float(tmp(1,2))
z0 = float(tmp(1,3))
n0 = float(tmp(1,4))
atmos = tmp(1,5)
ofile = tmp(1,6)

print,'b0     = ',b0
print,'u0     = ',u0
print,'theta0 = ',theta0
print,'z0     = ',z0
print,'n0     = ',n0
print,'planet = ',atmos

openw,2,'rhoa.txt'
openw,1,ofile
printf,1,'b0     = ',b0
printf,1,'u0     = ',u0
printf,1,'theta0 = ',theta0
printf,1,'z0     = ',z0
printf,1,'n0     = ',n0
printf,1,'planet = ',atmos

printf,1,'       z(m)          b(m)            u(m/s)         theta(K)       rhoB(kg/m^3)    rhoa(kg/m^3)'
;printf,1,'       z(m)          b(m)            u(m/s)         Ent Rate(kg/s)    rhoa(kg/m^3)'

rhos = 1000.
alpha = 0.09

Cm = 2000.	; water vapor
Cs = 920.	; pumice
Rm = 461.	; water vapor

if b0 ge 50 then delz = 5.0 else delz = b0/10

Case atmos of

	'earth':	Begin
        Cm = 1241. ; 42% H2O; 15% CO2; 43% SO2 (Gerlach, 1980)
        Rm = 279.  ; 42% H2O; 15% CO2; 43% SO2 (Gerlach, 1980)
				g = 9.8
				Ca = 998.
				Ra = 287.
				TP = atmos_earth(z0,g)
				End

	'hekla':	Begin
				g = 9.8
				Ca = 998.
				Ra = 287.
				TP = atmos_Jan_60(z0,g)
				End

  'decaan':  Begin
        g = 9.8
        Cd = 998.  ; same as for US Standard
        Rd = 286.6 ; increased CO2 by four times 
        TP = atmos_decaan(z0,g)
        humidity = 'null'
        End

   'miocene':  Begin
          g = 9.8
          Ca = 998.  ; same as for US Standard
          Ra = 287.  ; same as for US Standard
          TP = atmos_miocene(z0,g)
          humidity = 'null'
        End

  'venus-average':  Begin
        ;Cm = 1652. ; 70% H2O; 15% CO2; 8% N2; 7% SO2
        ;Rm = 385.  ; 70% H2O; 15% CO2; 8% N2; 7% SO2
        g = 8.4
        Ca = 835.
        Ra = 191.
        TP = atmos_venus_vira(z0,g)
        End

  'venus-low':  Begin
        ;Cm = 1652. ; 70% H2O; 15% CO2; 8% N2; 7% SO2
        ;Rm = 385.  ; 70% H2O; 15% CO2; 8% N2; 7% SO2
        g = 8.4
        Ca = 835.
        Ra = 191.
        TP = atmos_venus_low(z0,g)
        End

  'venus-high':  Begin
        ;Cm = 1652. ; 70% H2O; 15% CO2; 8% N2; 7% SO2
        ;Rm = 385.  ; 70% H2O; 15% CO2; 8% N2; 7% SO2
        g = 8.4
        Ca = 835.
        Ra = 191.
        TP = atmos_venus_high(z0,g)
        End

  'venus-vega':  Begin
        ;Cm = 1652. ; 70% H2O; 15% CO2; 8% N2; 7% SO2
        ;Rm = 385.  ; 70% H2O; 15% CO2; 8% N2; 7% SO2
        g = 8.4
        Ca = 835.
        Ra = 191.
        TP = atmos_venus_vega2(z0,g)
        End

	'mars':		Begin
				;Cm = 1400.	; same as atmos
				;Rm = 350.	; same as atmos
				g = 3.7
				Ca = 835.
				Ra = 191.
				TP = atmos_mars(z0,g)
				End

	'paleomars':Begin
				g = 3.7
				Ca = 835.
				Ra = 191.
				TP = atmos_paleomars(z0,g)
				End

	'io':		Begin
				Cm = 1317.	; 50/50 SO2 & water vapor
				Rm = 296.	; 50/50 SO2 & water vapor
				g = 1.8
				Ca = 634.
				Ra = 130.
				TP = atmos_io(z0,g)
				End

EndCase

P0 = double(TP(1))
Ta0 = TP(0)
z = z0

rhom0 = P0/(Rm*double(theta0))
rhoe0 = P0/(Ra*double(theta0))
rhoa0 = P0/(Ra*double(Ta0))

phis0 = double(((1-n0)*rhom0)/(n0*rhos+(1-n0)*rhom0))
phie0 = double(0)
phim0 = double(1-phis0-phie0)

rhoB0 = rhom0*phim0 + rhoe0*phie0 + rhos*phis0

MsF = double(rhos*phis0*u0*b0)
MmF = double(rhom0*phim0*u0*b0)
MeF = double(rhoe0*phie0*u0*b0)

m0 = double(MsF + MmF + MeF)
CB0 = (MsF*Cs + MmF*Cm + MeF*Ca)/m0

nuF = double(m0*u0)
TF = double(m0*CB0*theta0)

y = dblarr(3)

y(0) = MeF
y(1) = nuF
y(2) = TF

defs = {z:z, y:y}
return, defs
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro lin_plumes, ifile

; This is the core of the code for solving the system of
; three differential equations describing the thermal
; dynamics of buoyant eruption columns from linear vents. 
; In this case, the plumes are comprised of magmatic gas 
; and solids. The ambient atmosphere may have a different 
; composition to the magmatic gas.  The program takes as
; input the initial values for gas content, plume half-width,
; bulk column velocity, and column temperature, and solves
; three differnetial equations, for plume half-width,
; velocity and temperature as functions of z, using the
; Runge-Kutta method. This IDL code has been translated
; from the IDL code 'plumes' for central vent eruptions.
; That code was translated from the fortran code 'mars.f'
; on puuoo.  The mars.f fortran code was modified from the
; venus.f code, which worked fine as of 5/4/94.  As of 12/18/00,
; the basic code for calculating r, u, and theta as 
; functions of z works for conditions on Earth. The code
; is also capable of detecting the plume top.  As of 12/29/00,
; the code is working fine for conditions on Earth, Venus and Mars.
; Modifications for linear plumes was made 10/29/08.
;
; The variables and constants used in the program are:
; g 		gravity [m s^(-2)] = 9.8 e; 3.7 m; 8.4 v
; rhoa0	ambient density 0 km above mean radius, Ta = T0, P = P0)
; alpha	entrainment constant = 0.09 (Woods, 1988; Morton
;		    et al., 1956)
; Ca	  specific heat for ambient [J K^(-1) kg^(-1)] = 998 e, 835 m,
;		835 v (m & v are a guess based on a mainly CO2 atmosphere)
; Cm	  specific heat of the magmatic gas = Ca, for
;		similar composition, = 2000 J K^(-1) kg^(-1) for
;		water vapor
; Cs	specific heat of ash particles = 920 J K^(-1) kg^(-1)
;		(Riehle, 1973)
; rhog	density of the plume gas, defined by the ideal gas law
; rhos	density of pumice = 1000 kg m^(-3) (Wilson, 1976)
; Ra	gas constant for ambient [J K^(-1) kg^(-1)], 287 e, 191 m, 191 v
;		(m & v are a guess based on a mainly CO2 atmosphere)
; Rm	gas constant for magmatic gas = Ra for similar
;		composition, = 461 J K^(-1) kg^(-1) for water vapor
; T0	temperature of atmosphere 0 km above mean radius [K]
; P0	pressure of atmospehre 0 km above mean radius [Pa]
; delz	step size increment = 5 m
; phi	volume fraction of solids in the control volume

common block,Ra,Rm,Ca,Cm,Cs,rhos,MsF,MmF,delz,g,alpha,atmos

defs = define(ifile)

top = 0
if delz le 1.0 then u_cutoff = 0.1 else u_cutoff = 10
z = defs.z
y = defs.y
outpro,z,y

repeat begin

	dydz = plumederivs(z, y)
	yout = rk4(y, dydz, z, delz, 'plumederivs',/double)
	top = toptest(z, y, yout, u_cutoff)
	z = z + delz
	y = yout
	outpro,z,y

endrep until (top eq 1)

print,'Plume Height (km):  ',z/1000

close,1
close,2

end