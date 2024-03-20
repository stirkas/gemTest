import os

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

import scriptUtils as su
import pydiag.utils.comm as cm
import pydiag.utils.geom as gm

def gemGeom():
    #Normalizations/definitions.
    pi = np.arctan(1.0)*4
    qref = 1.6e-19
    mref = 1.67e-27
    Bref = 2.0
    nref = 4.66e19
    Tref = 2.14*1000*qref #keV to Joules.
    oref = qref*Bref/mref
    vref = np.sqrt(Tref/mref)
    xref = mref*vref/(qref*Bref)

    beta = 4*pi*1e-7*nref*Tref/Bref**2
    teti = 1.0
    isgnf = 1
    isgnq = 1

    Lref  = 1.67/xref
    Rmaj0 = Lref
    a     = Rmaj0*0.36

    #Start making grid.
    r0a   = 0.5
    lxa   = 0.6
    rina  = r0a - lxa/2
    routa = r0a + lxa/2
    rin   = rina*a
    rout  = routa*a
    r0    = r0a*a
    nr  = 300
    nth = 100
    dr  = (rout-rin)/nr
    dth = 2*pi/nth

    #Compute Miller parameters.
    elon0 = 1.0  #Elongation: up-down stretch.
    tria0 = 0.0  #Triangularity: radial outward stretch.
    selon0 = 0.0 #d(elon)/dr
    stria0 = 0.0 #d(tria)/dr
    Rmaj0p = 0.0 #Shafranov shift.
    elonp0 = elon0*selon0/r0
    triap0 = np.sqrt(1-tria0**2)*stria0/r0
    q0    = 1.41
    shat0 = 0.78
    f     = np.zeros(nr+1)
    Rmaj  = np.zeros(nr+1)
    Rmajp = np.zeros(nr+1)
    elon  = np.zeros(nr+1)
    selon = np.zeros(nr+1)
    tria  = np.zeros(nr+1)
    stria = np.zeros(nr+1)
    for i in range(nr+1):
        r    = rin + i*dr
        f[i] = isgnf*Rmaj0
        Rmaj[i]  = Rmaj0 - (Rmaj0p/2)*(1-(r/a)**2)*a
        Rmajp[i] = Rmaj0p*(r/a)
        elon[i]  = (elon0-0.1) + 0.1*(r/a)**4
        selon[i] = r*0.4*(r/a)**3/(a*elon[i])
        tria[i]  = tria0*(r/a)**2
        stria[i] = r*tria0*2*(r/a)/(a*np.sqrt(1-tria[i]**2))

    if (elon0 == 1):
        for i in range(nr+1):
            elon[i]  = 1.0
            selon[i] = 0.0

    #Compute R(r,th) and Z(r,th)
    rad  = np.zeros((nr+1,nth+1))
    hght = np.zeros((nr+1,nth+1))
    for i in range(nr+1):
        r = rin + i*dr
        for j in range(nth+1):
            th = -pi + dth*j
            rad[i,j]  = Rmaj[i] + r*np.cos(th + np.arcsin(tria[i])*np.sin(th))
            hght[i,j] = r*elon[i]*np.sin(th)

    #Compute grad(r),grad(th),grdgt,grcgt
    elonp = np.zeros(nr+1)
    triap = np.zeros(nr+1)
    gr    = np.zeros((nr+1,nth+1))
    gth   = np.zeros((nr+1,nth+1))
    grdgt = np.zeros((nr+1,nth+1))
    grcgt = np.zeros((nr+1,nth+1))
    for i in range(nr+1):
        r  = rin + i*dr
        x  = np.arcsin(tria[i])
        xp = stria[i]/r
        elonp[i] = selon[i]*elon[i]/r
        triap[i] = stria[i]*np.sqrt(1-tria[i]**2)/r
        for j in range(nth+1):
            th = -pi + j*dth
            c1 = elon[i]*r*np.cos(th)                        #dZ/dth
            c2 = -r*np.sin(th+x*np.sin(th))*(1+x*np.cos(th)) #dR/dth
            c3 = elonp[i]*r*np.sin(th)+elon[i]*np.sin(th)    #dZ/dr
            c4 = Rmajp[i]+np.cos(th+x*np.sin(th))-r*np.sin(th+x*np.sin(th))*xp*np.sin(th) #dR/dr
            denom = c1*c4-c2*c3
            gr[i,j]  = np.sqrt(c1**2+c2**2)/denom
            gth[i,j] = np.sqrt(c3**2+c4**2)/denom
            grdgt[i,j] = (-c1*c3-c2*c4)/denom**2
            grcgt[i,j] = 1/denom

    #Make q,T,n profiles.
    kappat = 6.96
    kappan = 2.23
    wt = wn = 0.3
    q   = np.zeros(nr+1)
    t0i = np.zeros(nr+1)
    t0e = np.zeros(nr+1)
    n0i = np.zeros(nr+1)
    n0e = np.zeros(nr+1)
    for i in range(nr+1):
        r = rin + i*dr
        s = r/a
        q[i] = (2.52*s**2 - 0.16*s + 0.86)*isgnq
        t0i[i] = np.exp(-(kappat*wt*a/Lref)*np.tanh((r-r0)/(wt*a)))
        t0e[i] = t0i[i]
        n0i[i] = np.exp(-(kappan*wn*a/Lref)*np.tanh((r-r0)/(wn*a)))
        n0e[i] = n0i[i]
    q0    = q[nr//2]
    q0abs = np.abs(q0)

    #Compute psi(r) and psip(r)
    psi  = np.zeros(nr+1)
    psip = np.zeros(nr+1)
    for i in range(nr+1):
        var = 0
        for j in range(nth):
            var += dth/(rad[i,j]*grcgt[i,j])
        psip[i] = f[i]*var/(2*pi*q[i])
    for i in range(1,nr+1):
        psi[i] = psi[i-1] + dr*(psip[i-1] + psip[i])/2

    #Compute B(r,th) and J(r,th)
    bfld = np.zeros((nr+1,nth+1))
    jac  = np.zeros((nr+1,nth+1))
    for i in range(nr+1):
        r = rin + i*dr
        x = np.arcsin(tria[i])
        xp = stria[i]/r
        for j in range(0,nth+1):
            th = -pi + j*dth
            bfld[i,j] = np.sqrt((f[i]/rad[i,j])**2 + (psip[i]*gr[i,j]/rad[i,j])**2)
            jac[i,j]  = 1/(r0*Rmaj0*grcgt[i,j]/rad[i,j])

    #fig1 = plt.figure()
    #rArr  = np.linspace(rin,rout,nr+1)
    #thArr = np.linspace(-np.pi,np.pi,nth+1)
    #plt.plot(thArr,jac[nr//2,:])
    #plt.savefig('jz.pdf')
    #fig2 = plt.figure()
    #plt.plot(thArr,bfld[nr//2,:])
    #plt.savefig('bfld.pdf')
    
    return bfld,jac

def gaulag(x1, x2, n):
    """
    Compute the knots and weights for Gauss-Laguerre quadrature of order n
    on the interval [x1, x2].
    """
    import numpy as np
    
    hn = 1.0/n
    x = np.zeros(n)
    w = np.zeros(n)
    
    for nr in range(1, n+1):
        z = hn
        if nr > 1:
            z = x[nr-2] + hn*(nr**1.27)
        it = 0
        while True:
            it += 1
            z0 = z
            p = 1.0
            for i in range(nr-1):
                p *= (z - x[i])
            f0 = 1.0
            f1 = 1.0 - z
            for k in range(2, n+1):
                pf = ((2.0*k - 1.0 - z)*f1 - (k - 1.0)*f0)/k
                pd = k/z*(pf - f1)
                f0 = f1
                f1 = pf
            fd = pf/p
            q = 0.0
            for i in range(nr-1):
                wp = 1.0
                for j in range(nr-1):
                    if j != i:
                        wp *= (z - x[j])
                q += wp
            gd = (pd - q*fd)/p
            z = z - fd/gd
            if it <= 40 and abs((z - z0)/z) > 1.0E-15:
                continue
            x[nr-1] = z
            w[nr-1] = 1.0/(z*pd*pd)
            break
    
    w = w*np.exp(x)
    fac = x2/np.sum(w)
    x = x*fac
    w = w*fac
    
    return x, w

nz = 16
nw = 8
nv = 32
lv = 3
lw = 9
dv = 2*lv/(nv-1)
dz = 2*np.pi/nz
T0i = T0e = 1
#Get Gauss-Laguerre weights for mu integration.
knots,weights = gaulag(0,lw,nw)

#z-int first
gamzi  = np.zeros((nw,nv,nz))
gamze  = np.zeros((nw,nv,nz))
heatzi = np.zeros((nw,nv,nz))
heatze = np.zeros((nw,nv,nz))
gami   = np.zeros((nw,nv))
game   = np.zeros((nw,nv))
heati  = np.zeros((nw,nv))
heate  = np.zeros((nw,nv))
int1gami  = np.zeros(nv)
int1game  = np.zeros(nv)
int1heati = np.zeros(nv)
int1heate = np.zeros(nv)

geneIn_ge  = np.zeros((nw,nv))
geneIn_qe  = np.zeros((nw,nv))
geneInz_ge = np.zeros((nw,nv,nz))
geneInz_qe = np.zeros((nw,nv,nz))
geneIn_gi  = np.zeros((nw,nv))
geneIn_qi  = np.zeros((nw,nv))
geneInz_gi = np.zeros((nw,nv,nz))
geneInz_qi = np.zeros((nw,nv,nz))
gamze2     = np.zeros((nw,nv,nz))

#Read 3d vsp data files (mu,vpar,z). Both fluxes for both species.
for z in range(nz):
    filename = "geneFiles/el/gamma/gene" + str(z) + ".in"
    f = open(filename, "r")

    i = 0
    for line in f:
        if (line[0] == "#"):
            continue
        else:
            gamze[:,i,z] = line.split()[:]
            i += 1
    f.close()
    
    filename = "geneFiles/el/qheat/gene" + str(z) + ".in"
    f = open(filename, "r")

    i = 0
    for line in f:
        if (line[0] == "#"):
            continue
        else:
            heatze[:,i,z] = line.split()[:]
            i += 1
    f.close()

    filename = "geneFiles/ion/gamma/gene" + str(z) + "_gi.in"
    f = open(filename, "r")

    i = 0
    for line in f:
        if (line[0] == "#"):
            continue
        else:
            gamzi[:,i,z] = line.split()[:]
            i += 1
    f.close()

    filename = "geneFiles/ion/qheat/gene" + str(z) + "_qi.in"
    f = open(filename, "r")

    i = 0
    for line in f:
        if (line[0] == "#"):
            continue
        else:
            heatzi[:,i,z] = line.split()[:]
            i += 1
    f.close()

    #f = open('gene.in', 'r')

    #i = 0
    #for line in f:
    #    if (line[0] == "#"):
    #        continue
    #    else:
    #        gamze2[:,i,0] = line.split()[:]
    #        i += 1
    #f.close()

origDir = os.getcwd()
os.chdir('geneFiles')
common = cm.CommonData('_0008',-1,-2)
geomData = gm.Geometry(common)
jac = su.getJacobian(common)
os.chdir(origDir)
B0 = geomData.Bfield
jacInt = np.sum(jac*dz) #Don't integrate with trapz. Less accuracy but same as GENE.
vArea = np.sum(weights)*6

gem = 1
f0flag = 0

print("GENE integration of G and Q at r/a=0.5 only.")
#Use gamma to check q
for w in range(nw):
    for v in range(nv):
        for z in range(nz):
            #Remove vspace integration factors from data since they are pre-included.
            gamzi[w,v,z] /= (B0[z]*np.pi*weights[w]*dv) 
            gamze[w,v,z] /= (B0[z]*np.pi*weights[w]*dv) 
            heatzi[w,v,z] /= (B0[z]*np.pi*weights[w]*dv)
            heatze[w,v,z] /= (B0[z]*np.pi*weights[w]*dv)
            
            #Store temporarily for z-avging.
            geneInz_ge[w,v,z] = gamze[w,v,z]
            geneInz_qe[w,v,z] = heatze[w,v,z]
            geneInz_gi[w,v,z] = gamzi[w,v,z]
            geneInz_qi[w,v,z] = heatzi[w,v,z]

            #Try to remove average gamma so what is input into gem doesn't change particle fluxes for now. When just adding G_ES to GEM.
            #gamze[w,v,z] -= -0.5737773811447934/(B0[z]*np.pi*weights[w]*dv*nv*nw) 
            #gamzi[w,v,z] -= -0.5738052985178931/(B0[z]*np.pi*weights[w]*dv*nv*nw)

            #Use gamma to make q if wanted.
            vpar = -3 + dv*v
            heatzi[w,v,z] = gamzi[w,v,z]*(vpar**2 + knots[w]*B0[z])*T0i
            heatze[w,v,z] = gamze[w,v,z]*(vpar**2 + knots[w]*B0[z])*T0e

for w in range(nw):
    for v in range(nv):
        for z in range(nz):
            gamzi[w,v,z]  *= np.pi*B0[z]*weights[w]*dv*jac[z]*dz
            gamze[w,v,z]  *= np.pi*B0[z]*weights[w]*dv*jac[z]*dz
            heatzi[w,v,z] *= np.pi*B0[z]*weights[w]*dv*jac[z]*dz
            heatze[w,v,z] *= np.pi*B0[z]*weights[w]*dv*jac[z]*dz

gami  = np.sum(gamzi)/jacInt
game  = np.sum(gamze)/jacInt
heati = np.sum(heatzi)/jacInt
heate = np.sum(heatze)/jacInt

print("Gi_ES: " + str(gami))
print("Ge_ES: " + str(game))
print("Qi_ES: " + str(heati))
print("Qe_ES: " + str(heate))

#Output just z-average for GEM.
for z in range(nz):
    geneInz_ge[:,:,z] = geneInz_ge[:,:,z]*jac[z]*dz
    geneInz_qe[:,:,z] = geneInz_qe[:,:,z]*jac[z]*dz
    geneInz_gi[:,:,z] = geneInz_gi[:,:,z]*jac[z]*dz
    geneInz_qi[:,:,z] = geneInz_qi[:,:,z]*jac[z]*dz

geneIn_ge = np.sum(geneInz_ge,axis=2)/jacInt
geneIn_qe = np.sum(geneInz_qe,axis=2)/jacInt
geneIn_gi = np.sum(geneInz_gi,axis=2)/jacInt
geneIn_qi = np.sum(geneInz_qi,axis=2)/jacInt

for w in range(nw):
    for v in range(nv):
        for z in range(nz):
            vpar = -3 + v*dv
            gamzi[w,v,z]  = geneIn_gi[w,v]*np.pi*B0[z]*weights[w]*dv*jac[z]*dz
            gamze[w,v,z]  = geneIn_ge[w,v]*np.pi*B0[z]*weights[w]*dv*jac[z]*dz
            heatzi[w,v,z] = geneIn_gi[w,v]*np.pi*B0[z]*weights[w]*dv*jac[z]*dz*(vpar**2 + knots[w]*B0[z])
            heatze[w,v,z] = geneIn_ge[w,v]*np.pi*B0[z]*weights[w]*dv*jac[z]*dz*(vpar**2 + knots[w]*B0[z])

gami2  = np.sum(gamzi)/jacInt
game2  = np.sum(gamze)/jacInt
heati2 = np.sum(heatzi)/jacInt
heate2 = np.sum(heatze)/jacInt

print("Note: z-avg changes fluxes!!!")
print("<Gi_ES>_z: " + str(gami2))
print("<Ge_ES>_z: " + str(game2))
print("<Qi_ES>_z: " + str(heati2))
print("<Qe_ES>_z: " + str(heate2))

geneIn_ge *= (heate/heate2)
geneIn_qe *= (heate/heate2)
geneIn_gi *= (heate/heate2)
geneIn_qi *= (heate/heate2)

geneIn = np.transpose(geneIn_ge)
geneIn2 = np.transpose(geneIn_gi)
with open('gene.in', 'w') as file:
    file.truncate(0)
    for row in geneIn:
        file.write(' '.join([str(a) for a in row]) + '\n')
file.close()
with open('gene2.in', 'w') as file:
    file.truncate(0)
    for row in geneIn2:
        file.write(' '.join([str(a) for a in row]) + '\n')
file.close()

for w in range(nw):
    for v in range(nv):
        for z in range(nz):
            vpar = -3 + v*dv
            gamzi[w,v,z]  = geneIn_gi[w,v]*np.pi*B0[z]*weights[w]*dv*jac[z]*dz
            gamze[w,v,z]  = geneIn_ge[w,v]*np.pi*B0[z]*weights[w]*dv*jac[z]*dz
            heatzi[w,v,z] = geneIn_gi[w,v]*np.pi*B0[z]*weights[w]*dv*jac[z]*dz*(vpar**2 + knots[w]*B0[z])
            heatze[w,v,z] = geneIn_ge[w,v]*np.pi*B0[z]*weights[w]*dv*jac[z]*dz*(vpar**2 + knots[w]*B0[z])

gami3  = np.sum(gamzi)/jacInt
game3  = np.sum(gamze)/jacInt
heati3 = np.sum(heatzi)/jacInt
heate3 = np.sum(heatze)/jacInt

print("Note: z-avg changes fluxes!!!")
print("<Gi_ES>_z: " + str(gami3))
print("<Ge_ES>_z: " + str(game3))
print("<Qi_ES>_z: " + str(heati3))
print("<Qe_ES>_z: " + str(heate3))

vArr = np.linspace(-3,3,32)
interp = RegularGridInterpolator((knots,vArr),geneIn_ge)

pts = np.array([2.0140615,0.92003])
intp = interp(pts)
print(intp)