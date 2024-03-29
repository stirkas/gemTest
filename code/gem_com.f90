module gem_com

  !common data used for gem

  use mpi
  use gem_pputil
  use iso_c_binding

  implicit none

  INTERFACE
     real function revers(num,n)
     end function revers

     real function ran2(i)
     end function ran2

     real function en3(s)
       real :: s
     end function en3
  END INTERFACE

  integer, bind(c) :: imx,jmx,kmx,mmx,mmxe,nmx,nsmx,nsubd=8,&
       modemx,ntube,nxpp,ngdx=5,nb=6, &
       negrd=8,nlgrd=8

  character(len=70) outname
  REAL :: endtm,begtm,pstm
  REAL :: starttm,lasttm,tottm
  real :: aux1(50000),aux2(20000)
  real,dimension(:),allocatable :: workx,worky,workz
  complex,dimension(:),allocatable :: tmpx
  complex,dimension(:),allocatable :: tmpy
  complex,dimension(:),allocatable :: tmpz

  !          imx,jmx,kmx = max no. of grid pts in x,y,z
  !          mmx         = max no. of particles
  !          nmx         = max. no. of time steps
  !          nsmx        = max. no. of species (including tracer particles
  !                                             as a species)
  !          ntrmx       = max. no. of tracer particles
  !          modemx      = max. no. of mode history plots

  integer :: mme,mmb
  REAL, dimension(:,:),allocatable :: rwx,rwy
  INTEGER,target,dimension(:),allocatable :: mm,tmm,lr
  type(c_ptr), bind(c) :: mm_ptr, lr_ptr, tmm_ptr

  integer :: micell,mecell !jycheng
  integer :: nonlin1,nonlin2 !jycheng
  real :: mims3,q3 !jycheng
  REAL, target, dimension(:),allocatable :: tets,mims,q
  REAL,dimension(:),allocatable :: kapn, kapt
  type(c_ptr), bind(c) :: mims_ptr, q_ptr

  INTEGER :: timestep,im,jm,km,mykm,iseed,nrst,nfreq,isft,mynf,ifskp,iphbf,iapbf,idpbf
  real,dimension(:),allocatable :: time
  REAL, bind(c):: dx,dy,dz,pi,pi2,dt,dte,totvol,n0,n0e,tcurr,rmpp,rmaa,eprs
  REAL, bind(c) :: lx,ly,lz,xshape,yshape,zshape,pzcrit(5),pzcrite,encrit,tot_field_e,tot_joule,tot_joule1
  INTEGER, bind(C) :: nm,nsm,kcnt,jcnt,ncurr,llk,mlk,onemd,iflr,iorb
  integer :: izonal,adiabatic_electron,ineq0,iflut,nlow,ntor0,mstart
  REAL, bind(c):: cut,amp,tor,amie,isg,rneu,rneui,emass,qel,mbeam,qbeam,teth,vexbsw,vparsw
  REAL :: c4,fradi,kxcut,kycut,bcut,ftrap,adwn,adwe,adwp,frmax
  INTEGER, bind(c) :: iput,iget,idg,kzlook,ision,isiap,peritr,iadi,ipred,icorr,jpred,jcorr
  REAL,target, DIMENSION(:,:),allocatable :: yyamp,yyre,yyim
  complex,dimension(:,:),allocatable :: camp,campf
  REAL, bind(c) :: br0,lr0,qp,width,e0,vwidth,vwidthe,vcut,vpp,vt0,yd0
  integer, bind(c) :: nonlin(5),nonline,ipara,isuni,ifluid,ishift,nopz,nopi(5),noen,nowe
  complex :: IU
  real,dimension(:),allocatable :: coefx,coefy,coefz
  complex,dimension(1:8) :: apk,ptk,dpdtk
  integer,dimension(1:8) :: lapa,mapa,napa
  real :: mrtio(0:1),aven,avptch
  integer :: icrs_sec,ipg,isphi
  integer,dimension(0:255) :: isgnft,jft

  !          im,jm,km = max no. of grid pts in x,y,z
  !          mm       = max no. of particles
  !          nm       = max. no. of time steps
  !          nsm      = max. no. of species (including the set of
  !                           tracer particles as a species)
  !          ntrm     = max. no. of tracer particles
  !          modem    = max. no. of mode history plots
  !          iput     = 1 save for restart into dump.b, =0 no save
  !          iget     = 1 restart from dump.b, =0 use loader

  !            field or grid quantities

  REAL,DIMENSION(:,:,:,:),allocatable :: den
  REAL,DIMENSION(:,:,:,:),allocatable :: dnidt,jpar,jpex,jpey,dti
  REAL,DIMENSION(:,:,:),allocatable :: rho,jion,jionx,jiony
  real,dimension(:,:,:),allocatable :: phi
  real,dimension(:,:,:),allocatable :: drhodt,dnedt,dphidt,drhoidt
  REAL,target,DIMENSION(:,:,:),allocatable :: ex
  REAL,target, DIMENSION(:,:,:),allocatable :: ey
  REAL,target, DIMENSION(:,:,:),allocatable :: ez
  REAL,target, DIMENSION(:,:,:),allocatable :: dpdz,dadz
  REAL,target, DIMENSION(:,:,:),allocatable :: delbx,delby
  REAL,DIMENSION(:),allocatable :: xg,yg,zg
  real,target, dimension(:,:,:),allocatable :: apar,dene
  real,dimension(:,:,:),allocatable :: upar,upart,delte
  real,dimension(:,:,:),allocatable :: upex,upey,upa0,den0,upazd,upa00,upa0t,den0apa
  real,dimension(:,:),allocatable :: cfx,cfy,jac,bmag,bdgxcgy,bdgrzn,ggxdgy,ggy2,ggx
  real,target, dimension(:),allocatable :: gn0e,gt0e,gt0i,avap 
  real,dimension(:,:),allocatable :: gn0s
  type(c_ptr), bind(c) :: phi_ptr, ex_ptr, ey_ptr, ez_ptr, dpdz_ptr, dadz_ptr, delbx_ptr, delby_ptr, apar_ptr
  !          particle array declarations
  REAL,target,DIMENSION(:,:),allocatable :: mu,xii,pzi,eki,z0i,u0i
  REAL,target,DIMENSION(:,:),allocatable :: x2,y2,z2,u2
  REAL,target,DIMENSION(:,:),allocatable :: x3,y3,z3,u3
  REAL,target,DIMENSION(:,:),allocatable :: w2,w3
  type(c_ptr), bind(c) :: x2_ptr, z2_ptr, u2_ptr, mu_ptr, y2_ptr, w3_ptr, w2_ptr, pzi_ptr, xii_ptr, z0i_ptr, uoi_ptr, u3_ptr, x3_ptr, y3_ptr, z3_ptr, eki_ptr


  REAL,DIMENSION(:),allocatable :: mue,xie,pze,eke,z0e,u0e
  REAL,DIMENSION(:),allocatable :: x2e,y2e,z2e,u2e,mue2
  REAL,DIMENSION(:),allocatable :: x3e,y3e,z3e,u3e,mue3
  REAL,DIMENSION(:),allocatable :: w2e,w3e
  real,dimension(:),allocatable :: ipass, index
  REAL,DIMENSION(:),allocatable :: w000,w001,w010,w011,w100,w101,w110,w111

  !              Various diagnostic arrays and scalars
  !    plotting constants

  INTEGER, bind(c) :: nplot,xnplt,imovie=1000000,nzcrt,npze,npzi,npzc,npzb
  REAL :: contu,wmax

  !    energy diagnostic arrays

  REAL, target, DIMENSION(:,:),allocatable :: ke
  REAL,DIMENSION(:),allocatable :: fe,te
  REAL,target, DIMENSION(:),allocatable :: rmsphi,rmsapa,avewe
  REAL,target, DIMENSION(:,:),allocatable :: nos,avewi
  type(c_ptr), bind(c) :: avewi_ptr, ke_ptr, nos_ptr

  !    flux diagnostics
  REAL,target, DIMENSION(:),allocatable :: vol
  REAL,target, DIMENSION(:,:),allocatable :: efle_es,efle_em,pfle_es,pfle_em
  REAL,target, DIMENSION(:,:,:),allocatable :: pfl_es,pfl_em,efl_es,efl_em
  REAL,DIMENSION(:,:),allocatable :: chii, chie, ddi
  REAL,DIMENSION(:),allocatable :: achii, achie, addi
  type(c_ptr), bind(c) :: efle_es_prt, pfle_es_ptr
  !   mode diagnositics
  INTEGER :: modem
  INTEGER,dimension(:),allocatable :: lmode,mmode,nmode
  complex,dimension(:,:),allocatable :: pmodehis
  real,target, dimension(:),allocatable :: mdhis,mdhisa,mdhisb,mdhisc,mdhisd
  complex,target, dimension(:,:),allocatable :: aparhis,phihis

  !   kr, ktheta spectrum plots
  REAL,DIMENSION(:,:),allocatable :: phik

  !     weighty variables
  INTEGER,dimension(:),allocatable :: deljp,deljm
  INTEGER,dimension(:,:),allocatable :: jpl
  INTEGER,dimension(:,:),allocatable :: jpn
  INTEGER,dimension(:,:),allocatable :: jmi
  INTEGER,dimension(:,:),allocatable :: jmn
  REAL,DIMENSION(:),allocatable :: weightp,weightm
  REAL,DIMENSION(:),allocatable :: weightpn,weightmn

  !blending variable
  complex,dimension(:,:,:,:),allocatable :: pol,pmtrx,pmtrxi
  complex,dimension(:,:),allocatable :: pfac

  !      MPI variables
  !  include '/usr/include/mpif.h'

  integer, bind(c) :: Master=0!,parameter
  integer :: numprocs
  INTEGER, bind(c):: MyId,Last,cnt,ierr
  INTEGER :: GRID_COMM,TUBE_COMM
  INTEGER, bind(c) :: GCLR,TCLR,GLST,TLST
  INTEGER :: stat(MPI_STATUS_SIZE)
  INTEGER :: lngbr,rngbr,idprv,idnxt

  character(len=*) directory
  parameter(directory='./dump/')

  character(len=*) outdir
  parameter(outdir='./out/')

  !real :: ran2,revers
  !integer :: mod
  !real :: amod
  save

  !this is the pointer definition for ()
  type(c_ptr), bind(c) :: vol_ptr,rmsapa_ptr,rmsphi_ptr,avewe_ptr,mdhis_ptr,mdhisa_ptr
  type(c_ptr), bind(c) :: mdhisb_ptr,mdhisc_ptr,mdhisd_ptr
  
  type(c_ptr), bind(c) :: efle_es_ptr,yyre_ptr,yyim_ptr
  type(c_ptr), bind(c) :: yyamp_ptr, efle_em_ptr, pfle_em_ptr

  type(c_ptr), bind(c) :: phihis_ptr,aparhis_ptr

  type(c_ptr), bind(c) :: pfl_es_ptr, efl_es_ptr, efl_em_ptr, pfl_em_ptr


contains

  subroutine new_gem_com()
    nxpp = imx !/ntube
    allocate(workx(4*imx),worky(4*jmx),workz(4*kmx))
    allocate(tmpx(0:imx-1))
    allocate(tmpy(0:jmx-1))
    allocate(tmpz(0:kmx-1))

    allocate(rwx(nsmx,4),rwy(nsmx,4))
    allocate(mm(nsmx),tmm(nsmx),lr(nsmx))
    allocate(tets(nsmx),mims(nsmx),q(nsmx))
    allocate(kapn(nsmx),kapt(nsmx))
    allocate(time(0:nmx))
    allocate(yyamp(jmx,0:4),yyre(jmx,0:4),yyim(jmx,0:4),camp(0:6,0:50000),campf(0:6,0:nfreq-1))
    allocate(aparhis(0:6,0:jcnt-1),phihis(0:6,0:jcnt-1))
    allocate(mdhis(0:100),mdhisa(0:100),mdhisb(0:100))
    allocate(mdhisc(0:100),mdhisd(0:100))
    allocate(coefx(100+8*imx),coefy(100+8*jmx),coefz(100+8*kmx))

    ALLOCATE( den(nsmx,0:nxpp,0:jmx,0:1),dti(nsmx,0:nxpp,0:jmx,0:1), &
         delte(0:nxpp,0:jmx,0:1))
    ALLOCATE( rho(0:nxpp,0:jmx,0:1),drhoidt(0:nxpp,0:jmx,0:1), &
         jion(0:nxpp,0:jmx,0:1),jionx(0:nxpp,0:jmx,0:1), &
         jiony(0:nxpp,0:jmx,0:1))
    allocate( phi(0:nxpp,0:jmx,0:1))
    allocate( drhodt(0:nxpp,0:jmx,0:1),dnedt(0:nxpp,0:jmx,0:1))
    allocate( dnidt(nsmx,0:nxpp,0:jmx,0:1),jpar(nsmx,0:nxpp,0:jmx,0:1),  &
         jpex(nsmx,0:nxpp,0:jmx,0:1),jpey(nsmx,0:nxpp,0:jmx,0:1))
    allocate( dphidt(0:nxpp,0:jmx,0:1))
    ALLOCATE( ex(0:nxpp,0:jmx,0:1)) 
    ALLOCATE( ey(0:nxpp,0:jmx,0:1)) 
    ALLOCATE( ez(0:nxpp,0:jmx,0:1))
    ALLOCATE( dpdz(0:nxpp,0:jmx,0:1),dadz(0:nxpp,0:jmx,0:1))
    ALLOCATE( delbx(0:nxpp,0:jmx,0:1),delby(0:nxpp,0:jmx,0:1))
    ALLOCATE( xg(0:nxpp),yg(0:jmx),zg(0:1))
    allocate( apar(0:nxpp,0:jmx,0:1),dene(0:nxpp,0:jmx,0:1))
    allocate( upar(0:nxpp,0:jmx,0:1),upart(0:nxpp,0:jmx,0:1),upex(0:nxpp,0:jmx,0:1), &
         upey(0:nxpp,0:jmx,0:1),upa0(0:nxpp,0:jmx,0:1), &
         den0(0:nxpp,0:jmx,0:1),upazd(0:nxpp,0:jmx,0:1),&
         upa00(0:nxpp,0:jmx,0:1),upa0t(0:nxpp,0:jmx,0:1),den0apa(0:nxpp,0:jmx,0:1))
    allocate( cfx(0:nxpp,0:1),cfy(0:nxpp,0:1),jac(0:nxpp,0:1))
    allocate( bmag(0:nxpp,0:1),bdgxcgy(0:nxpp,0:1),bdgrzn(0:nxpp,0:1))
    allocate( ggxdgy(0:nxpp,0:1),ggy2(0:nxpp,0:1),ggx(0:nxpp,0:1))
    allocate (gn0e(0:nxpp),gt0e(0:nxpp),gt0i(0:nxpp),avap(0:nxpp))
    allocate (gn0s(1:5,0:nxpp))
    !          particle array declarations
    allocate( mu(nsmx,1:mmx),xii(nsmx,1:mmx),pzi(nsmx,1:mmx), &
         eki(nsmx,1:mmx),z0i(nsmx,1:mmx),u0i(nsmx,1:mmx))
    allocate( x2(nsmx,1:mmx),y2(nsmx,1:mmx),z2(nsmx,1:mmx),u2(nsmx,1:mmx))
    allocate( x3(nsmx,1:mmx),y3(nsmx,1:mmx),z3(nsmx,1:mmx),u3(nsmx,1:mmx))
    allocate( w2(nsmx,1:mmx),w3(nsmx,1:mmx))

    allocate( mue(1:mmxe),xie(1:mmxe),pze(1:mmxe),eke(1:mmxe),z0e(1:mmxe),u0e(1:mmxe))
    allocate( x2e(1:mmxe),y2e(1:mmxe),z2e(1:mmxe),u2e(1:mmxe),mue2(1:mmxe))
    allocate( x3e(1:mmxe),y3e(1:mmxe),z3e(1:mmxe),u3e(1:mmxe),mue3(1:mmxe))
    allocate( w2e(1:mmxe),w3e(1:mmxe))
    allocate( ipass(1:mmxe), index(1:mmxe))
    allocate(w000(1:mmxe),w001(1:mmxe),w010(1:mmxe),w011(1:mmxe),&
         w100(1:mmxe),w101(1:mmxe),w110(1:mmxe),w111(1:mmxe))

    !              Various diagnostic arrays and scalars
    !    plotting constants

    !    energy diagnostic arrays

    ALLOCATE( ke(nsmx,0:nmx),fe(0:nmx),te(0:nmx))
    ALLOCATE( rmsphi(0:nmx),rmsapa(0:nmx),avewi(1:3,0:nmx),avewe(0:nmx))
    ALLOCATE( nos(nsmx,0:nmx))

    !    flux diagnostics
    ALLOCATE( vol(1:nsubd),efle_es(1:nsubd,0:nmx),pfle_es(1:nsubd,0:nmx), &
         pfl_es(1:3,1:nsubd,0:nmx),efl_es(1:3,1:nsubd,0:nmx), &
         pfle_em(1:nsubd,0:nmx),efle_em(1:nsubd,0:nmx), &
         pfl_em(1:3,1:nsubd,0:nmx),efl_em(1:3,1:nsubd,0:nmx))
    ALLOCATE( chii(1:nsubd,0:nmx),chie(1:nsubd,0:nmx),ddi(1:nsubd,0:nmx), &
         achii(1:nsubd),achie(1:nsubd),addi(1:nsubd))

    !   mode diagnositics
    allocate(lmode(modemx),mmode(modemx),nmode(modemx))
    allocate(pmodehis(modemx,0:nmx))

    !   kr, ktheta spectrum plots
    ALLOCATE( phik(imx,jmx))

    !     weighty variables
    ALLOCATE( deljp(0:nxpp),deljm(0:nxpp))
    ALLOCATE( jpl(0:nxpp,0:jmx))
    ALLOCATE( jpn(0:nxpp,0:jmx))
    ALLOCATE( jmi(0:nxpp,0:jmx))
    ALLOCATE( jmn(0:nxpp,0:jmx))
    ALLOCATE( weightp(0:nxpp),weightm(0:nxpp))
    ALLOCATE( weightpn(0:nxpp),weightmn(0:nxpp))

    !Blending variable
    ALLOCATE(pol(1:nb,0:imx-1,0:jmx-1,0:kmx),pfac(0:imx-1,0:jmx-1), &
         pmtrx(0:imx-1,0:jmx-1,1:nb,1:nb), &
         pmtrxi(0:imx-1,0:jmx-1,1:nb,1:nb))

     !spec 1d arrays
     vol_ptr    = c_loc(vol(1))
     rmsapa_ptr = c_loc(rmsapa(0))
     rmsphi_ptr = c_loc(rmsphi(0))
     avewe_ptr  = c_loc(avewe(0))
     mdhis_ptr  = c_loc(mdhis(0))
     mdhisa_ptr = c_loc(mdhisa(1))
     mdhisb_ptr = c_loc(mdhisb(1))
     mdhisc_ptr = c_loc(mdhisc(1))
     mdhisd_ptr = c_loc(mdhisd(1))
     lr_ptr     = c_loc(lr(1))
     mm_ptr     = c_loc(mm(1))
     tmm_ptr    = c_loc(tmm(1))
     mims_ptr   = c_loc(mims(1))

     !spec 2d arrays
     pfle_es_ptr = c_loc(pfle_es(1,0))
     avewi_ptr   = c_loc(avewi(1,0))
     efle_es_ptr = c_loc(efle_es(1,0))
     yyre_ptr    = c_loc(yyre(1,0))
     yyim_ptr    = c_loc(yyim(1,0))
     yyamp_ptr   = c_loc(yyamp(1,0))
     phihis_ptr  = c_loc(phihis(0,0))
     aparhis_ptr = c_loc(aparhis(0,0))
     efle_em_ptr = c_loc(efle_em(1,0))
     pfle_em_ptr = c_loc(pfle_em(1,0))

     x2_ptr  = c_loc(x2(1,1))
     z2_ptr  = c_loc(z2(1,1))
     u2_ptr  = c_loc(u2(1,1))
     mu_ptr  = c_loc(mu(1,1))
     y2_ptr  = c_loc(y2(1,1))
     w3_ptr  = c_loc(w3(1,1))
     w2_ptr  = c_loc(w2(1,1))
     pzi_ptr = c_loc(pzi(1,1))
     xii_ptr = c_loc(xii(1,1))
     uoi_ptr = c_loc(u0i(1,1))
     u3_ptr  = c_loc(u3(1,1))
     x3_ptr  = c_loc(x3(1,1))
     y3_ptr  = c_loc(y3(1,1))
     z3_ptr  = c_loc(z3(1,1))
     eki_ptr = c_loc(eki(1,1))
     nos_ptr = c_loc(nos(1,0))
     ke_ptr  = c_loc(ke(1,0))




     !spec 3d arrays
     pfl_es_ptr = c_loc(pfl_es(1,1,0))
     efl_es_ptr = c_loc(efl_es(1,1,0))
     efl_em_ptr = c_loc(efl_em(1,1,0))
     pfl_em_ptr = c_loc(pfl_em(1,1,0))

     ex_ptr = c_loc(ex(0,0,0))
     ey_ptr = c_loc(ey(0,0,0))
     ez_ptr = c_loc(ez(0,0,0))

     delbx_ptr = c_loc(delbx(0,0,0))
     delby_ptr = c_loc(delby(0,0,0))
     dpdz_ptr  = c_loc(dpdz(0,0,0))
     dadz_ptr  = c_loc(dadz(0,0,0))
     apar_ptr  = c_loc(ex(0,0,0))






     call new_gem_com_c()

  end subroutine new_gem_com

end module gem_com










