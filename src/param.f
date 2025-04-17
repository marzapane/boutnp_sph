c     parameter (m1=97,m2=101,m3=1025,ndd=2,ndv=3)
c     parameter (m1=257,m2=193,m3=1537,ndd=2,ndv=3)
c     parameter (m1=193,m2=129,m3=513,ndd=2,ndv=3)
c     parameter (m1=385,m2=257,m3=1025,ndd=2,ndv=3)
c     parameter (m1=1025,m2=129,m3=513,ndd=2,ndv=3)   ! huge
c     parameter (m1=769,m2=129,m3=581,ndd=2,ndv=3)    ! bigger
      parameter (m1=513,m2=65,m3=513,ndd=2,ndv=3)     ! big
c     parameter (m1=257,m2=33,m3=257,ndd=2,ndv=3)     ! small
c     parameter (m1=129,m2=17,m3=129,ndd=2,ndv=3)     ! little
c     parameter (m1=257,m2=193,m3=770,ndd=2,ndv=3)
c     parameter (m1=97,m2=129,m3=257,ndd=2,ndv=3)
c     parameter (m1=33,m2=33,m3=64,ndd=2,ndv=3)       ! solo per assuref
c     parameter (m1=385,m2=385,m3=385,ndd=2,ndv=3)
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      parameter (m3mh=m3m/2+1)
      parameter (m12=2*m1m,m3c=2*m3m,m32=2*m3c)
      parameter (m3p=m3c+2)
      parameter (m1p=m1m+2)
      parameter (m1c=m1m)
      real   q1(m1,m2,m3),  q2(m1,m2,m3),  q3(m1,m2,m3)
      real   ru1(m1,m2,m3), ru2(m1,m2,m3), ru3(m1,m2,m3)
      real   dq(m1,m2,m3)
      real   pr(m1,m2,m3)
      real   dens(m1,m2,m3), ruro(m1,m2,m3), hro(m1,m2,m3)
      common / qarr1/ q1
      common / qarr2/ q2
      common / qarr3/ q3
      common /ruarr1/ ru1
      common /ruarr2/ ru2
      common /ruarr3/ ru3
      common /vaarr1/ dq
      common /vaarr2/ pr
      common /vaarr3/ dens
      common /vaarr4/ ruro
      common /vaarr5/ hro
      real       dph(m1+1,m2,m3)
      real       qcap(m1,m2,m3)
      real       rhs(m1*m2*m3)
      real       time
      real       epsil,re
      real     xr(m1m+2,m2m),work(m1,m2m)
      real     trigx1(3*m1m/2+1),trigxv(3*m1m/2+1),ak1(m1+1),ap(m1+1)
      integer  ifx1(13),ifxv(13)
      real     amphj(m2),acphj(m2),apphj(m2)
      real     amphk(m3),acphk(m1+1,m3),apphk(m3),acphkk(m3)
      real     yfis(m3,m2),w(m2*m3)
      real     dx1,dx1q,dx2,dx2q,dx3,dx3q
      real     rc(m2),rm(m2),zz(m3),zm(m3)
      real     g2rm(m2),g2rc(m2)
      real     g3rc(m3),g3rm(m3),str3
c     real     qqmax,qqtot
      real     rmed31,istr3,etdp3,strb3
      real     uaxN(3,m2),uaxS(3,m2),qaxN(3,m1,m2),qaxS(3,m1,m2)
      real     mode(m3,m1),spec(m3,m1)
      integer  nspec,ispec
      real     forc_scale
      integer  iforc
      common /specf/  mode,spec,nspec,ispec
      common /uqax/   uaxN,uaxS,qaxN,qaxS
      common /avesta/ istat,iresta,timeav,timerms
      common /axsym/  iaxsy
      common /cft/    nx3fft
      common /cordvo/ rc,rm,zz,zm
      common /corrt/  thetac(m1),ragc(m2),thetam(m1),ragm(m2)
      common /ctrdph/ amph(m2),acph(m2),apph(m2)
      common /d1/     re
      common /conv/   ray,pra,tfin,eps
      common /d13/    alx3
      common /d2/     nstop,ntst,npstf,ireset,ntime
      common /dd2/    tprint,tpin,tmax
      common /dim/    n1,n1m,n2,n2m,n3,n3m
      common /dimens/ alx3d
      common /diverg/ qcap
      common /dphc/   dph
      common /enerw/  nini,nfin,nstri
      common /fftcm1/ ifx1,ifx3(13)
      common /fftcm3/ trigx1,trigx3(m32),dism
      common /fftvel/ trigxv,ifxv
      common /force/  ifugo,icorio,ibuo,isca,ivisc
      common /favis/  fvisc1(m2,m3),fvisc2(m2,m3),fvisc3(m2,m3)
     %               ,fvisc(m2,m3)
      common /favise/  fvisc1e(m2,m3),fvisc2e(m2,m3),fvisc3e(m2,m3)
      common /idtopt/ cflmax,dtmax,cfllim,idtv
      common /indat/  h0,h1,delro
      common /indbo/  imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
     %               ,isym(m1),ksym(m3),ksyc(m3)
      common /inqca/  imxq,jmxq,kmxq
      common /inslip/ inslws,inslwn,inslwr,inslwri
      common /invoma/ iazm,jazm,kazm,irrm,jrrm,krrm
     %               ,irzm,jrzm,krzm
      common /inwal/  jpc(m2),jup(m2),jmc(m2),jum(m2)
      common /inwalk/ kpc(m3),kup(m3),kmc(m3),kum(m3)
      common /maxdi/  qqmax,qqtot
      common /mesh/   dx1,dx1q,dx2,dx2q,dx3,dx3q
      common /mmdens/ denmax,denmin,densm
      common /njump/  n1p,n2p,n3p
      common /parcoo/ r0
      common /peckle/ pec
      common /posma/  radm(m1),zetm(m1)
      common /posmax/ imaxv(ndv),jmaxv(ndv),kmaxv(ndv)
      common /rando/  rnd(10000)
      common /residue/resid
      common /rhsc/   rhs
      common /richar/ ri
      common /rossby/ ros
      common /smooth/ nlev,nn
      common /sonda/  coson(400,3),mason(400,3),nson
      common /stfft/  coecos(m3),coesin(m3)
      common /trici/  ami(m1,m2),aci(m1,m2),api(m1,m2),
     %                fei(m2,m1),q(m2,m1),s(m2,m1)
      common /tricj/  amj(m2),acj(m2),apj(m2)
      common /trick/  amk(m3),ack(m3),apk(m3)
      common /tscoe/  ga,ro,al,nsst
      common /tscoev/ gam(3),rom(3),alm(3)
      common /tstep/  dt,beta,ren
      common /velmax/ vmax(ndv),vmaxo(ndv)
      common /vopa/   vmx,sig,yc2mo,yc3mo
      common /vperin/ vper,rper,lamb,nsect
      common /perin/  epsil
      common /wrre/   nwrit,nread
      common /waves/  an(m3),ap,ak3(m3),ak1
C      NUOVE AGGIUNTE
      common /metrr/  g2rm,g2rc
      common /nonunif/strr,rext,rint,rmed,etdp,strb,istr
      common /corrk/  aprk(m2,m3),acrk(m2,m3),amrk(m2,m3)
      common /cor1j/  ap1j(m2),ac1j(m2),am1j(m2)
      common /cor1k/  ap1k(m2,m3),ac1k(m2,m3),am1k(m2,m3)
      common /cor2j/  ap2j(m2),ac2j(m2),am2j(m2)
      common /cor2je/ ap2je(m2),ac2je(m2),am2je(m2)
      common /cor2k/  ap2k(m2,m3),ac2k(m2,m3),am2k(m2,m3)
      common /cor3j/  ap3j(m2),ac3j(m2),am3j(m2)
      common /cor3k/  ap3k(m2,m3),ac3k(m2,m3),am3k(m2,m3)
      common /corscj/ apscj(m2),acscj(m2),amscj(m2)
      common /bcdens/ denbs(m1,m3),denbn(m1,m3)
      common /anusse/ anusslow,anussupp
      common/quafi/yfis,w
      common/cofisn/mw,np,mp
      common/phcok/amphk,acphk,apphk
      common/phcoj/amphj,acphj,apphj
      common /cor3ck/  ap3ck(m3),ac3ck(m3),am3ck(m3)
      common /cor3sk/  ap3sk(m3),ac3sk(m3),am3sk(m3)
      common /cor3ss/  ap3ssk(m3),ac3ssk(m3),am3ssk(m3)
      common /ametrz/  g3rc,g3rm,str3
      common/nonunz/   rmed31,istr3,etdp3,strb3
      common/factors/  udx3m(m3),udx2rm(m2),udx2rmq(m2),udx1vis1(m2,m3),
     %                 udx2c(m2),udx1vis3(m2,m3),
     %                 udx1vis2(m2),udx3c(m3),usrm(m2),usrmq(m2),
     %                 usrc(m2),
     %                 udx1qc(m2),udx1qm(m2),udx2m(m2),
     %                 rmsg2rm(m2),usg2rc(m2)
      common/sphe1/    usinpm(m3), usinpc(m3), sinpm(m3), sinpc(m3),
     %                 cospc(m3),cospm(m3),utanpc(m3),
     %                 udx3vis2(m3),utanpm(m3)
      common /forcing/ forc_scale,iforc
      common /icond/ indexForInitial
