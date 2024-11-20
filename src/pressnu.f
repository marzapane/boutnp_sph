c
************************************************************************
c  this subroutine perform the calculation of trigz for temperton fft
c
      subroutine fftqua
      include'param.f'
      pi=2.*asin(1.)
      n1mh=n1m/2+1
      n1mp=n1mh+1
c
c     wave number definition
c
      do 18 i=1,n1mh
        ap(i)=(i-1)*2.*pi
        ap(n1mh+i)=(i-1)*2.*pi
   18 continue
c     do 19 i=n1mp,n1+1
c       ap(i)=-(n1m-i+1)*2.*pi
c  19 continue
      nx1fft=n1m
      if(n1m.ne.1) call fftfax(nx1fft,ifx1,trigx1)
c
c   modified wave number
c
c     do 28 i=1,n1m
      do 28 i=1,n1+1
        ak1(i)=2.*(1.-cos(ap(i)/n1m))*dx1q
   28 continue
c
      return
      end
************************************************************************
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 and x1to use the real fourier transform
c
      subroutine phcalc
      include'param.f'
      dimension yfisl(m3,m2),wl1(m3*m2),wl2(m3*m2)
      dimension amjl(m2),acjl(m2),apjl(m2)
      dimension amkl(m3),ackl(m3),apkl(m3)
      my=m3
c
      pi=2.*asin(1.)
c
c
c   fft applied to the x1 direction to the
c   complex coeff. from cos fft
c   from physical to wave number space
c
      n1mh=n1m/2+1
      n1md=n1m+2
!$OMP  PARALLEL DO
!$OMP$ SHARED(dph,n1m,n2m,n3m,n1mh,n1md,trigx1,ifx1)
!$OMP$ PRIVATE(xr,work,is,ip,id,ii,i,j,k)
      do k=1,n3m
        do j=1,n2m
          xr(1,j)=dph(n1m,j,k)
          xr(n1md,j)=dph(1,j,k)
          do i=1,n1m
            is=i+1
            xr(is,j)=dph(i,j,k)
          end do
        end do
        call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n2m,-1)
        do j=1,n2m
          do i=1,n1mh
            ip=2*i
            id=2*i-1
            ii = i + n1mh
            dph(ii,j,k)=xr(id,j)
            dph(i,j,k)=xr(ip,j)
          end do
        end do
      end do
!$OMP  END PARALLEL DO
      go to 9874
c
c  North pole
      k=1
      do j=1,n2m
c       dph(1,j,k) = 0.
c       dph(n1mh,j,k) = 0.
        s1 = 0.
c       s2 = 0.
        do i=1,n1m/2
          s1 = s1 + float(i-1)*dph(i,j,k)
c         s2 = s2 + float(i-1)*dph(n1mh+i,j,k)
        end do
        dph(n1mh,j,k) = -s1/float(n1m/2)
c       dph(n1+1,j,k) = -s2/float(n1m/2)
      end do
c
c  South pole
      k=n3m
      do j=1,n2m
c       dph(1,j,k) = 0.
c       dph(n1mh,j,k) = 0.
        s1 = 0.
c       s2 = 0.
        do i=1,n1m/2
          s1 = s1 + float(i-1)*dph(i,j,k)
c         s2 = s2 + float(i-1)*dph(n1mh+i,j,k)
        end do
        dph(n1mh,j,k) = -s1/float(n1m/2)
c       dph(n1+1,j,k) = -s2/float(n1m/2)
      end do
c
 9874 continue
c    
c     inversion of the matrix in the x2 and x3 directions (real part)
c    
c
c     dph_max = 0.
!$OMP  PARALLEL DO
!$OMP$ SHARED(acphk,amphk,apphk,amphj,acphj,apphj)
!$OMP$ SHARED(dph,n1,n1mh,n2m,n3m,w)
!$OMP$ PRIVATE(yfisl,amkl,ackl,apkl,ierror)
!$OMP$ PRIVATE(wl1,ii,i,j,k,kl,jl,ll)
CC!$OMP$ PRIVATE(wl2)
CC    do i=1,n1mh
      do i=1,n1+1
CC      ii = i + n1mh
        do kl=1,m3
          do jl=1,m2
            ll = kl + (jl-1)*m3
            wl1(ll)=w(ll)
          end do
        end do
        do k=1,n3m
          amkl(k)=amphk(k)
          apkl(k)=apphk(k)
          ackl(k)=acphk(i,k)
          do j=1,n2m
CC          yfisl(k,j)=dph(ii,j,k)
            yfisl(k,j)=dph(i,j,k)
          end do
        end do
        call blktri(1,1,n2m,amphj,acphj,apphj
     %             ,1,n3m,amkl,ackl,apkl,
     %              m3,yfisl,ierror,wl1)
c       write(6,*) ' 2 IERROR = ',ierror
        do j=1,n2m
          do k=1,n3m
CC          dph(ii,j,k)=yfisl(k,j)
            dph(i,j,k)=yfisl(k,j)
c           dph_max = max(abs(dph(i,j,k)),dph_max)
          end do
        end do
      end do
!$OMP END PARALLEL DO
c     write(6,*) 'dph_max=',dph_max
c
c    
c   inverse fft applied to the phi x1 direction
c   from wave number space to physical space
c    
      n1mu=n1m-1
!$OMP  PARALLEL DO
!$OMP$ SHARED(dph,n1m,n1mu,n3m,n2m,n1mh,trigx1,ifx1)
!$OMP$ PRIVATE(xr,work,is,ip,id,i,j,k,ii)
      do k=1,n3m
        do j=1,n2m
          do i=1,n1mh
            ii = i + n1mh
            ip=2*i
            id=2*i-1
            xr(id,j)=dph(ii,j,k)
            xr(ip,j)=dph(i,j,k)
          enddo
        enddo
        call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n2m,+1)
        do j=1,n2m
          dph(n1m,j,k)=xr(1,j)
          do i=1,n1mu
            is=i+1
            dph(i,j,k)=xr(is,j)
          enddo
        enddo
      enddo
!$OMP  END PARALLEL DO
c
      pmed = dph(1,1,1)
c     pmed = 0.

!$OMP  PARALLEL DO
!$OMP$ SHARED(dph,pmed,n1m,n3m,n2m)
!$OMP$ PRIVATE(i,j,k)
      do i=1,n1m
        do j=1,n2m
          do k=1,n3m
            dph(i,j,k) = dph(i,j,k) - pmed
          end do
        end do
      end do
!$OMP  END PARALLEL DO
      return
      end
************************************************************************
c   in this subr the coefficients of the poisson eq. for dph
c   are calculated this subr. is called only at the beginning
c
      subroutine phini
      include'param.f'
      call fftqua
      n1mh=n1m/2+1
c
c   tridiagonal matrix coefficients at each k and i
c   x1 and x3 cartesian coordinates
c
      do jc=1,n2m
        jp=jc+1
        a22icc=rc(jc)**2*jmc(jc)*dx2q/g2rc(jc)
        a22icp=rc(jp)**2*jpc(jc)*dx2q/g2rc(jp)
        ugmmm=1./g2rm(jc)
        amphj(jc)=a22icc*ugmmm
        apphj(jc)=a22icp*ugmmm
        acphj(jc)=-(amphj(jc)+apphj(jc))
      end do
c
      do kc=1,n3m
        km=kmv(kc)
        kp=kc+1
        a33icc=kmc(kc)*dx3q*sin(zz(kc))/g3rc(kc)
        a33icp=kpc(kc)*dx3q*sin(zz(kp))/g3rc(kp)
        ugmmm=1./(g3rm(kc)*sin(zm(kc)))
        amphk(kc)=a33icc*ugmmm
        apphk(kc)=a33icp*ugmmm
        acphkk(kc)=-(amphk(kc)+apphk(kc))
      enddo
c
c     do ic=1,n1mh
      do ic=1,n1+1
        do kc=1,n3m
          acphk(ic,kc)=acphkk(kc)-ak1(ic)/(sin(zm(kc)))**2 
        enddo
      enddo
c     do ic=1,n1mh
c       ii = ic + n1mh
c       do kc=1,n3m
c         acphk(ii,kc)=acphk(ic,kc)
c       enddo
c     enddo
      an3=n2m
      ax=log(an3)/log(2.)
      k=ax+1
      L=2**(K+1)
      mw=(K-2)*L+K+5+MAX(2*n2m,6*n3m)
      if(mw.gt.(m3*m2)) then
        write(6,*) ' WARNING THE ARRAY W HAS INSUFFICIENT '
        write(6,*) ' DIMENSIONS : mw = ',mw
        pause
      end if
     
      np=1
      mp=1
      do j=1,n2m
        do k=1,n3m
          yfis(k,j)=1.
        enddo
      enddo
      my=m3
        call blktri(0,np,n2m,amphj,acphj,apphj
     %             ,mp,n3m,amphk,acphkk,apphk,
     %              my,yfis,ierror,w)
      write(6,*) ' 1 IERROR = ',ierror 
      return
      end
