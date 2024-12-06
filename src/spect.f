c***********************************************************************
c                                                                      *
c            ************** SPECT  ************************            *
c                                                                      *
c***********************************************************************
      subroutine spect(time)
      include 'param.f'
      dimension spec(m3,m1),rwc(m2),rwm(m2)
      real      rtr,rsc,rsm,mode(m3,m1)
      real   xxr(m1m+2,m3m),wwork(m1m+1,m3m)
      integer pos
      character*30 filsp,buffer
      character*14 filth,filte,filba,filen,filso,filve,filnu,filbaper
      common/names/filth,filte,filba,filen,filso,filve,filnu,filbaper
c
c***********************************************************************
c
      call fftv
      dvol = 1./(dx2 * dx3)
      pi = 2. * asin (1.)
      n1mh=n1m/2+1
      n1md=n1m+2
      rsc = 0.0d0
      rsm = 0.0d0
      do j = 1,n2
        rwc(j) = rc(j)*rc(j) * g2rc(j)/dx2
        rsc = rsc + rwc(j)
        rwm(j) = rm(j)*rm(j) * g2rm(j)/dx2
        rsm = rsm + rwm(j)
      end do
c initialize xxr
      do k = 1, n3m
        do i = 1, n1md
          xxr(i,k) = 0.0d0
        end do
      end do
c average v_r in xxr
      do j = 1,n2
        rtr = rwc(j)/rc(j)
        do k = 1,n3m
          xxr(1,k) = xxr(1,k) + q2(n1m,j,k)*rtr
          do i = 1,n1m
            is = i+1
            xxr(is,k) = xxr(is,k) + q2(i,j,k)*rtr
          end do
          xxr(n1md,k) = xxr(n1md,k) + q2(1,j,k)*rtr
        end do
      end do
      do k = 1, n3m
        do i = 1, n1md
          xxr(i,k) = xxr(i,k)/rsc
        end do
      end do
c
c     fft applied to the (q_2) by fft99
c   from physical to wave number space
c
      call fft99(xxr,wwork,trigxv,ifxv,1,m1+1,n1m,n3m,-1)
      do 31 i=1,n1mh
        ip=2*i
        id=2*i-1
        do 41 k=1,n3m
          spec(k,i) = real( cmplx(xxr(id,k),xxr(ip,k))
     1               *conjg(cmplx(xxr(id,k),xxr(ip,k))) )
  41    continue
  31  continue
c
c
      do k = 1, n3m
        do i = 1, n1md
          xxr(i,k) = 0.
        end do
      end do
      do 2 j=1,n2m
        do 12 k=1,n3m
          xxr(1,k) = xxr(1,k) + q1(n1m,j,k)*rwm(j)
          do 22 i=1,n1m
            is=i+1
            xxr(is,k) = xxr(is,k) + q1(i,j,k)*rwm(j)
  22      continue
          xxr(n1md,k) = xxr(nimd,k) + q1(1,j,k)*rwm(j)
  12    continue
  2   continue
      do k = 1, n3m
        do i = 1, n1md
          xxr(i,k) = xxr(i,k)/rsm
        end do
      end do
c
c     fft applied to the (q_1) by fft99
c   from physical to wave number space
c
      call fft99(xxr,wwork,trigxv,ifxv,1,m1+1,n1m,n3m,-1)
      do 32 i=1,n1mh
        ip=2*i
        id=2*i-1
        do 42 k=1,n3m
          spec(k,i) = spec(k,i)
     1              + real( cmplx(xxr(id,k),xxr(ip,k))
     1              * conjg(cmplx(xxr(id,k),xxr(ip,k))) )
  42      continue
  32  continue
c
c
      do k = 1, n3m
        do i = 1, n1md
          xxr(i,k) = 0.
        end do
      end do
      do 3 j=1,n2m
        do 13 k=1,n3m
          xxr(1,k) = xxr(1,k) + q3(n1m,j,k)*rwm(j)
          do 23 i=1,n1m
            is=i+1
            xxr(is,k) = xxr(is,k) + q3(i,j,k)*rwm(j)
   23     continue
          xxr(n1md,k) = xxr(n1md,k) + q3(1,j,k)*rwm(j)
   13   continue
   3  continue
      do k = 1, n3m
        do i = 1, n1md
          xxr(i,k) = xxr(i,k)/rsm
        end do
      end do
c
c     fft applied to the (q_3) by fft99
c   from physical to wave number space
c
      call fft99(xxr,wwork,trigxv,ifxv,1,m1+1,n1m,n3m,-1)
      do 33 i=1,n1mh
        ip=2*i
        id=2*i-1
        do 43 k=1,n3m
          spec(k,i) = spec(k,i)
     1              + real( cmplx(xxr(id,k),xxr(ip,k))
     1              * conjg(cmplx(xxr(id,k),xxr(ip,k))) )
  43    continue
  33  continue
c
c***********************************************************************
      do i=1,n1mh
        do k=1,n3m
          mode(k,i) = real(i-1)/sin(zm(k))
        end do
      end do
c***********************************************************************
c
      write(filsp,'(A6,I8.8,A4)') 'spctr_', ntime, '.bin'
      open(39,file=filsp,form='unformatted',status='unknown')
      write(39) n3m, n1mh
      write(39) (((mode(k,i)),k=1,n3m),i=1,n1mh)
      write(39) (((spec(k,i)),k=1,n3m),i=1,n1mh)
      write(39) time
      close(39)
      return
      end
