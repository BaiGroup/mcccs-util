      implicit none

      integer i,j,k,ibin
      integer npoints,nbins
      double precision a,b,c
      double precision x(500),y(500),dx(500),dy(500)
      double precision wx(500),wy(500),w(500)
      double precision mupper,mlower,mtrial,merror
      double precision mbest,cbest,mbesterror
      double precision xmean,ymean,u(500),v(500)
      double precision xsum,ysum,wsum
      double precision aa,bb,cc,dd,ee,ff,gg,hh,jj
      double precision dmdx(500),dmdy(500)
      double precision dcdx(500),dcdy(500)
      double precision bigs,unweighteds
      double precision mvariance,cvariance

      mbesterror = 1d20

      write(6,*)'How many data points are there?'
      read(5,*)npoints

      do i = 1,npoints
         read(10,*)x(i),dx(i),y(i),dy(i)
         wx(i) = 1d0/dx(i)**2d0
         wy(i) = 1d0/dy(i)**2d0
      enddo

      write(6,*)'input slope lower and upper bounds'
      read(5,*)mupper,mlower
      write(6,*)'input number of grid points'
      read(5,*)nbins
      
      do ibin = 1,nbins
         mtrial = mlower + ibin*((mupper-mlower)/nbins)
         xsum = 0d0
         ysum = 0d0
         wsum = 0d0
         do i = 1,npoints
            w(i) = (wx(i)*wy(i)) / (wy(i)*mtrial**2d0 + wx(i))
            xsum = xsum + w(i)*x(i)
            ysum = ysum + w(i)*y(i)
            wsum = wsum + w(i)
         enddo
         xmean = xsum / wsum
         ymean = ysum / wsum
         do i = 1,npoints
            u(i) = x(i) - xmean
            v(i) = y(i) - ymean
         enddo
         a = 0d0
         b = 0d0
         c = 0d0
         do i = 1,npoints
            a = a + (u(i)*v(i)*w(i)**2d0) / wx(i)
            b = b + w(i)**2d0 * ((u(i)**2d0/wy(i)) - 
     *                           (v(i)**2d0/wx(i)))
            c = c + (u(i)*v(i)*w(i)**2d0) / wy(i)
         enddo
         a = a * mtrial**2d0
         b = b * mtrial
         merror = a + b - c
         if (dabs(merror) .lt. mbesterror) then
            mbest = mtrial
            cbest = ymean - mtrial*xmean
            mbesterror = dabs(merror)
         endif
      enddo

c      mbest = 2.7d0
c      cbest = -0.63d0

      mvariance = 0d0
      cvariance = 0d0
      bigs = 0d0
      unweighteds = 0d0

      xsum = 0d0
      ysum = 0d0
      wsum = 0d0
      do i = 1,npoints
         w(i) = (wx(i)*wy(i)) / (wy(i)*mbest**2d0 + wx(i))
         xsum = xsum + w(i)*x(i)
         ysum = ysum + w(i)*y(i)
         wsum = wsum + w(i)
      enddo
      xmean = xsum / wsum
      ymean = ysum / wsum
      do i = 1,npoints
         u(i) = x(i) - xmean
         v(i) = y(i) - ymean
      enddo
      a = 0d0
      b = 0d0
      c = 0d0
      do i = 1,npoints
         a = a + (u(i)*v(i)*w(i)**2d0) / wx(i)
         b = b + w(i)**2d0 * ((u(i)**2d0/wy(i)) -
     &                        (v(i)**2d0/wx(i)))
         c = c + (u(i)*v(i)*w(i)**2d0) / wy(i)
      enddo

      hh = 0d0
      jj = 0d0
      do i = 1,npoints
         hh = hh + (v(i)*w(i)**2d0)/wx(i)
         jj = jj + (u(i)*w(i)**2d0)/wx(i)
      enddo
      hh = hh * (-2d0*mbest/wsum)
      jj = jj * (-2d0*mbest/wsum)  

      do j = 1,npoints
         aa = 0d0
         bb = 0d0
         cc = 0d0
         dd = 0d0
         ee = 0d0
         ff = 0d0
         gg = 0d0
         do i = 1,npoints
            aa = aa + ( (u(i)*v(i)*w(i)**3d0)/(wx(i)**2d0) ) 
            bb = bb + w(i)**2d0 * ((4d0*mbest * ((w(i)/wx(i)) *
     &                ( (u(i)**2d0/wy(i)) - (v(i)**2d0/wx(i))))) - 
     &                (2d0*v(i)*hh / wx(i)) + (2d0*u(i)*jj / wy(i)))
            cc = cc + (w(i)**2d0 / wy(i)) * 
     &                (((4*mbest*w(i)*u(i)*v(i)) / wx(i) ) +
     &                 v(i)*jj + u(i)*hh)
            if (i .eq. j) then
               dd = dd + ((v(i)*w(i)**2d0)/wx(i)) *
     &                   (1d0 - (w(j)/wsum))
               ee = ee + ((u(i)*w(i)**2d0)/wy(i)) *
     &                   (1d0 - (w(j)/wsum))
               ff = ff + ((v(i)*w(i)**2d0)/wy(i)) *
     &                   (1d0 - (w(j)/wsum))
               gg = gg + ((u(i)*w(i)**2d0)/wx(i)) *
     &                   (1d0 - (w(j)/wsum))
            else
               dd = dd + ((v(i)*w(i)**2d0)/wx(i)) *
     &                   (-1d0*(w(j)/wsum))
               ee = ee + ((u(i)*w(i)**2d0)/wy(i)) *
     &                   (-1d0*(w(j)/wsum))
               ff = ff + ((v(i)*w(i)**2d0)/wy(i)) *
     &                   (-1d0*(w(j)/wsum))
               gg = gg + ((u(i)*w(i)**2d0)/wx(i)) *
     &                   (-1d0*(w(j)/wsum))
            endif
         enddo

         aa = aa * 4d0 * mbest
         aa = aa - ( (wsum*hh*jj) / mbest )
         bb = bb * (-1d0)
         cc = cc * (-1d0)
         ee = ee * 2d0

         dmdx(j) = ((dd*mbest**2d0) + (ee*mbest) - ff) /
     &             ((2d0*mbest*a) + b - 
     &              (aa*mbest**2d0) + (bb*mbest) - cc)
         dmdx(j) = dmdx(j) * (-1d0)

         dmdy(j) = ((gg*mbest**2d0) - (2d0*dd*mbest) - (ee/2d0)) /
     &             ((2d0*mbest*a) + b - 
     &              (aa*mbest**2d0) + (bb*mbest) - cc)
         dmdy(j) = dmdy(j) * (-1d0)

         dcdx(j) = (hh - (mbest*jj) - xmean) * dmdx(j) -
     &             (mbest*w(j)/wsum)
         dcdy(j) = (hh - (mbest*jj) - xmean) * dmdy(j) +
     &             (w(j)/wsum)

         mvariance = mvariance + (1d0/wy(j))*dmdy(j)**2d0 + 
     &                           (1d0/wx(j))*dmdx(j)**2d0
         cvariance = cvariance + (1d0/wy(j))*dcdy(j)**2d0 +
     &                           (1d0/wx(j))*dcdx(j)**2d0

         bigs = bigs + w(j)*(y(j) - cbest - mbest*x(j))**2d0
         unweighteds = unweighteds + 
     &                 (y(j) - cbest - mbest*x(j))**2d0

      enddo

      mvariance = mvariance * (bigs/(npoints-2))
      cvariance = cvariance * (bigs/(npoints-2))

      write(6,*)mbesterror
      write(6,*)mbest,mvariance
      write(6,*)cbest,cvariance
      write(6,*)bigs,unweighteds

      end
