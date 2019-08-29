      
      SUBROUTINE shellsb(ear, ne, param, nenrgy, nlow, nvel,
     &           efil, ebin, low, vel, shelln, bin, start, end, 
     &           fstart, fend, photar)


      INTEGER ne, nenrgy, nlow, nvel

      REAL param(5), ear(0:ne), photar(ne)
      REAL efil(0:*), ebin(*), low(*), vel(*), shelln(*), bin(*)
      REAL start(*), end(*), fstart(*), fend(*)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  model to calculate the line shape for an expanding 
c  (partial) shell.
c
c  parameters :
c       1        "energy"
c       2        opening angle lower (0-90 deg)
c       3        opening angle upper (0-90 deg)
c       4        line-of-sight angle (0-90 deg)
c       5        maximum velocity (km/s)
c
c Franz Bauer and 24601
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL smo(3),x,y,tmp1,spm

      INTEGER i,ismo,j

      REAL aopenl,aopenu,alos,vmax,b,dera
      INTEGER ilowi,ilowo,ivel,ilos,theta,fila
      REAL whgh_los,wlow_los,whgh_vel,wlow_vel
      REAL angulo
 
      DATA dera /0.017453293/
      aopenl=param(2)
      aopenu=param(3) 
      IF (aopenl .GT. aopenu) THEN
        tmp1=aopenu
        aopenu=aopenl
        aopenl=tmp1
      ENDIF
      alos=param(4)
      vmax=param(5)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      efil(0) = 0.
      efil(1) = (3*ebin(1)-ebin(2)) * param(1) / 2.
      DO i = 2, nenrgy
         efil(i) = (ebin(i-1)+ebin(i)) * param(1) / 2.
      ENDDO
      efil(nenrgy+1)=(3*ebin(nenrgy)-ebin(nenrgy-1))*param(1)/2.
      efil(nenrgy+2) = 1.e6
      DO i = 1, nenrgy+2
         bin(i) = 0.
      ENDDO
      DO i=1,nlow
        print *,low(i)
      ENDDO

c Locate theta angles

      ilowi=1
      ilowo=nlow
      DO i=1,nlow-1
        IF(aopenl.GE.low(i).AND.aopenl.LE.low(i+1)) ilowi=i
        IF(aopenu.GE.low(i).AND.aopenu.LE.low(i+1)) ilowo=i+1
        IF(alos.GE.low(i).AND.alos.LE.low(i+1)) ilos=i
      ENDDO   
c weight due to line-of-sigh angle
      whgh_los = (low(ilos+1)-alos)/(low(ilos+1)-low(ilos))
      wlow_los = (alos-low(ilos))/(low(ilos+1)-low(ilos))

c locate velocities
      DO i=1,nvel-1
        IF(vmax.GE.vel(i).AND.vmax.LE.vel(i+1)) ivel=i
      ENDDO

c weight associated to the maximum velocity
      whgh_vel = (vel(ivel+1)-vmax)/(vel(ivel+1)-vel(ivel))
      wlow_vel = (vmax-vel(ivel))/(vel(ivel+1)-vel(ivel))

      b=1.-cos((aopenu-aopenl)*dera)
c      open (unit = 93, file = "photar.txt")
c      print *,low(ilowi),low(ilowo),low(ilos),vel(ivel)
      DO theta=ilowi+1,ilowo-1
c        print *,theta,low(theta)        
        IF(theta.EQ.ilowi+1) THEN
          angulo=abs(b-(1.-cos((aopenu-low(theta))*dera)))
        ELSE
          IF(theta.EQ.ilowo-1)THEN
            angulo=1.-cos((aopenu-low(theta))*dera)
          ELSE     
            x= cos((aopenu-((low(theta+1)+low(theta))/2.))*dera)
            y= cos((aopenu-((low(theta)+low(theta-1))/2.))*dera)
            angulo=abs(x-y)
          ENDIF
        ENDIF

        fila=INT((((theta-1)*46*60)+((ilos-1)*60)+ivel)*nenrgy)
c        print *,wlow_los,whgh_los,wlow_vel,whgh_vel
c        print *,fila,low(theta-1),low(ilos-1),vel(ivel)
        DO i=1,nenrgy
            bin(i+1)=bin(i+1)+wlow_los*shelln(fila+i)*wlow_vel*angulo    
        ENDDO
        fila=fila+nenrgy
c        print *,fila,low(theta-1),low(ilos-1),vel(ivel+1)
        DO i=1,nenrgy
            bin(i+1)=bin(i+1)+wlow_los*shelln(fila+i)*whgh_vel*angulo
        ENDDO

        fila=INT((((theta-1)*46*60)+((ilos)*60)+ivel)*nenrgy)
c        print *,fila,low(theta-1),low(ilos),vel(ivel)
        DO i=1,nenrgy
            bin(i+1)=bin(i+1)+whgh_los*shelln(fila+i)*wlow_vel*angulo
        ENDDO
        fila=fila+nenrgy
c        print *,fila,low(theta-1),low(ilos),vel(ivel+1)
        DO i=1,nenrgy
            bin(i+1)=bin(i+1)+whgh_los*shelln(fila+i)*whgh_vel*angulo
        ENDDO
      ENDDO

c  Perform 3-pt smoothing to reduce the small wiggles which result
c  from the finite number of radial bins

      DO i = 1, 2
         smo(i) = bin(i+1)
      ENDDO
      ismo = 3
      DO i = 2, nenrgy-1
         smo(ismo) = bin(i+2)
         ismo = ismo + 1
         IF(ismo.EQ.4) ismo = 1
         bin(i+1) = 0.
         DO j = 1, 3
            bin(i+1) = bin(i+1) + smo(j)
         ENDDO
         bin(i+1) = bin(i+1)/3.
      ENDDO

c  Rebin onto passed energies

      CALL inibin((nenrgy+2), efil, ne, ear, start, end,
     &            fstart, fend, 0)
      CALL erebin((nenrgy+2), bin, ne, start, end, fstart, fend, photar)

c normalise values to total

      spm = 0.
      DO i = 1, ne
         spm = spm + photar(i)
      ENDDO
      IF ( (spm .NE. 0.) .AND. (spm .NE. 1.) ) THEN
         DO i = 1, ne
            photar(i) = photar(i)/spm
         ENDDO
      ENDIF

      RETURN
      END

