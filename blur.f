      SUBROUTINE blur(earb, photarb, maxb, ear, photar, ne,
     &           vmax, enh, tmp)

      IMPLICIT NONE
c General routine for blurring by a line given by earb, photarb.

      INTEGER maxb, ne

      REAL earb(0:maxb), photarb(maxb)
      REAL ear(0:ne), photar(ne), photer(ne), tmp(ne)
      REAL enh, vmax

c earb and photarb  are energies and fluxes for line
c ear and photar are input energies and fluxes to be blurred

      REAL energy, sum
      REAL r1, r2, f1, f2
      REAL vrel, gamma, transmiss
      REAL Emin, Emax

      INTEGER halfmaxb
      INTEGER ie, je, i
      INTEGER if, if1, if2
      INTEGER ilow, ihigh, jlow, jhigh

      halfmaxb = maxb/2

      sum=0.
      do i=1,maxb
         sum=sum+photarb(i)
      enddo
      IF ( sum .NE. 0.0 .AND. sum .NE. 1.0 ) THEN
         do i=1,maxb
            photarb(i)=photarb(i)/sum
         enddo
      ENDIF

      vrel=vmax/300000.
      gamma=(1/(1-vrel**2))**0.5

c trap out the pathological case of sum = 0. In this case don't alter the
c input spectrum

      IF ( sum .EQ. 0.0 ) RETURN


c Initialize the convolution output

      do ie=1,ne
         tmp(ie)=0.0
      enddo

c Loop over energy ranges for the current dataset

      jlow = 1
      jhigh = 1

      do ie=1,ne

         energy=0.5*(ear(ie-1)+ear(ie))

c find first and last bins required in the convolution
         Emin=energy/(1.+vrel)/gamma
         Emax=energy/(1.-vrel)/gamma

         DO WHILE ( ear(jlow) .LT. Emin .AND. jlow .LT. ne )
            jlow = jlow + 1
         ENDDO
         IF ( jlow .GT. 1 ) jlow = jlow - 1
         DO WHILE ( ear(jhigh) .LT. Emax .AND. jhigh .LT. ne )
            jhigh = jhigh + 1
         ENDDO

         DO je = jlow, jhigh
            f1=ear(je-1)/energy*float(halfmaxb)*
     &           (1.0+1.0/float(halfmaxb))
            f2=ear(je)/energy*float(halfmaxb)*
     &           (1.0+1.0/float(halfmaxb))
            if1=INT(f1)
            r1=f1-if1
            if2=INT(f2)
            r2=f2-if2
            sum=0.0

            DO if=if1+1,if2-1
               IF(if.GE.1.AND.if.LE.maxb)THEN
                  sum=sum+photarb(if)
               ENDIF
            ENDDO
             
            IF(if1.GE.1.AND.if1.LE.maxb)THEN
               sum=sum+(1.0-r1)*photarb(if1)
            ENDIF
             
            IF(if2.GE.1.AND.if2.LE.maxb)THEN
               sum=sum+(r2)*photarb(if2)
            ENDIF
             
            IF(if1 .EQ. if2)THEN
               IF(if1.GE.1.AND.if1.LE.maxb)THEN
                  sum=sum-photarb(if1)
               ENDIF
            ENDIF

            IF (ear(je) .LE. Emin) THEN
               jlow=je
            ENDIF
            IF (ear(je+1) .GT. Emax) THEN
               jhigh=je+2
            ENDIF
            IF (jhigh .GT. ne-1) THEN
               jhigh=ne-1
            ENDIF


c       determine amount of absorption the backside will incur                    
c       enh refers to absorption value across full 2*r distance                    
c       since we presume the absorption stems from full
c       line of sight distance interior to the shell
            IF (ear(je) .LT. ear(ie)) THEN
              transmiss=EXP(-1.*(2.*ABS((ear(ie)-ear(je))/ear(je)/vrel)*
     &               enh*2.4*ear(je)**(-2.6)))
            ELSE
              transmiss = 1.
            ENDIF
            tmp(je)=tmp(je)+photar(ie)*sum*transmiss
         enddo
      enddo

      DO ie = 1, ne
         photar(ie) = tmp(ie)
         photer(ie) = photer(ie)
      ENDDO

      RETURN
      END

