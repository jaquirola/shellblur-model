      SUBROUTINE doshblur (ear, ne, param, ifl, photar, tmp)

      IMPLICIT NONE

      INTEGER ne, ifl
      REAL ear(0:ne), param(5), photar(ne), tmp(ne)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Subroutine to smooth a model spectrum by an expanding shell geometry
c FEB March 2012

c Arguments :
c    ear      r        i: energy range bins
c    ne       i        i: total number of elements in photar array
c    param    r        i: model parameters (defined below)
c    ifl      i        i: data set
c    photar   r      i/r: model flux

c internal variables :
c    tmp      r        i: Dynamically allocated workspace array

c input parameters :
c       1        opening angle lower (0-90 deg)
c       2        opening angle upper (0-90 deg)
c       3        line-of-sight angle (0-90 deg)
c       4        maximum velocity (km/s)
c       5        maximum uniform line-of-sight ejecta nh (1e22)
c Franz Bauer and 24601
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      INTEGER MAXB, HALFMAXB
      PARAMETER (MAXB=1000,HALFMAXB=MAXB/2)

      REAL energy,sum,earc(0:MAXB)
      REAL parama(5), photara(MAXB), photera(MAXB)
      REAL paramb(5), photarb(MAXB), photerb(MAXB)
      REAL photarc(MAXB)
      REAL enh, vmax

      INTEGER ie, i, pind

      LOGICAL first

      save first,pind,earc,parama,paramb,photara,photarb
      DATA first/.true./
      DATA pind/0/
      IF (first) THEN
         do i=0,MAXB
           earc(i)=20.0*float(i)/float(MAXB)
         end do
         first=.false.
         pind=1
      ENDIF
      parama(1)=10.
      parama(2)=param(1)
      parama(3)=param(2)
      parama(4)=param(3)
      parama(5)=param(4)
      call xsshell(earc,MAXB,parama,ifl,photara,photera)
      sum=0.
      do i=1,maxb
        sum=sum+photara(i)
      enddo
      IF ( sum .NE. 0.0 .AND. sum .NE. 1.0 ) THEN
        do i=1,maxb
          photara(i)=photara(i)/sum
        enddo
      END IF
      photarc=photara
      vmax=param(4)
      enh=param(5)
      CALL blur(earc, photarc, MAXB, ear, photar, ne, 
     &         vmax, enh, tmp)

      RETURN
      END

