!------------------------------------------------------------
!------------------------------------------------------------
!                 PROGRAMA UPSCALING         
!                                       
!
! NAME: MARCIO RENTES BORGES
! DATA: 17/02/2012
! 
!
!------------------------------------------------------------
!------------------------------------------------------------
!   PROGRAMA UTILIZADO PARA A REALIZAR O UPSCALING DO
!    CAMPO DE PERMEABILIDADES EM MALHA GROSSA
!
!------------------------------------------------------------
!------------------------------------------------------------
!
MODULE VARIAVEIS
  CHARACTER(LEN=128):: FILEIN,FILEOUT,FILEINI
  CHARACTER(LEN=128):: FILE_OUTX,FILE_OUTY
  DOUBLE PRECISION  :: FDX,FDY,FDZ,KG,RHO
  DOUBLE PRECISION  :: VARIAN,FBETA,VAR,MED,U
  INTEGER           :: NX,NY,NXG,NYG
  INTEGER           :: NI,NF,NTYPE,NTIPO
  DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:) :: PERM,PERMT
  DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:) :: PERMG,PERMGT
END MODULE VARIAVEIS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM MAIN
!
  USE VARIAVEIS
!
  IMPLICIT NONE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LIST OF LOCAL VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER                    :: I,J,K,NS,NERROREAD
  CHARACTER*128              :: COMMAND
  REAL                       :: TIME
  REAL, DIMENSION(2)         :: TARRAY
  INTEGER                    :: CONT,ISTAT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ INPUT DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NERROREAD = 0
  CALL READ_INPUT(NERROREAD)
!
  IF(NERROREAD.NE.1)THEN
     WRITE(*,*)'#### LEITURA ERRADA ####'
     WRITE(*,*)'### DO ARQUIVO INPUT ###'
     STOP
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ALOCACAO DE MEMORIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(PERM(NX*NY))
  ALLOCATE(PERMT(NX*NY))
  ALLOCATE(PERMG(NXG*NYG))
  ALLOCATE(PERMGT(NXG*NYG))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LOOP PARA OS CAMPOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CONT = NF-NI
  I=0
  DO J=0,CONT
     IF(NTIPO.EQ.0)THEN
        CALL READPERM(I)
        CALL UPSCALE(0)
        CALL PRINT_OLD(FILE_OUTX,PERMGT,NXG,NYG,I)
        CALL TOMBA()
        CALL UPSCALE(0)
        CALL DTOMBA()
        CALL PRINT_OLD(FILE_OUTY,PERMGT,NXG,NYG,I)
     END IF
     IF(NTIPO.EQ.1)THEN
        CALL READPERM(I)
        CALL UPSCALE(1)
        CALL PRINT_OLD(FILE_OUTX,PERMGT,NXG,NYG,I)
     END IF
     I=I+1
  ENDDO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ INPUT DATA OF RANDOM FIELDS !!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DEALLOCATE(PERM,PERMT,PERMG,PERMGT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END PROGRAM MAIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE READ_INPUT(NERROREAD)
!
        USE VARIAVEIS
        IMPLICIT NONE
!
        INTEGER :: IN_FILE,ISTAT,I,J
        INTEGER :: CONT,OUT_FILE
        CHARACTER(LEN=128):: FILEF
        INTEGER, INTENT(OUT) :: NERROREAD
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FILEF = './entrada.in'
        WRITE(*,*)'###################################'
        WRITE(*,*)'## READ INPUT DATA FROM: ##########'
        WRITE(*,*)FILEF
        WRITE(*,*)'###################################'
!
        IN_FILE = 27
        OPEN(UNIT=IN_FILE,FILE=FILEF,STATUS='UNKNOWN',&
             FORM='FORMATTED',IOSTAT=ISTAT)
        IF(ISTAT.NE.0)THEN
           WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILEF
           STOP
        END IF
!
!! NUMERO DE CAMPOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        READ(IN_FILE,103)NTIPO
        WRITE(*,104)NTIPO
        IF(NTIPO.EQ.0) WRITE(*,*)'DURLOFSKY'
        IF(NTIPO.EQ.1) WRITE(*,*)'MEDIA'
!! NUMERO DE CAMPOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        READ(IN_FILE,100)NI,NF
        CONT = NF-NI+1
        WRITE(*,101)NI,NF,CONT
!
!! NOME DO ARQUIVO DE PERMEABILIDADES MALHA FINA !!!!!!!!
        READ(IN_FILE,200)FILEIN
        FILEIN=ADJUSTL(TRIM(FILEIN))
        WRITE(*,201)FILEIN
!
!! LEITURA DE KG E RHO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        READ(IN_FILE,*)KG,RHO
        WRITE(*,300)KG,RHO
        IF(RHO.LE.0.D0)THEN
           WRITE(*,*)'######## CAMPO HOMOGENEO #########'
           WRITE(*,*)'##### VERIFICAR VALOR DE RHO #####'
           STOP
        END IF
!
!! NOME DO ARQUIVO DE PERMEABILIDADES MALHA GROSSA !!!!!!
        READ(IN_FILE,200)FILEOUT
        FILEOUT=ADJUSTL(TRIM(FILEOUT))
        WRITE(*,202)FILEOUT
        FILE_OUTX=TRIM(FILEOUT)//TRIM('x_')
        FILE_OUTX=ADJUSTL(TRIM(FILE_OUTX))
        FILE_OUTY=TRIM(FILEOUT)//TRIM('y_')
        FILE_OUTY=ADJUSTL(TRIM(FILE_OUTY))
!
!! DADOS DA MALHA GROSSA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        READ(IN_FILE,100)NXG,NYG
        WRITE(*,102)NXG,NYG
!       
!! LEITURA DAS INFORMACOES DO CAMPOS DE ENTRADA !!!!!!!!!
        CALL READ_INFO()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        CLOSE(UNIT=IN_FILE)
!
        NERROREAD=1
!
103 FORMAT(I10)
104 FORMAT(/,'TIPO DE UPSCALING',I10)
100 FORMAT(2I10)
101 FORMAT(/,'NUMERO DA PRIMEIRA AMOSTRA:',I10,/,&
             'NUMERO FINAL DA AMOSTRA   :',I10,/,&
             'NUMERO TOTAL DE AMOSTRAS  :',I10)
102 FORMAT(/,'NUMERO DE ELEMENTOS EM (X):',I10,/,&
             'NUMERO DE ELEMENTOS EM (Y):',I10)
200 FORMAT(A)
201 FORMAT(/,'NOME DO ARQUIVO DE PERMEABILIDADES (MALHA FINA)  : ',A)
202 FORMAT('NOME DO ARQUIVO DE PERMEABILIDADES (MALHA GROSSA): ',A)
300 FORMAT('VALOR DE KG =',F12.5,/,&
           'VALOR DE RHO=',F12.5,/)
!
      END SUBROUTINE READ_INPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ROTINA PARA LEITURA DAS INFORMACOES DOS CAMPOS !!!!!!!
      SUBROUTINE READ_INFO()
!
        USE VARIAVEIS
        IMPLICIT NONE
!
        INTEGER            :: ISTAT,IN_FILE,I
        CHARACTER(len=128) :: NAME
        CHARACTER(len=4)   :: EXT
        CHARACTER(len=5)   :: C
!
!! ABERTURA DO PRIMEIRO CAMPO !!!!!!!!!!!!!!!!!!!!!!!!!!!
        WRITE(C,113)0
        C=ADJUSTL(C)
        EXT='.dat'
        IN_FILE=212
        NAME=TRIM(FILEIN)//TRIM(C)//TRIM(EXT)
        NAME=ADJUSTL(TRIM(NAME))
        !WRITE(*,111)0,NAME(1:LEN_TRIM(NAME))
        OPEN(UNIT=IN_FILE,FILE=NAME,ACTION='READ',IOSTAT=ISTAT)
        IF(ISTAT.NE.0)THEN
           WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',NAME
           STOP
        END IF
!
!! LEITURA DAS INFORMACOES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        READ(IN_FILE,*)FDX
        READ(IN_FILE,*)FDY
        READ(IN_FILE,*)NX
        READ(IN_FILE,*)NY
!        WRITE(*,10)FDX,FDY,NX,NY
!
        CLOSE(IN_FILE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
111     FORMAT('NAME OF INPUT FILE (RANDOM FIELD)',I5,': ',A,&
             /,'###### INFORMACOES ######')
113     FORMAT(I5)
10      FORMAT(/,'DIMENSION X =',F12.5,/,&
                 'DIMENSION Y =',F12.5,/,&
                 'N. OF ELEMENTS IN X =',I10,/,&
                 'N. OF ELEMENTS IN Y =',I10)
!
      END SUBROUTINE READ_INFO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ROTINA PARA LEITURA DO CAMPO DE PERMEABILIDADES !!!!!!
!========================================================
!                                              
      SUBROUTINE READPERM(K)
!
        USE VARIAVEIS
        IMPLICIT NONE
!
        INTEGER :: I,J,K,NLINE,NFLAG
        INTEGER :: IN_FILE,ISTAT
        CHARACTER(len=128) :: NAME
        CHARACTER(len=4)   :: EXT
        CHARACTER(len=1)   :: tipo
        CHARACTER(len=5)   :: C
        DOUBLE PRECISION   :: LIXO
!
!! ABERTURA DO CAMPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        WRITE(C,113)K
        C=ADJUSTL(C)
        EXT='.dat'
        IN_FILE=213
        NAME=TRIM(FILEIN)//TRIM(C)//TRIM(EXT)
        NAME=ADJUSTL(TRIM(NAME))
        WRITE(*,111)K,NAME(1:LEN_TRIM(NAME))
        OPEN(UNIT=IN_FILE,FILE=NAME,IOSTAT=ISTAT)
        IF(ISTAT.NE.0)THEN
           WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',NAME
           STOP
        END IF
!
!! LEITURA DAS INFORMACOES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        READ(IN_FILE,*)FDX
        READ(IN_FILE,*)FDY
        READ(IN_FILE,*)NX
        READ(IN_FILE,*)NY
        !WRITE(*,10)FDX,FDY,NX,NY
        READ(IN_FILE,*)NTYPE
        READ(IN_FILE,*)FBETA
        READ(IN_FILE,*)LIXO
        READ(IN_FILE,*)LIXO
!
!! LEITURA DOS CAMPOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO J=1,NY
!
           READ(IN_FILE,*) NLINE
           IF(NLINE+1.NE.J)THEN
              WRITE(*,*) 'Erro na leitura do campo de PERMeabilidade',NLINE
              STOP
           END IF
!
           READ(IN_FILE,*) (PERM(I+(J-1)*NX),I=1,NX)
!      
           read(IN_FILE,*) NFLAG
           IF(NFLAG.ne.192837465)THEN
              WRITE(*,*) 'Erro na leitura do campo de PERMeabilidade, nflag'
              STOP
           END IF
!      
        END DO    ! NY     
!     
        CLOSE(IN_FILE)
!     calculando a PERMeabilidade lognormal
!
!
        tipo='Y'
        call estatistica(PERM,NX,NY,tipo)
!
        do j=1,NY*NX
           PERM(j)=kg*dexp(rho*PERM(j))
        end do
!
111     FORMAT('######################################',/, &
               '######################################',/, &
               'NAME OF INPUT FILE (RANDOM FIELD)',I5,': ',A)
113     FORMAT(I5)
115     FORMAT(                                      &
       '######################################',/, &
       'PROBLEMA NO TAMANHO DO VETOR PARA A   ',/, &
       'LEITURA DO CAMPO DE PERMEABILIDADES   ',/, &
       'NUMERO DO CAMPO:',I5,/,                    &
       '######################################',/)
10      FORMAT(/,'DIMENSION X =',F12.5,/,&
                 'DIMENSION Y =',F12.5,/,&
                 'N. OF ELEMENTS IN X =',I10,/,&
                 'N. OF ELEMENTS IN Y =',I10)
!
      END SUBROUTINE READPERM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ROTINA PARA CALCULAR MEDIA E VARIANCIA !!!!!!!!!!!!!!!
!========================================================
      subroutine estatistica(prm,nmx,nmy,sc)
!
!     calculando a media e variancia
        USE VARIAVEIS
        implicit none
!
        real(8), dimension(*) :: prm
        integer :: nmx,nmy,i,j
        real(8) :: xlx,xly,xmm,XMAX,XMIN
        character*1 sc
!
        MED  = 0.d0
        xmm = 0.d0
        VAR  = 0.d0
        XMAX = -1.d14
        XMIN = 1.d14
!
        do i=1,nmy*nmx
           MED = MED+prm(i)
           xmm= xmm+prm(i)*prm(i)
           if(prm(i).gt.XMAX) XMAX = prm(i)
           if(prm(i).lt.XMIN) XMIN = prm(i)
        end do
        MED  = MED/(nmx*nmy)
        xmm = xmm/(nmx*nmy)
        VAR  = xmm - MED*MED
!
!
        write(*,123)sc,MED,VAR,XMAX,XMIN
123     FORMAT(A,                          &
             & ' ######################### ',/,  &
             & 'MEDIA     =',f10.5,/,            &
             & 'VARIANCIA =',f10.5,/,            &
             & 'MAXIMO    =',f10.5,/,            &
             & 'MINIMO    =',f10.5,/,            &
             & '##############################')
!
      end subroutine estatistica
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ROTINA PARA TOMADA DE MEDIA DO CAMPO - UPSCALING !!!!!
!========================================================
      SUBROUTINE UPSCALE(TIPO)
!
        USE VARIAVEIS
        IMPLICIT NONE
!
        INTEGER :: TIPO,I,J,EL,NO,K,M,L
        INTEGER :: MX,MY,N,LOCAL
        DOUBLE PRECISION :: AUX
        REAL             :: TIME,T_START,T_FINAL
        REAL,DIMENSION(2):: TARRAY
        DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:) :: MAT
!
        MX=NX/NXG
        MY=NY/NYG
        N=MX*MY
        AUX=1.D0/N
        ALLOCATE(MAT(N))
        MAT=0.D0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        CALL ETIME(TARRAY,T_START)
!        T_START = TARRAY(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        PERMG =0.d0
        PERMGT=0.d0
!
!! MEDIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(TIPO.EQ.1)THEN
           EL=0d0
           DO J=0,NYG-1
              DO I=0,NXG-1
                 EL=EL+1
                 DO K=0,MY-1
                    DO M=0,MX-1
                       NO=1+MX*I+MY*NX*J
                       LOCAL=NO+M+K*NX
                       PERMGT(EL)=PERMGT(EL)+AUX*PERM(LOCAL)
                    END DO
                 END DO
                 PERMGT(EL)=(LOG(PERMGT(EL))-LOG(KG))/RHO
              END DO
           END DO
        END IF
!! DURLOFSKY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(TIPO.EQ.0)THEN
           L=0D0
           DO J=0,NYG-1
              DO I=0,NXG-1
                 EL=0D0
                 DO K=0,MY-1
                    DO M=0,MX-1
                       EL=EL+1
                       NO=1+MX*I+MY*NX*J
                       LOCAL=NO+M+K*NX
                       MAT(EL)= PERM(LOCAL)
                    END DO
                 END DO
                 L=L+1
                 CALL BLOCOPOISSON2(MAT,MX,MY)
                 PERMGT(L)=(LOG(U)-LOG(KG))/RHO
              END DO
           END DO
        END IF
        DEALLOCATE(MAT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        CALL ETIME(TARRAY,T_FINAL)
!        T_FINAL = TARRAY(1)
!        TIME = T_FINAL-T_START
!        WRITE(*,334)TIME,TARRAY(1),TARRAY(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
334     FORMAT('TIME TO UPSCALING(s) =',F12.2,/,&
               'USER TIME        (s) =',F12.2,/,&
               'SYSTEM TIME      (s) =',F12.2,/)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END SUBROUTINE UPSCALE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CLEAR(VEC,NN)
        INTEGER, INTENT(IN):: NN
        INTEGER            :: I,J
        DOUBLE PRECISION   :: VEC(1:NN)
!
!$OMP PARALLEL PRIVATE(I)
!$OMP DO
        DO I=1,NN
           VEC(I)=0.d0
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
      END SUBROUTINE CLEAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PRINT_OLD(NAME,VEC,NNX,NNY,K)
!
        USE VARIAVEIS
        IMPLICIT NONE
!
        INTEGER            :: J,I,BAND,OUT_FILE,MMY,K
        INTEGER            :: NNX,NNY
        CHARACTER(len=128) :: NAME
        CHARACTER(len=128) :: NAMEL
        CHARACTER(len=4)   :: EXT
        CHARACTER(len=5)   :: C
        DOUBLE PRECISION   :: VEC(NNX*NNY)
!
        BAND=192837465
        OUT_FILE=321
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        WRITE(C,113)K
        C=ADJUSTL(C)
        EXT='.dat'
        NAMEL=TRIM(NAME)//TRIM(C)//TRIM(EXT)
        NAMEL=ADJUSTL(TRIM(NAMEL))
        WRITE(*,111)K,NAMEL(1:LEN_TRIM(NAMEL))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        OPEN(UNIT=OUT_FILE,FILE=NAMEL,&
             STATUS='REPLACE',ACTION='WRITE')
        WRITE(OUT_FILE,103)FDX
        WRITE(OUT_FILE,103)FDY
        WRITE(OUT_FILE,100)NNX
        WRITE(OUT_FILE,100)NNY
        WRITE(OUT_FILE,100)NTYPE
        WRITE(OUT_FILE,103)FBETA
        WRITE(OUT_FILE,100)2
        WRITE(OUT_FILE,100)2
!      
        DO J=1,NNY
           WRITE(OUT_FILE,100)J-1
           WRITE(OUT_FILE,101)(VEC(I+(J-1)*NNX),I=1,NNX)
           WRITE(OUT_FILE,102)BAND
        ENDDO
        call estatistica(VEC,NNX,NNY,'K')
!
        CLOSE(OUT_FILE)
!
100     FORMAT(I7)
101     FORMAT(10000F12.4)
102     FORMAT(I9)
103     FORMAT(1F12.4)
111     FORMAT('######################################',/, &
               '######################################',/, &
               'NAME OF OUTPUT FILE (RANDOM FIELD)',I5,': ',A)
113     FORMAT(I5)
!
      END SUBROUTINE PRINT_OLD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE TOMBA()
!
        USE VARIAVEIS
        IMPLICIT NONE
!
        INTEGER :: I,J,K,N,M
!
        K=0d0
        N=(NY-1)*NX+1
        DO I=0,NX-1
           DO J=0,NY-1
              K=K+1
              M=N-J*NX
              PERMT(K)=PERM(M)
           END DO
           N=N+1
        END DO
!
        DO I=1,NX*NY
           PERM(I)=PERMT(I)
        END DO
        M=NX
        NX=NY
        NY=M
        M=NXG
        NXG=NYG
        NYG=M
!
      END SUBROUTINE TOMBA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE DTOMBA()
!
        USE VARIAVEIS
        IMPLICIT NONE
!
        INTEGER :: I,J,K,N,M
!
        M=NXG
        NXG=NYG
        NYG=M
        K=0d0
        N=NYG
        DO I=0,NYG-1
           DO J=0,NXG-1
              K=K+1
              M=N+J*NYG
              PERMG(K)=PERMGT(M)
           END DO
           N=N-1
        END DO
!
        DO I=1,NXG*NYG
           PERMGT(I)=PERMG(I)
        END DO
!
      END SUBROUTINE DTOMBA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE BLOCOPOISSON(MAT,MX,MY)
!
        USE VARIAVEIS
        IMPLICIT NONE
!
        INTEGER :: NT,NB,I,J,IP,IM,JP,JM,K
        INTEGER :: KIM,KIP,KJM,KJP,MX,MY
        DOUBLE PRECISION :: XL,YL,HX,HY,Q,RXY,RYX,HXY
        DOUBLE PRECISION :: PLEFT,PRIGHT,XK,XKIM,XKIP
        DOUBLE PRECISION :: XKJM,XKJP,SUM,AUX
        DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:) :: A
        DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:)   :: F
        INTEGER         , ALLOCATABLE,DIMENSION(:)   :: IPIV
        INTEGER          :: LDA,LDB,INFO,NRHS
        DOUBLE PRECISION :: MAT(MX*MY)
        EXTERNAL         :: DGESV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DIMENSOES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        XL =1.D0
        YL =1.D0
        NT =MX*MY
        HX =XL/MX
        HY =YL/MY
        HXY=HX*HY
        RXY=2.D0*HX/HY
        RYX=2.D0*HY/HX
        NB =1
        Q  =0.D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INICIALIZANDO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        ALLOCATE(A(NT,NT))
        ALLOCATE(F(NT))
        ALLOCATE(IPIV(NT))
!
        A   =0.d0
        F   =0.d0
        IPIV=0.d0
!
        PLEFT =1.d0
        PRIGHT=0.d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MONTAGEM DO SISTEMA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        DO J=1,MY
           JP=J+1
           JM=J-1
!
           DO I=1,MX
!
              IP=I+1
              IM=I-1
!
              K=I+MX*(J-1)
              KIP=K+1 
              KIM=K-1
              KJP=K+MX
              KJM=K-MX
!
!     coNdIcoes
!
              IF(I.EQ. 1)  KIM = K
              IF(I.EQ.MX) KIP = K
!
!     FluXo Nulo
!
              IF(J.EQ.  1)KJM = K
              IF(J.EQ.MY)KJP = K
!
!     PERMGeAbIlIdAdes
!
              XK  =MAT(K)
              XKIP=MAT(KIP)
              XKIM=MAT(KIM)
              XKJP=MAT(KJP)
              XKJM=MAT(KJM)
!
!     MedIA HARMoNIcA
!
              XKIP=RYX*(XKIP*XK)/(XKIP+XK)
              XKIM=RYX*(XKIM*XK)/(XKIM+XK)
              XKJP=RXY*(XKJP*XK)/(XKJP+XK)
              XKJM=RXY*(XKJM*XK)/(XKJM+XK)
!
!     dIAgoNAl
!
              A(K,K) = A(K,K)-(XKJP+XKJM+XKIP+XKIM)
!
!     vIzINHos eM X
!
              A(KIP,K) = A(KIP,K)+XKIP
              A(KIM,K) = A(KIM,K)+XKIM
!
              IF(I.EQ.1) THEN
                 A(KIM,K) = A(KIM,K)-2.d0*XKIM
                 F(K)=F(K)-2.d0*XKIM*PleFt
              END IF
!
              IF(I.EQ.MX) THEN
                 A(KIP,K) = A(KIP,K)-2.d0*XKIP
                 F(K)=F(K)-2.d0*XKIP*PRIgHt
              END IF
!
!     vIzINHos eM Y
!
              A(K,KJP) = A(K,KJP)+XKJP
              A(K,KJM) = A(K,KJM)+XKJM
!
!     teRMo de FoNTe
!
              F(K)=F(K)+HXY*Q
!
           END DO    ! I
!
        END DO       ! J
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESOLUCAO DO SISTEMA LINEAR !!!!!!!!!!!!!!!!!!!!!!!!!!!
        NRHS=1
        LDA=NT
        LDB=MAX(1,NT)
        CALL DGESV(NT,NRHS,A,LDA,IPIV,F,LDB,INFO)
        !IF(INFO.EQ.0)WRITE(*,*)'successful exit'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INTEGRAL NO BORDO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        U=0.D0
        AUX=2.D0/HX
        DO J=1,MY
           I=J*MX
           U=U+F(I)*AUX*MAT(I)*HY
        END DO
        DEALLOCATE(A,F,IPIV)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    END SUBROUTINE BLOCOPOISSON
