!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------
!------------------------------------------------------------
!                 2D RANDOM FIELDS GENERATOR
!                USING KARHUNEN-LOEVE METHOD
! #################  CONDITIONAL VERSION  ###################
!
! NAME: MARCIO RENTES BORGES
! DATA: 03/11/2011
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------
!     Metodo baseado na expansao de Karhunen-Loeve         
!     
!     ESTA VERSAO PERMITE CONDICIONAR VALORES EM ALGUNS
!                   PONTOS DO DOMINIO
!                                               
!    
!------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MODULO COM AS VARIAVEIS GLOBAIS DO PROGRAMA !!!!!!!!!!!!!!
MODULE VARIAVEIS
  CHARACTER(LEN=128):: FILE_IN,FILE_OUT,FILE_COND
  CHARACTER(LEN=128):: FILE_VAL,FILE_VET,FILE_THE
  DOUBLE PRECISION  :: FDX,FDY,FDZ,BASE,SIG
  DOUBLE PRECISION  :: FBETA,FVAR,FMEAN,DIMX,DIMY
  DOUBLE PRECISION  :: VARIAN,XMEAN,XVAR,FATCOND
  INTEGER           :: NX,NY,NFLAG,NERROR,INIF,NFILES
  INTEGER           :: NFILE,NDIM,NS
  LOGICAL           :: TFCOND
  INTEGER           :: NCOND,MKL,NPROPOSAL
  INTEGER         , ALLOCATABLE,DIMENSION(:)   :: NSEED,PCOND
  DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:) :: PSI,VET,AVET
  DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:)   :: X,Y,AVAL,THETA
  INTEGER         , ALLOCATABLE,DIMENSION(:,:) :: NVET
END MODULE VARIAVEIS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM MAIN
!
  USE VARIAVEIS
  USE RANDOM
!
  IMPLICIT NONE
! LIST OF LOCAL VARIABLES
  INTEGER                    :: I,J,ERROR_READ,NUM_ELEM
  INTEGER, ALLOCATABLE       :: SEED(:)
  DOUBLE PRECISION, EXTERNAL :: POWERLAW,EXPONENTIAL
  REAL                       :: TIME
  REAL, DIMENSION(2)         :: TARRAY
!
  ERROR_READ = 0
! READ INPUT DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL READ_DATA(ERROR_READ)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MEMORY ALLOCATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NUM_ELEM = NX*NY
  ALLOCATE(PSI(0:NX-1,0:NY-1))
  ALLOCATE(X(0:NX))
  ALLOCATE(Y(0:NY))
  ALLOCATE(AVAL(1:MKL))
  ALLOCATE(AVET(1:NUM_ELEM,1:MKL))
  ALLOCATE(PCOND(1:MKL))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CLEAR VECTORS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL CLEAR2D(PSI,NX,NY)
  CALL CLEAR(X,NX)
  CALL CLEAR(Y,NY)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONSTRUCT VECTORS OF POSITIONS !!!!!!!!!!!!!!!!!!!!!!!!
!$OMP PARALLEL PRIVATE(I)
!$OMP DO
  DO I=0,NX
     X(I)=(I*FDX)+FDX*.5
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
!$OMP PARALLEL PRIVATE(J)
!$OMP DO
  DO J=0,NY
     Y(J)=(J*FDY)+FDY*.5
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GETS THE SEEDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL RANDOM_SEED(put=NSEED)      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RANDOM FIELDS GENERATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WRITE(*,*)'LOADING EINGENPAIRS'
  CALL ETIME(TARRAY,TIME)
!
  CALL LOAD_AUTOV(AVAL,AVET,NUM_ELEM,MKL)
!
  CALL LOAD_THETA()
!
  CALL ETIME(TARRAY,TIME)
  WRITE(*,333)TIME,TARRAY(1),TARRAY(2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RANDOM CONDITIONING AND GENERATION !!!!!!!!!!!!!!!!!!!!
  WRITE(*,*)'CONSTRUCTING       '
  CALL ETIME(TARRAY,TIME)
!
  CALL KLCOND()
!
  CALL ETIME(TARRAY,TIME)
  WRITE(*,331)TIME,TARRAY(1),TARRAY(2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IMPRESSAO DA SEMENTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CALL PRINTSEED()
! MEMORY FREE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DEALLOCATE(AVAL,AVET,PSI)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NERROR = 1
!
333 FORMAT('TIME TO LOAD EIGENPAIRS (s) =',F8.2,/,&
           'USER TIME               (s) =',F8.2,/,&
           'SYSTEM TIME             (s) =',F8.2,/)
!
331 FORMAT('TIME TO GENERATE        (s) =',F8.2,/,&
           'USER TIME               (s) =',F8.2,/,&
           'SYSTEM TIME             (s) =',F8.2,/)
!
END PROGRAM MAIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ THE INPUT DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE READ_DATA(NERROREAD)
!
        USE VARIAVEIS
        USE RANDOM
!
        INTEGER :: IN_FILE,ISTAT
        INTEGER, INTENT(OUT) :: NERROREAD
!
        IN_FILE = 123
        FILE_IN = './entrada.in'
        OPEN(UNIT=IN_FILE,FILE=FILE_IN,STATUS='UNKNOWN',&
             FORM='FORMATTED',IOSTAT=ISTAT)
        IF(ISTAT.NE.0)THEN
           WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_IN
           STOP
        END IF
!
        READ(IN_FILE,8000)INIF,NFILES
        WRITE(*,8001)INIF,NFILES
        NFILE=NFILES-INIF+1
        WRITE(*,8002)NFILE
!
        READ(IN_FILE,1000)DIMX,DIMY
        WRITE(*,1001)DIMX,DIMY
!
        READ(IN_FILE,2000)NX,NY
        WRITE(*,2001)NX,NY
!
        FDX = DIMX/NX
        FDY = DIMY/NY
        WRITE(*,5001)FDX,FDY
!
        READ(IN_FILE,3000)NFLAG
!
        IF(NFLAG==2)THEN
           WRITE(*,6002)NFLAG
        END IF
        IF(NFLAG==1)THEN
           WRITE(*,6001)NFLAG
        END IF
!
        READ(IN_FILE,9000)FBETA,VARIAN
        IF(NFLAG==1)THEN
           WRITE(*,9001)FBETA,VARIAN
        END IF
        IF(NFLAG==2)THEN
           WRITE(*,9002)FBETA,VARIAN
        END IF
!
        CALL RANDOM_SEED(size=NS)
        ALLOCATE(NSEED(NS))
        READ(IN_FILE,97)(NSEED(I),I=1,NS)
        WRITE(*,96)NS,(NSEED(I),I=1,NS)
!
        READ(IN_FILE,3000)MKL
        WRITE(*,2005)MKL
!
        READ(IN_FILE,7000)FILE_OUT
        WRITE(*,7001)FILE_OUT
!
        READ(IN_FILE,7000)FILE_VAL
        WRITE(*,7002)FILE_VAL
!
        READ(IN_FILE,7000)FILE_VET
        WRITE(*,7003)FILE_VET
!
        READ(IN_FILE,7000)FILE_THE
        WRITE(*,7004)FILE_THE
!
        READ(IN_FILE,12)NPROPOSAL,SIG
        WRITE(*,13)NPROPOSAL,SIG
        IF(NPROPOSAL.EQ.1)WRITE(*,*)'RANDOM WALK'
        IF(NPROPOSAL.EQ.2)WRITE(*,*)'pCN'
        IF(NPROPOSAL.EQ.4)WRITE(*,*)'DE or DREAM'
!
        READ(IN_FILE,3000)NCOND
        WRITE(*,3002)NCOND
!
        ALLOCATE(VET(0:2,0:NCOND-1))
        DO I=0,NCOND-1
           READ(IN_FILE,1000)VET(0,I),VET(1,I),VET(2,I)
           WRITE(*,1002)I+1,VET(0,I),VET(1,I),VET(2,I)
           DO J=0,1
              VET(J,I)=VET(J,I)+1E-6
           ENDDO
        ENDDO
!
        NERROREAD = 1d0
        CLOSE(UNIT=IN_FILE)
!
        RETURN

96      FORMAT('NUMERO DE SEMENTES:',I2,/,'SEMENTES:',40I12)
97      FORMAT(40I12)
12      FORMAT(I7,E10.3)
13      FORMAT('TIPO DE PROPOSAL         =',I7,/,&
               'PARAMETRO DO RANDOM WALK =',E10.3)
1000    FORMAT(3F12.5)
1001    FORMAT('DIMENSION X =',F12.5,/,&
               'DIMENSION Y =',F12.5,/)
1002    FORMAT('POINT    =',I6,/,&
               'COORD. X =',F12.5,/,&
               'COORD. Y =',F12.5,/,&
               'VALUE    =',F12.5,/)
2000    FORMAT(3I7)
2001    FORMAT('N. OF ELEMENTS IN X =',I7,/,&
               'N. OF ELEMENTS IN Y =',I7,/)
2005    FORMAT('M (terms) =',I7,/)
3000    FORMAT(I7)
3001    FORMAT('FLAG  =',I7)
3002    FORMAT('NCOND =',I7)
4000    FORMAT(3I12)
4001    FORMAT(/,'SEEDS =',I12,I12,I12,I12,I12,I12,I12,I12,/)
5001    FORMAT('DELTA X =',F12.4,/,&
               'DELTA Y =',F12.4,/)
6001    FORMAT(I7,': EXPONENTIAL')
6002    FORMAT(I7,': FRACTAL')
7000    FORMAT(A)
7001    FORMAT('NAME OF OUTPUT FILE = ',A)
7002    FORMAT('NAME OF INPUT AUTOVALUES FILE  = ',A)
7003    FORMAT('NAME OF INPUT AUTOVECTORS FILE = ',A)
7004    FORMAT('NAME OF INPUT THETA (NORMAL RV)= ',A)
8000    FORMAT(2I7)
8001    FORMAT('N. OF INTIAL FILE =',I7,/,&
               'N. OF FINAL FILE  =',I7)
8002    FORMAT('NUMBER OF FILES   =',I7,/)
9000    FORMAT(2F12.4)
9002    FORMAT('HURST COEF.=',F12.4,/,&
               'VARIANCE   =',F12.4)
9001    FORMAT('CORRELATION LENGTH =',F12.4,/,&
               'VARIANCE   =',F12.4) 
!
      END SUBROUTINE READ_DATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CLEAR2D(VEC,NNX,NNY)
        INTEGER, INTENT(IN):: NNX,NNY
        INTEGER            :: I,J
        DOUBLE PRECISION   :: VEC(0:NNX-1,0:NNY-1)
!
!$OMP PARALLEL PRIVATE(I,J)
!$OMP DO
        DO I=0,NNX-1
           DO J=0,NNY-1
             VEC(I,J)=0.d0
           ENDDO
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
      END SUBROUTINE CLEAR2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CLEAR(VEC,NN)
        INTEGER, INTENT(IN):: NN
        INTEGER            :: I,J
        DOUBLE PRECISION   :: VEC(0:NN)
!
!$OMP PARALLEL PRIVATE(I)
!$OMP DO
        DO I=0,NN
          VEC(I)=0.d0
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
      END SUBROUTINE CLEAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE PRINT_OLD(VEC,VX,VY,NNX,NNY,NAME)
      USE VARIAVEIS
      INTEGER,            INTENT(IN) :: NNX,NNY
      DOUBLE PRECISION,   INTENT(IN) :: VEC(0:NNX-1,0:NNY-1),VX(0:NNX),VY(0:NNY)
      CHARACTER(LEN=128)             :: NAME
      INTEGER                        :: J,I,BAND,OUT_FILE,MMY
!
      BAND=192837465
      OUT_FILE=321
!
      OPEN(UNIT=OUT_FILE,FILE=NAME,STATUS='REPLACE',ACTION='WRITE')
      WRITE(OUT_FILE,103)DIMX
      WRITE(OUT_FILE,103)DIMY
      WRITE(OUT_FILE,100)NX
      WRITE(OUT_FILE,100)NY
      WRITE(OUT_FILE,100)NFLAG
      WRITE(OUT_FILE,103)FBETA
      WRITE(OUT_FILE,100)2
      WRITE(OUT_FILE,100)2
      
      MMY=NNY-1
      DO J=0,MMY
         WRITE(OUT_FILE,100)J
         WRITE(OUT_FILE,101)(VEC(I,J),I=0,NNX-1)
         WRITE(OUT_FILE,102)BAND
      ENDDO
!
      CLOSE(OUT_FILE)
!
100   FORMAT(I7)
101   FORMAT(10000F12.4)
102   FORMAT(I9)
103   FORMAT(1F12.4)
!
      END SUBROUTINE PRINT_OLD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PRINTASC(VEC,VX,VY,NNX,NNY,NAME)
        INTEGER,            INTENT(IN) :: NNX,NNY
        DOUBLE PRECISION,   INTENT(IN) :: VEC(0:NNX-1,0:NNY-1),VX(0:NNX),VY(0:NNY)
        CHARACTER(LEN=128)             :: NAME
        INTEGER                        :: J,I,BAND,OUT_FILE
!
        BAND=192837465
        OUT_FILE=321
!
        OPEN(UNIT=OUT_FILE,FILE=NAME,STATUS='REPLACE',ACTION='WRITE')
          DO J=0,NNY-1
             WRITE(OUT_FILE,100)J+1
             DO I=0,NNX-1
                WRITE(OUT_FILE,101)VX(I),VY(J),VEC(I,NNY-(1+J))
             ENDDO
             WRITE(OUT_FILE,102)BAND
          ENDDO
        CLOSE(OUT_FILE)
!
 100  FORMAT(I7)
 101  FORMAT(F12.4,2X,F12.4,2X,F12.4)
 102  FORMAT(I9)
      END SUBROUTINE PRINTASC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PRINTBIN(VEC,VX,VY,NNX,NNY,NAME)
        INTEGER,            INTENT(IN) :: NNX,NNY
        DOUBLE PRECISION,   INTENT(IN) :: VEC(0:NNX-1,0:NNY-1),VX(0:NNX),VY(0:NNY)
        CHARACTER(LEN=128)             :: NAME
        INTEGER                        :: J,I,BAND,OUT_FILE
!
        BAND=192837465
        OUT_FILE=21
!
        OPEN(UNIT=OUT_FILE,FILE=NAME,STATUS='REPLACE',ACTION='WRITE')!,FORM='UNFORMATTED')
        WRITE(OUT_FILE,*)VEC
        CLOSE(OUT_FILE)
!
 100  FORMAT(I7)
 101  FORMAT(F12.4,2X,F12.4,2X,F12.4)
 102  FORMAT(I9)
      END SUBROUTINE PRINTBIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION MENOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER FUNCTION MENOR(A,B)
        INTEGER :: A,B
!
        IF(A<B)THEN
           MENOR = A
        ELSE
           MENOR = B
        END IF
!
      END FUNCTION MENOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION MAIOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER FUNCTION MAIOR(A,B)
        INTEGER :: A,B
!
        IF(A>B)THEN
           MAIOR = A
        ELSE
           MAIOR = B
        END IF
!
      END FUNCTION MAIOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION NARRED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER FUNCTION NARRED(LD,C)
!
        DOUBLE PRECISION, INTENT(IN) :: C,LD
!        DOUBLE PRECISION             :: AUX,AUX2
!        INTEGER                      :: N
!!
!        N=INT(LD/C)
!        AUX=(N*C)
!        AUX2=AUX+C*0.5
!        write(*,*)LD,REAL(NINT(LD))
!        IF(LD.LT.AUX2)THEN
!           NARRED=AUX
!        ELSE
!           NARRED=AUX+C
!        END IF
        NARRED = (NINT(LD))
!
      END FUNCTION NARRED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION MEAN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DOUBLE PRECISION FUNCTION VMEAN(VEC,NNX,NNY)
!
        DOUBLE PRECISION, INTENT(IN) :: VEC(0:NNX-1,0:NNY-1)
        DOUBLE PRECISION             :: SUM
        INTEGER, INTENT(IN)          :: NNX,NNY
        INTEGER                      :: I,J
!
        SUM=0.0
        DO I=0,NNX-1
           DO J=0,NNY-1
              SUM=SUM+VEC(I,J)
           ENDDO
        ENDDO
!
        VMEAN=SUM/DBLE(NNX*NNY)
!
      END FUNCTION VMEAN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION VARIANCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DOUBLE PRECISION FUNCTION VAR(VEC,NNX,NNY)
!
        DOUBLE PRECISION, INTENT(IN) :: VEC(0:NNX-1,0:NNY-1)
        DOUBLE PRECISION             :: SUM,AUX
        DOUBLE PRECISION, EXTERNAL   :: VMEAN
        INTEGER, INTENT(IN)          :: NNX,NNY
        INTEGER                      :: I,J
!
        AUX=VMEAN(VEC,NNX,NNY)
        SUM=0.0
        DO I=0,NNX-1
           DO J=0,NNY-1
              SUM=SUM+VEC(I,J)*VEC(I,J)
           ENDDO
        ENDDO
!
        VAR=SUM/DBLE(NNX*NNY)-AUX*AUX
!
      END FUNCTION VAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION VARIANCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DOUBLE PRECISION FUNCTION VARIANCIA(VEC,NNX)
!
        DOUBLE PRECISION, INTENT(IN) :: VEC(NNX)
        DOUBLE PRECISION             :: SUM,ME
        INTEGER, INTENT(IN)          :: NNX
        INTEGER                      :: I
!
        SUM=0.0
        ME =0.0
        DO I=1,NNX
           ME =ME+VEC(I)
        ENDDO
        ME=ME/DBLE(NNX-1)
!
        DO I=1,NNX
           SUM=SUM+(VEC(I)-ME)**2
        END DO
        VARIANCIA=SUM/DBLE(NNX-1)
!
      END FUNCTION VARIANCIA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE LOAD_AUTOV(VAL,VE,NELEM,M)
!
      USE VARIAVEIS
      IMPLICIT NONE
!
      DOUBLE PRECISION :: VAL(1:M)
      DOUBLE PRECISION :: VE(1:NELEM,1:M)
      INTEGER          :: NELEM,M
      INTEGER          :: ISTAT,IN_FILEA,IN_FILEV
      INTEGER          :: I,J
!
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IN_FILEA = 321
      OPEN(UNIT=IN_FILEA,FILE=FILE_VAL,STATUS='UNKNOWN', &
 FORM='FORMATTED',IOSTAT=ISTAT)
      IF(ISTAT.NE.0)THEN
         WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_VAL
         STOP
      END IF
!
      IN_FILEV = 322
      OPEN(UNIT=IN_FILEV,FILE=FILE_VET,STATUS='UNKNOWN', &
 FORM='FORMATTED',IOSTAT=ISTAT)
      IF(ISTAT.NE.0)THEN
         WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_VET
         STOP
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,M
         READ(IN_FILEA,111)VAL(I)
      ENDDO
!
      DO I=1,NELEM
         !READ(IN_FILEV,112)(VE(I,J),J=1,M)
         READ(IN_FILEV,*)(VE(I,J),J=1,M)
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111   FORMAT(e25.8)
112   FORMAT(40000e25.8)
!
      END SUBROUTINE LOAD_AUTOV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE KLCOND()
!
        USE VARIAVEIS
        USE RANDOM
        IMPLICIT NONE
!
        INTEGER                       :: PNODE(1:NCOND),K,ISTAT
        INTEGER                       :: NRHS,LDA,LDB,INFO,IN_FILE
        INTEGER                       :: IPIV(NCOND,1)
        DOUBLE PRECISION              :: B(NCOND,1),MAT( NCOND,NCOND)
        DOUBLE PRECISION              :: AUX1,VARR,AUXSQ
        DOUBLE PRECISION              :: XI(1:NX*NY)
        INTEGER                       :: I,J,M,MX,MY
        DOUBLE PRECISION              :: POSIX,POSIY
        DOUBLE PRECISION              :: AUX
        INTEGER                       :: MAX,MIN,MIN2
        CHARACTER(LEN=256)            :: NAME
        CHARACTER(LEN=4)              :: EXT
        CHARACTER(LEN=5)              :: C
        REAL                          :: TIME,T_START,T_FINAL
        REAL            , DIMENSION(2):: TARRAY
        DOUBLE PRECISION, EXTERNAL    :: MAIOR,MENOR
        DOUBLE PRECISION, EXTERNAL    :: VMEAN,VAR,VARIANCIA
        EXTERNAL                      :: DGESV
        DOUBLE PRECISION              :: WT
!
!        EXTERNAL :: ETIME
! SET VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IDENTIFICACAO DOS PONTOS NO VETOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL ETIME(TARRAY,T_START)
        T_START = TARRAY(1)
        DO M=0,NCOND-1
           K=0
           DO J=NY,1,-1
              POSIY = REAL(NY-J)*FDY
              DO I=1,NX
                 K=K+1
                 POSIX = REAL(I)*FDX
                 IF((VET(0,M).LT.POSIX.AND.VET(0,M).GT.POSIX-FDX).AND.&
                      (VET(1,M).LT.POSIY+FDY.AND.VET(1,M).GT.POSIY))THEN
                    PNODE(M+1)=K
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO I=1,NX*NY
           DO J=1,MKL
!              AVET(I,J)=SQRT(AVAL(J))*AVET(I,J)
              AVET(I,J)=AVET(I,J)
           ENDDO
        ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        AUX=0.D0
        DO I=1,MKL
           AUX1=AVET(1,I)
           AUX=AUX+0.5*AUX1*AUX1
        ENDDO
        AUXSQ=AUX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL ETIME(TARRAY,T_FINAL)
        T_FINAL = TARRAY(1)
        TIME = T_FINAL-T_START
        WRITE(*,334)TIME,TARRAY(1),TARRAY(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        MX=NX-1
        MY=NY-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        EXT='.dat'
        FILE_OUT=ADJUSTL(FILE_OUT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONDITIONING ==> BEST CONDITIONING !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(NCOND.GT.0)THEN
           FILE_COND = 'in/cond_seq.in'
           FILE_COND = ADJUSTL(TRIM(FILE_COND))
           INQUIRE(FILE=FILE_COND,EXIST=TFCOND)
           WRITE(*,*)'#########################'
           IF(TFCOND)THEN
              WRITE(*,*)'READING FILE:',ADJUSTL(TRIM(FILE_COND))
              IN_FILE = 22
!
              OPEN(UNIT=IN_FILE,FILE=FILE_COND,STATUS='OLD',&
                   ACTION='READ',FORM='FORMATTED',IOSTAT=ISTAT)
              IF(ISTAT.NE.0)THEN
                 WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_COND
                 STOP
              END IF
              DO I=1,MKL
                 READ(IN_FILE,*)PCOND(I)
              END DO
              CLOSE(IN_FILE)
           ELSE
              WRITE(*,*)'CONVENTIONAL CONDITIONING'
              DO I=1,MKL
                 PCOND(I)=I
              END DO
           END IF
           WRITE(*,*)'#########################'
        END IF
! MAIN LOOP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NUMBER OF FILES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        DO M=INIF,NFILES+INIF
!
           CALL ETIME(TARRAY,T_START)
           T_START = TARRAY(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GERACAO DA VA GAUSSIANA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           IF(NPROPOSAL.NE.4)THEN
              IF(NPROPOSAL.EQ.1)THEN
                 DO I=1,MKL
                    THETA(I)=THETA(I)+SIG*RANDOM_NORMAL()
                 ENDDO
! NORMALIZACAO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 AUX=SQRT(1.d0/VARIANCIA(THETA,MKL))
                 DO I=1,MKL
                    THETA(I)=AUX*THETA(I)
                 END DO
              END IF
              IF(NPROPOSAL.EQ.2)THEN
                 DO I=1,MKL
                    THETA(I)=SQRT(1.0-SIG*SIG)*THETA(I)+SIG*RANDOM_NORMAL()
!                         +MEANTHETA*(1.D0-SQRT(1.0-SIG*SIG))
                 END DO
              END IF
           END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONDICIONAMENTO DO CAMPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
           VARR=1.d0
           NRHS=1
           LDA =NCOND
           LDB =MAX(1,NCOND)
           IF(NCOND>0)THEN
              DO I=1,NCOND
                 B(I,1)=0.D0
              ENDDO
              DO I=1,NCOND
                 DO J=NCOND+1,MKL
                    B(I,1)=B(I,1)-&
                         AVET(PNODE(I),PCOND(J))*THETA(PCOND(J))
                 ENDDO
                 B(I,1)=B(I,1)+VET(2,I-1)+AUXSQ-VARR*0.5
                 DO J=1,NCOND
                    MAT(I,J)=AVET(PNODE(I),PCOND(J))
                 ENDDO
              ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! resolucao do sitema linear ! Aa=b !!!!!!!!!!!!!!!!!!!!!!!!!!!!
              CALL DGESV(NCOND,NRHS,MAT,LDA,IPIV,B,LDB,INFO)
              IF(INFO.EQ.0)THEN
                 WRITE(*,*)'successful exit'
              ELSE
                 CALL INFORMACAO(MAT,INFO)
                 STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              DO I=1,NCOND
                 THETA(PCOND(I))=B(I,1)
              ENDDO
           ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SAVE THE NEW THETA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           CALL WRITE_THETA()
!
           DO I=1,NX*NY
              AUX=0.D0
              DO J=1,MKL
                 AUX=AUX+AVET(I,J)*THETA(J)
              ENDDO
              XI(I)=0.5*VARR+AUX-AUXSQ
           ENDDO
!$OMP PARALLEL PRIVATE(I,J)
!$OMP DO
           K=0
           DO J=0,MY
              DO I=0,MX
                 K=K+1
                 PSI(I,J)=SQRT(VARIAN)*XI(K)
              ENDDO
           ENDDO
!$OMP END DO
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           WRITE(C,113)M
           C=ADJUSTL(C)
           NAME=TRIM(FILE_OUT)//TRIM(C)//TRIM(EXT)
           NAME=ADJUSTL(TRIM(NAME))
           WRITE(*,111)M,NAME(1:LEN_TRIM(NAME))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PRINT FIELD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	   IF(NFLAG==1)CALL CORRETOR_VAR(1.D0/CORRETOR)
!           CALL PRINTASC(PSI,X,Y,NX,NY,NAME)
           CALL PRINT_OLD(PSI,X,Y,NX,NY,NAME)
!           CALL PRINT_OLDUW(PSI,X,Y,NX,NY,NAME)
!
!          CALL PRINTBIN(PSI,X,Y,NX,NY,NAME)
           WRITE(*,112)VMEAN(PSI,NX,NY),VAR(PSI,NX,NY)
           CALL ETIME(TARRAY,T_FINAL)
           T_FINAL = TARRAY(1)
           TIME = T_FINAL-T_START
           WRITE(*,333)TIME,TARRAY(1),TARRAY(2)
        ENDDO
!
333     FORMAT('TIME TO GENERATE (s) =',F12.2,/,&
               'USER TIME        (s) =',F12.2,/,&
               'SYSTEM TIME      (s) =',F12.2,/)
334     FORMAT('TIME TO LOAD AUT (s) =',F12.2,/,&
               'USER TIME        (s) =',F12.2,/,&
               'SYSTEM TIME      (s) =',F12.2,/)
111     FORMAT('NAME OF OUTPUT FILE',I5,': ',A)
112     FORMAT('MEAN=',F12.4,'   VARIANCE=',F12.4)
113     FORMAT(I5)
!
      END SUBROUTINE KLCOND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE PRINTSEED()
!
      USE VARIAVEIS
      IMPLICIT NONE
!
      CHARACTER(LEN=128)             :: NAME
      INTEGER                        :: I,OUT_FILE
      INTEGER, DIMENSION(8)          :: SEED 
!
      OUT_FILE=341
      NAME='./out/seed.dat'
      NAME=ADJUSTL(TRIM(NAME))
      WRITE(*,111)NAME(1:LEN_TRIM(NAME))

      OPEN(UNIT=OUT_FILE,FILE=NAME,STATUS='REPLACE',ACTION='WRITE')
      CALL RANDOM_SEED(get=NSEED)
      WRITE(OUT_FILE,103)(SEED(I),I=1,8)
!
      CLOSE(OUT_FILE)
!
103   FORMAT(I10)
111   FORMAT('NAME OF SEED OUTPUT FILE:',A)
!
      END SUBROUTINE PRINTSEED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE LOAD_THETA()
!
        USE VARIAVEIS
        IMPLICIT NONE
!
        INTEGER          :: ISTAT,IN_FILE
        INTEGER          :: I,J
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ALOCACAO DA MEMORIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ALLOCATE(THETA(MKL))
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IN_FILE = 325
        OPEN(UNIT=IN_FILE,FILE=FILE_THE,STATUS='UNKNOWN', &
             FORM='FORMATTED',IOSTAT=ISTAT)
        IF(ISTAT.NE.0)THEN
           WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_THE
           STOP
        END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO I=1,MKL
           READ(IN_FILE,*)THETA(I)
        ENDDO
!
        CLOSE(UNIT=IN_FILE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111     FORMAT(e15.8)
!
      END SUBROUTINE LOAD_THETA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE WRITE_THETA()
!
        USE VARIAVEIS
        IMPLICIT NONE
!
        INTEGER            :: ISTAT,OUT_FILE
        INTEGER            :: I,J
        CHARACTER(LEN=128) :: FILETHE
        CHARACTER(LEN=5)   :: FEXP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO I=1,128-4
           FEXP = FILE_THE(I:I+4)
           IF(FEXP.EQ.'theta')GOTO 123
        END DO
        123 CONTINUE
        FILETHE = TRIM(ADJUSTL(FILE_THE(1:I+4)))//('new')//&
             TRIM(ADJUSTL(FILE_THE(I+5:128)))
        OUT_FILE = 326
        FILETHE=TRIM(ADJUSTL(FILETHE))
        OPEN(UNIT=OUT_FILE,FILE=FILETHE,STATUS='REPLACE',&
             ACTION='WRITE',IOSTAT=ISTAT)
        IF(ISTAT.NE.0)THEN
           WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILETHE
           STOP
        END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO I=1,MKL
           WRITE(OUT_FILE,111)THETA(I)
        ENDDO
!
        CLOSE(UNIT=OUT_FILE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111     FORMAT(e15.8)
!
      END SUBROUTINE WRITE_THETA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PRINT_OLDUW(VEC,VX,VY,NNX,NNY,NAME)
        USE VARIAVEIS
        INTEGER,            INTENT(IN) :: NNX,NNY
        DOUBLE PRECISION,   INTENT(IN) :: VEC(0:NNX-1,0:NNY-1),VX(0:NNX),VY(0:NNY)
        CHARACTER(LEN=128)             :: NAME
        INTEGER                        :: J,I,BAND,OUT_FILE,MMY
!
        OUT_FILE=321
!
        OPEN(UNIT=OUT_FILE,FILE=NAME,STATUS='REPLACE',ACTION='WRITE')
        WRITE(OUT_FILE,*)'# FIELD'
        WRITE(OUT_FILE,100)NX,NY
      
        DO J=0,NNY-1
           DO I=0,NNX-1
              WRITE(OUT_FILE,101)VEC(I,J)
           END DO
        END DO
!
        CLOSE(OUT_FILE)
!
100     FORMAT(I4,' x ',I4,' x 1')
101     FORMAT(e15.8)
!
      END SUBROUTINE PRINT_OLDUW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=============================================================================
!=============================================================================
!=============================================================================
      SUBROUTINE INFORMACAO(A,inform)
        USE VARIAVEIS, ONLY: NCOND
        IMPLICIT NONE
        DOUBLE PRECISION :: A(NCOND,NCOND)
        INTEGER          :: I,J,inform
        WRITE(*,*)'###########################################################'
        WRITE(*,*)'###########################################################'
        WRITE(*,*)'### PROBLEMA NA RESOLUCAO DO SISTEMA LINEAR - INFO= ',inform
        WRITE(*,*)'A positive value of INFO on return from an LAPACK routine'
        WRITE(*,*)'indicates a failure in the course of the algorithm. Common'
        WRITE(*,*)'causes are:'
        WRITE(*,*)'--> a matrix is singular (to working precision);'
        WRITE(*,*)'--> a symmetric matrix is not positive definite;'
        WRITE(*,*)'--> an iterative algorithm for computing eigenvalues or '
        WRITE(*,*)'    eigenvectors fails to converge in the permitted number '
        WRITE(*,*)'    of iterations.'
        WRITE(*,*)'###########################################################'
        WRITE(*,*)'### MATRIZ  ###############################################'
        DO I=1,NCOND
           WRITE(*,*)(A(I,J),J=1,NCOND)
        END DO
        WRITE(*,*)'###########################################################'
        WRITE(*,*)'DETALHES: http://www.netlib.org/lapack/lug/node138.html'
        WRITE(*,*)'###########################################################'
        WRITE(*,*)'###########################################################'
      END SUBROUTINE INFORMACAO
!=============================================================================
!=============================================================================
