 !!TODO extend to reading random buffer.cpp? 
!! TODO explist casting in Fortran?

!! http://stackoverflow.com/questions/2757424/discrepancy-between-the-values-computed-by-fortran-and-c

!!  -fdefault-real-8 can be used when compiling with gfortran in order to use double precision floating point literals
 
!! http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
!! http://softwareengineering.stackexchange.com/questions/62948/what-can-be-done-to-programming-languages-to-avoid-floating-point-pitfalls

PROGRAM F_test

   INTEGER*4, parameter :: test = 7
   INTEGER*4, parameter:: MUL = 1000
   
   INTEGER*4, parameter :: a = 5;
   INTEGER*4, parameter :: c  = 4; 
   
   REAL*4, parameter :: bF  = 1.012003;
   REAL*4, parameter :: dF  = 2.320421;
   REAL*8, parameter :: bD  = 1.012003;
   REAL*8, parameter :: dD  = 2.320421;


   REAL*4 :: FDiv=0.0
   REAL*4 :: FMul=0.0 
   REAL*4 :: FAdd=0.0 
   REAL*4 ::FSub=0.0  !signle precsion
   REAL*8 :: DDiv = 0.0
   REAL*8 ::  DMul = 0.0 
   REAL*8 :: DSub=0.0 
   REAL*8 :: DAdd =0.0 !double precsion
   INTEGER*4 :: Fcounter=0
   INTEGER*4 :: Dcounter=0
   REAL*4 :: FResult = 0
   REAL*8 :: DResult = 0

   OPEN(UNIT=11,ACCESS='SEQUENTIAL',FORM='FORMATTED',  &
        FILE='F_test.txt') 

   OPEN(UNIT=12,ACCESS='SEQUENTIAL',FORM='FORMATTED',  &
        FILE='F_P_test.txt') 

   WRITE(11,*) 'Integer test number is ', test
   WRITE(11,*) ''

   WRITE(11,*) 'Floats: precsion set to 7 decimal points ' 

   FResult=REAL(test)
   DO i=1000, 1, -1
      FDiv = FResult/REAL(a)
      FSub = FDiv - bF 
      FMul = FSub * REAL(c)
      FAdd = FMul + dF
      FResult = FAdd
      Fcounter = Fcounter + 1
   END DO
   
   WRITE(11,*) 'After ', Fcounter , 'iterations'
   
   WRITE(11,89) 'the result is ', FResult
89 FORMAT(A32,F16.7)
   WRITE(12,88) '', FResult
88 FORMAT(A32,F16.7)

   WRITE(11,*) 'Doubles: precsion set to 16 decimal points ' 

   DResult=DBLE(test)
   DO i=1000, 1, -1
      DDiv = DResult/DBLE(a)
      DSub = DDiv - bD 
      DMul = DSub * DBLE(c)
      DAdd = DMul + dD 
      Dcounter = Dcounter + 1
      DResult = DAdd
   END DO
   
   WRITE(11,*) 'After ', Dcounter , 'iterations' 
   WRITE(11, 90) 'the result is ', DResult
90 FORMAT(A32,F32.16)
   WRITE(12, 91) '', DResult
91 FORMAT(A32,F32.16)

   WRITE(11,*) ''
   WRITE(11,*) ''

   WRITE(11,*) 'Integer test number is ', test, 'multiplication factor is ', MUL 
   WRITE(11,*) ''

   WRITE(11,*) 'Floats: precsion set to 7 decimal points ' 

   FResult=REAL(test)*MUL
   Fcounter=0
   DO i=1000, 1, -1
      FDiv = FResult/REAL(a)
      FSub = FDiv - bF*MUL
      FMul = FSub * REAL(c)
      FAdd = FMul + dF*MUL
      Fcounter = Fcounter + 1
      FResult = FAdd
   END DO
    FResult = FResult/MUL
   WRITE(11,*) 'After ', Fcounter , 'iterations'
   
   WRITE(11,89) 'the result is ', FResult

   WRITE(12,88) '', FResult


   WRITE(11,*) 'Doubles: precsion set to 16 decimal points ' 

   DResult=DBLE(test)*MUL
   Dcounter=0
   DO i=1000, 1, -1
      DDiv = DResult/DBLE(a)
      DSub = DDiv - bD*MUL
      DMul = DSub * DBLE(c)
      DAdd = DMul + dD*MUL
      Dcounter = Dcounter + 1
      DResult = DAdd
   END DO
   DResult = DResult/MUL
   WRITE(11,*) 'After ', Dcounter , 'iterations' 
   WRITE(11, 90) 'the result is ', DResult

   WRITE(12, 91) '', DResult

   
   CLOSE (11)
   CLOSE(12)
END PROGRAM F_test
