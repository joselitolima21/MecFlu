PROGRAM TRABALHO
IMPLICIT NONE
	!VARIÁVEIS DO PROBLEMA
INTEGER n,i,j,G,c,o,ex,cap,mt,cas,casup,caspl,cash
REAL*8 k,r,x,l,TA,TB
REAL*8 oi,op,ap1,ap2,ap5,ae1,ae2,ae5,aw1,aw2,aw5,su1,su2,su5,sp1,sp2,sp5
REAL*8 q,z,rho,ux,D,F,phia,phib	
REAL*8 aplinha1,aplinha2,aplinha5,Bound1,Bound2,Bound5,Pe
	!ALOCÇAO DA MATRIZ
!REAL, DIMENSION(:,:), ALLOCATABLE :: A
REAL, DIMENSION(:), ALLOCATABLE :: U
REAL, DIMENSION(:), ALLOCATABLE :: S
REAL, DIMENSION(:), ALLOCATABLE :: V
	!APRESENTACAO!
PRINT*, "  ------------------------------------------------------ "
PRINT*, "  -------  Feito por Joselito Lima Reis Junior,  ------- "
PRINT*, "  -------    Bruno Claudio de Araujo Franco,     ------- "
PRINT*, "  -------        Erbert Barbosa Leite            ------- "
PRINT*, "  ------------------------------------------------------ "
PRINT*, "  -- Escolha o capitulo que voce quer avaliar: 4 ou 5 -- "
PRINT*, "  ------------------------------------------------------ "
READ*, cap
IF (cap==4) THEN
PRINT*, "  ------------------------------------------------------ "
PRINT*, "  --   Escolha o exemplo a ser avaliado: 1, 2 ou 3.   -- "
PRINT*, "  ------------------------------------------------------ "
READ*, ex
ELSE 
PRINT*, "  ------------------------------------------------------ "
PRINT*, "  --          D.D.(1), UPWIND(2), HIBRIDO(3)          -- "
PRINT*, "  ------------------------------------------------------ "
READ*, mt
ENDIF

IF (cap==4) THEN
PRINT*, "  ------------------------------------------------------ "
PRINT*, "  --  Digite o numero de volumes para discretizacao   -- "
PRINT*, "  ------------------------------------------------------ "
READ*, n
ELSE  IF (mt==1) THEN
PRINT*, "  ------------------------------------------------------ "
PRINT*, "  --       Escolha o caso a ser avaliado: 1, 2        -- "
PRINT*, "  ------------------------------------------------------ "
READ*, cas
PRINT*, "  ------------------------------------------------------ "
PRINT*, "  --  Digite o numero de volumes para discretizacao   -- "
PRINT*, "  ------------------------------------------------------ "
READ*, n
ELSE IF (mt==2) THEN
PRINT*, "  ------------------------------------------------------ "
PRINT*, "  --       Escolha o caso a ser avaliado: 1, 2        -- "
PRINT*, "  ------------------------------------------------------ "
READ*, casup
PRINT*, "  ------------------------------------------------------ "
PRINT*, "  --  Digite o numero de volumes para discretizacao   -- "
PRINT*, "  ------------------------------------------------------ "
READ*, n
ELSE IF (mt==3) THEN
PRINT*, "  ------------------------------------------------------ "
PRINT*, "  --       Escolha o caso a ser avaliado: 1, 2        -- "
PRINT*, "  ------------------------------------------------------ "
READ*, cash
PRINT*, "  ------------------------------------------------------ "
PRINT*, "  --  Digite o numero de volumes para discretizacao   -- "
PRINT*, "  ------------------------------------------------------ "
READ*, n
ENDIF


	!CHAMAR A SUBROTINA DO EXEMPLO ESCOLHIDO
IF (cap==4) THEN
	IF (ex == 1) THEN
		CALL EXEMPLO1()
	ELSE IF (ex == 2) THEN 
		CALL EXEMPLO2()	
	ELSE IF (ex == 3) THEN 
		CALL EXEMPLO3()
	END IF
ELSE IF (cap == 5) THEN
	IF (mt == 1) THEN
     		IF (cas == 1) THEN
			CALL EXEMPLOCC1()
		ELSE 
			CALL EXEMPLOCC2()
		END IF
	ELSE IF (mt == 2) THEN
		IF (casup == 1) THEN
			CALL EXEMPLOCCU1()
		ELSE
			CALL EXEMPLOCCU2()
		END IF
	ELSE IF (mt == 3) THEN
		IF (cash == 1) THEN
			CALL EXEMPLOCCH1()
		ELSE
			CALL EXEMPLOCCH2()
		END IF
END IF
END IF
	!CONSTANTES PARA RESOLVER POR TDMA
ALLOCATE (S(n))
DO i= 1, n, 1
S(1)=ae1/(ap1-aw1)
S(i)=ae2/(ap2-aw2*S(i-1))
S(n)=ae5/(ap5-aw5*S(i-1))
ENDDO

ALLOCATE (U(n))
DO i= 1, n, 1
U(1)=su1/ap1
U(i)=(aw2*U(i-1)+su2)/(ap2-aw2*S(i-1))
U(n)=(aw5*U(n-1)+su5)/(ap5-aw5*S(n-1))
ENDDO
	!SOLUÇAO POR SUBISTITUIÇÃO 
ALLOCATE (V(n))
DO i= n, 1,-1
V(i)=U(i)
V(i)=S(i)*V(i+1)+U(i)
ENDDO
	!ARQUIVO DE SAIDA DOS RESULTADOS
OPEN(903,FILE='Resultados.txt')

WRITE(903,*) "  -------------------- "
WRITE(903,*) "  --   Resultados   -- "
WRITE(903,*) "  -------------------- "
WRITE(903,*) "  --     Valor      -- "
WRITE(903,*) "  -------------------- "
DO i= 1,n
WRITE(903,111) V(i)
ENDDO
WRITE(903,*) "  -------------------- "

111 FORMAT ("   -",2x,F10.4,"      -")
	!BOUNDEDNESS
Bound1=(aw1+ae1)/ap1
Bound2=(aw2+ae2)/ap2
Bound5=(aw5+ae5)/ap5

Pe=F/D

OPEN(904,FILE='Criterios.txt')

WRITE(904,*) "  --------------------------------- "
WRITE(904,*) "  --         Boundedness         -- "
WRITE(904,*) "  --------------------------------- "
WRITE(904,110) "∑|aw1+ae1|/|ap'1| =",abs(Bound1)
WRITE(904,110) "∑|aw2+ae2|/|ap'2| =",abs(Bound2)
WRITE(904,110) "∑|aw5+ae5|/|ap'5| =",abs(Bound5)
WRITE(904,*) "  --------------------------------- "
WRITE(904,*) "  --       Transportiveness      -- "
WRITE(904,*) "  --------------------------------- "
WRITE(904,112) "Número de Peclet  =",Pe
WRITE(904,*) "  --------------------------------- "

110 FORMAT ("   -",2x,A21,F8.4,'  -')
112 FORMAT ("   -",2x,A20,F8.4,'  -')

IF (cap==4) THEN
	IF (ex == 1) THEN
		CALL GEXEMPLO1()
	ELSE IF (ex == 2) THEN 
		CALL GEXEMPLO2()	
	ELSE IF (ex == 3) THEN 
		CALL GEXEMPLO3()
	END IF
ELSE IF (cap == 5) THEN
	IF (mt == 1) THEN
     		IF (cas == 1) THEN
			CALL GEXEMPLOCC1()
		ELSE 
			CALL EXEMPLOCC2()
		END IF
	ELSE IF (mt == 2) THEN
		IF (casup == 1) THEN
			CALL GEXEMPLOCCU1()
		ELSE
			CALL GEXEMPLOCCU2()
		END IF
	ELSE IF (mt == 3) THEN
		IF (cash == 1) THEN
			CALL GEXEMPLOCCH1()
		ELSE
			CALL GEXEMPLOCCH2()
		END IF
END IF
END IF
	!FECHAR ARQUIVOS
	CALL SYSTEM('gnuplot SCRIPT.gnu')
	CLOSE(900)
	CLOSE(901)
	CLOSE(902)
	CLOSE(903)
	CLOSE(904)
	CALL SYSTEM('del TEMP100.dat')
	CALL SYSTEM('del TEMP101.dat')
	CALL SYSTEM('del SCRIPT.gnu')
	!SUBROTINAS
CONTAINS
!----------------------------------------------
	!EXEMPLO 1 
SUBROUTINE EXEMPLO1()
	!VALORES
k=1000
l=0.5
r=0.01
TA=100
TB=500
	!DECLARAÇÕES PARA DETERMIAÇÃO DA MATRIZ
x=l/n

aw1=0
aw2=(k*r)/x
aw5=aw2

ae1=aw2
ae2=aw2
ae5=0

ap1=3*ae1
ap2=2*aw2
ap5=3*ae1

su1=(2*k*r*TA)/x
su2=0
su5=(2*k*r*TB)/x
END SUBROUTINE EXEMPLO1
!----------------------------------------------
	!GEXEMPLO 1 
SUBROUTINE GEXEMPLO1()

OPEN(900,FILE='SCRIPT.gnu')
OPEN(901,FILE='TEMP100.dat')
OPEN(902,FILE='TEMP101.dat')

WRITE(901,*) 0, TA ,TA
DO i= 1,n,1
WRITE(901,*) (x*(i-1)+x/2), V(i)
ENDDO
WRITE(901,*) l, TB , TB

DO G= 0,1,1
WRITE(902,*) G, (800*G+100)
ENDDO

WRITE(900,*) '### GNUPLOT'
WRITE(900,*) 'set terminal pdf'					! FORMATO SAIDA
WRITE(900,*) 'set output "Grafico.pdf"'				! NOME SAIDA
WRITE(900,*) 'set title "Grafico"'	 			! TITULO

WRITE(900,*) '### EIXOS'
WRITE(900,*) 'set xlabel "Distancia" ; set ylabel "Temperatura"'
WRITE(900,*) 'set xrange [0:0.5] ; set yrange [0:500]'
WRITE(900,*) 'set xtics 0.05; set ytics 50'
WRITE(900,*) 'set xzeroaxis'
WRITE(900,*) 'set yzeroaxis'
WRITE(900,*) 'set grid'

WRITE(900,*) '### LEGENDA'
WRITE(900,*) 'set key title "legenda"'				! TITULO LEGENDA
WRITE(900,*) 'set key box'					! CAIXA LEGENDA
WRITE(900,*) 'set key right'					! POSIÇÃO HORIZONTAL
WRITE(900,*) 'set key top'					! POSIÇÃO VERTICAL
WRITE(900,*) 'set key outside'					! POSIÇÃO

WRITE(900,*) '### FONTE'
WRITE(900,*) 'set xlabel font "arial,12"'
WRITE(900,*) 'set ylabel font "arial,12"'
WRITE(900,*) 'set xtics font "arial,12"'
WRITE(900,*) 'set ytics font "arial,12"'
WRITE(900,*) 'set title font "arial,20"'
WRITE(900,*) 'set xtics font "arial,12"'

WRITE(900,*) '### PLOTAGEM'
WRITE(900,*) 'plot "TEMP100.dat" using 1:2 w p lc 1 lw 3 title "Numerica"'
WRITE(900,*) 'plot "TEMP101.dat" using 1:2 w l lc 3 lw 3 title "Analitica"'
END SUBROUTINE GEXEMPLO1
!----------------------------------------------
	!EXEMPLO 2 
SUBROUTINE EXEMPLO2()
	!VALORES
k=0.5
l=0.02
r=1
q=1000
TA=100
TB=200
	!PRECISÃO DO GRÁFICO
oi = 0
op = 0.0001
	!DECLARAÇÕES PARA DETERMIAÇÃO DA MATRIZ
x=l/n
q=q*1000

aw1=0
aw2=(k*r)/x
aw5=aw2

ae1=aw2
ae2=aw2
ae5=0

ap1=3*ae1
ap2=2*aw2
ap5=3*ae1

su2=q*r*x
su1=((2*k*r*TA)/x)+su2
su5=((2*k*r*TB)/x)+su2
END SUBROUTINE EXEMPLO2
!----------------------------------------------
	!GEXEMPLO 2 
SUBROUTINE GEXEMPLO2()
OPEN(900,FILE='SCRIPT.gnu')
OPEN(901,FILE='TEMP100.dat')
OPEN(902,FILE='TEMP101.dat')

WRITE(901,*) 0, TA
DO i= 1,n,1
WRITE(901,*) (x*(i-1)+x/2), V(i)
ENDDO
WRITE(901,*) l, TB

DO c = 1,200,1
WRITE(902,*) oi + op*c, ((((TB-TA)/l)+((q*l)/(2*k)))*(oi + op*c) - (q/(2*k))*(oi + op*c)*(oi + op*c) + TA)
ENDDO

WRITE(900,*) '### GNUPLOT'
WRITE(900,*) 'set terminal pdf'					! FORMATO SAIDA
WRITE(900,*) 'set output "Grafico.pdf"'				! NOME SAIDA
WRITE(900,*) 'set title "Grafico"'				! TITULO

WRITE(900,*) '### EIXOS'
WRITE(900,*) 'set xlabel "Distancia" ; set ylabel "Temperatura"'
WRITE(900,*) 'set xrange [0:0.02] ; set yrange [100:260]'
WRITE(900,*) 'set xtics 0.005; set ytics 50'
WRITE(900,*) 'set xzeroaxis'
WRITE(900,*) 'set yzeroaxis'
WRITE(900,*) 'set grid'

WRITE(900,*) '### LEGENDA'
WRITE(900,*) 'set key title "legenda"'				! TITULO LEGENDA
WRITE(900,*) 'set key box'					! CAIXA LEGENDA
WRITE(900,*) 'set key right'					! POSIÇÃO HORIZONTAL
WRITE(900,*) 'set key top'					! POSIÇÃO VERTICAL
WRITE(900,*) 'set key outside'					! POSIÇÃO

WRITE(900,*) '### FONTE'
WRITE(900,*) 'set xlabel font "arial,12"'
WRITE(900,*) 'set ylabel font "arial,12"'
WRITE(900,*) 'set xtics font "arial,12"'
WRITE(900,*) 'set ytics font "arial,12"'
WRITE(900,*) 'set title font "arial,20"'
WRITE(900,*) 'set xtics font "arial,12"'

WRITE(900,*) '### PLOTAGEM'
WRITE(900,*) 'plot "TEMP100.dat" using 1:2 w p lc 1 lw 3 title "Numerica" &
              & , "TEMP101.dat" using 1:2 w l lc 3 lw 3 title "Analitica"'

END SUBROUTINE GEXEMPLO2
!----------------------------------------------
	!EXEMPLO 3
SUBROUTINE EXEMPLO3()
	!VALORES
TA=100
TB=20
l=1
	!PRECISÃO DO GRÁFICO
oi=0
op=0.0001
	!DECLARAÇÕES PARA DETERMIAÇÃO DA MATRIZ
x=l/n
q=q*1000
z=25

aw1=0
aw2=1/x
aw5=aw2

ae1=aw2
ae2=aw2
ae5=0


ap1=3*aw2+z*x
ap2=2*aw2+z*x
ap5=aw2+z*x

su2=z*x*TB
su1=su2+2*aw2*TA
su5=su2
END SUBROUTINE EXEMPLO3
!----------------------------------------------
	!GEXEMPLO 3
SUBROUTINE GEXEMPLO3()
OPEN(900,FILE='SCRIPT.gnu')
OPEN(901,FILE='TEMP100.dat')
OPEN(902,FILE='TEMP101.dat')

WRITE(901,*) 0, TA
DO i= 1,n,1
WRITE(901,*) (x*(i-1)+x/2), V(i)
ENDDO
WRITE(901,*) l, TB

oi=0
op=0.0001
DO c = 1,10000,1
WRITE(902,*) oi + op*c, ((cosh(sqrt(z)*(l-(oi + op*c)))*(TA-TB))/cosh(sqrt(z)*l)+TB)
ENDDO

WRITE(900,*) '### GNUPLOT'
WRITE(900,*) 'set terminal pdf'					! FORMATO SAIDA
WRITE(900,*) 'set output "Grafico.pdf"'				! NOME SAIDA
WRITE(900,*) 'set title "Grafico"'				! TITULO

WRITE(900,*) '### EIXOS'
WRITE(900,*) 'set xlabel "Distancia" ; set ylabel "Temperatura"'
WRITE(900,*) 'set xrange [0:1] ; set yrange [20:100]'
WRITE(900,*) 'set xtics 0.1; set ytics 20'
WRITE(900,*) 'set xzeroaxis'


WRITE(900,*) 'set yzeroaxis'
WRITE(900,*) 'set grid'

WRITE(900,*) '### LEGENDA'
WRITE(900,*) 'set key title "legenda"'				! TITULO LEGENDA
WRITE(900,*) 'set key box'					! CAIXA LEGENDA
WRITE(900,*) 'set key right'					! POSIÇÃO HORIZONTAL
WRITE(900,*) 'set key top'					! POSIÇÃO VERTICAL
WRITE(900,*) 'set key outside'					! POSIÇÃO

WRITE(900,*) '### FONTE'
WRITE(900,*) 'set xlabel font "arial,12"'
WRITE(900,*) 'set ylabel font "arial,12"'
WRITE(900,*) 'set xtics font "arial,12"'
WRITE(900,*) 'set ytics font "arial,12"'
WRITE(900,*) 'set title font "arial,20"'
WRITE(900,*) 'set xtics font "arial,12"'

WRITE(900,*) '### PLOTAGEM'
WRITE(900,*) 'plot "TEMP100.dat" using 1:2 w p lc 1 lw 3 title "Numerica" & 
                   & , "TEMP101.dat" using 1:2 w l lc 3 lw 3 title "Analitica"'

END SUBROUTINE GEXEMPLO3
!----------------------------------------------
!EXEMPLO CONVECCAO CASO1 
SUBROUTINE EXEMPLOCC1()
	!VALORES
l=1
q=0.1
x=l/n
rho=1 
ux=0.1
D=q/x
F=rho*ux
phia=1
phib=0	
	!DECLARAÇÕES PARA DETERMIAÇÃO DA MATRIZ

aw1=0
aw2=D+F/2
aw5=D+F/2

ae1=D-F/2
ae2=D-F/2
ae5=0

sp1=-(2*D+F)
sp2=0
sp5=-(2*D-F)

ap1=aw1+ae1-sp1
ap2=aw2+ae2-sp2
ap5=aw5+ae5-sp5

su1=(2*D+F)*phia
su2=0
su5=(2*D-F)*phib
END SUBROUTINE EXEMPLOCC1
!----------------------------------------------
	!GEXEMPLOCONVEC CASO1
SUBROUTINE GEXEMPLOCC1()
!PRECISÃO DO GRÁFICO
oi = 0
op = 0.01

OPEN(900,FILE='SCRIPT.gnu')
OPEN(901,FILE='TEMP100.dat')
OPEN(902,FILE='TEMP101.dat')

WRITE(901,*) 0, phia
DO i= 1,n,1
WRITE(901,*) (x*(i-1)+x/2), V(i)
ENDDO
WRITE(901,*) l, phib

DO c = 1,100,1
WRITE(902,*) oi + op*c, (2.7183-exp(oi + op*c))/(1.7183)

ENDDO

WRITE(900,*) '### GNUPLOT'
WRITE(900,*) 'set terminal pdf'					! FORMATO SAIDA
WRITE(900,*) 'set output "Grafico.pdf"'				! NOME SAIDA
WRITE(900,*) 'set title "Grafico"'				! TITULO

WRITE(900,*) '### EIXOS'
WRITE(900,*) 'set xlabel "Distancia" ; set ylabel "Propiedade"'
WRITE(900,*) 'set xrange [0:1] ; set yrange [0:1]'
WRITE(900,*) 'set xtics 0.2; set ytics 0.2'	
WRITE(900,*) 'set xzeroaxis'
WRITE(900,*) 'set yzeroaxis'
WRITE(900,*) 'set grid'

WRITE(900,*) '### LEGENDA'
WRITE(900,*) 'set key title "legenda"'				! TITULO LEGENDA
WRITE(900,*) 'set key box'					! CAIXA LEGENDA
WRITE(900,*) 'set key right'					! POSIÇÃO HORIZONTAL
WRITE(900,*) 'set key top'					! POSIÇÃO VERTICAL
WRITE(900,*) 'set key outside'					! POSIÇÃO

WRITE(900,*) '### FONTE'
WRITE(900,*) 'set xlabel font "arial,12"'
WRITE(900,*) 'set ylabel font "arial,12"'
WRITE(900,*) 'set xtics font "arial,12"'
WRITE(900,*) 'set ytics font "arial,12"'
WRITE(900,*) 'set title font "arial,20"'
WRITE(900,*) 'set xtics font "arial,12"'

WRITE(900,*) '### PLOTAGEM'
WRITE(900,*) 'plot "TEMP100.dat" using 1:2 w p lc 1 lw 3 title "Numerica" &
              & , "TEMP101.dat" using 1:2 w l lc 3 lw 3 title "Analitica"'

END SUBROUTINE GEXEMPLOCC1
!----------------------------------------------
!EXEMPLOCONVEC CASO2 
SUBROUTINE EXEMPLOCC2()
	!VALORES
l=1
q=0.1
x=l/n
rho=1 
ux=2.5
D=q/x
F=rho*ux
phia=1
phib=0	
	!DECLARAÇÕES PARA DETERMIAÇÃO DA MATRIZ

aw1=0
aw2=D+F/2
aw5=D+F/2

ae1=D-F/2
ae2=D-F/2
ae5=0

sp1=-(2*D+F)
sp2=0
sp5=-(2*D-F)

ap1=aw1+ae1-sp1
ap2=aw2+ae2-sp2
ap5=aw5+ae5-sp5

su1=(2*D+F)*phia
su2=0
su5=(2*D-F)*phib
END SUBROUTINE EXEMPLOCC2
!----------------------------------------------
	!GEXEMPLOCONVEC CASO2
SUBROUTINE GEXEMPLOCC2()
!PRECISÃO DO GRÁFICO
oi = 0
op = 0.01

OPEN(900,FILE='SCRIPT.gnu')
OPEN(901,FILE='TEMP100.dat')
OPEN(902,FILE='TEMP101.dat')

WRITE(901,*) 0, phia
DO i= 1,n,1
WRITE(901,*) (x*(i-1)+x/2), V(i)
ENDDO
WRITE(901,*) l, phib

DO c = 1,100,1
WRITE(902,*) oi + op*c, 1+(1-exp(25*(oi + op*c)))/(7.20*10000000000)

ENDDO

WRITE(900,*) '### GNUPLOT'
WRITE(900,*) 'set terminal pdf'					! FORMATO SAIDA
WRITE(900,*) 'set output "Grafico.pdf"'				! NOME SAIDA
WRITE(900,*) 'set title "Grafico"'				! TITULO

WRITE(900,*) '### EIXOS'
WRITE(900,*) 'set xlabel "Distancia" ; set ylabel "Propiedade"'
WRITE(900,*) 'set xrange [0:1] ; set yrange [0:2.5]'
WRITE(900,*) 'set xtics 0.2; set ytics 0.5'	
WRITE(900,*) 'set xzeroaxis'
WRITE(900,*) 'set yzeroaxis'
WRITE(900,*) 'set grid'

WRITE(900,*) '### LEGENDA'
WRITE(900,*) 'set key title "legenda"'				! TITULO LEGENDA
WRITE(900,*) 'set key box'					! CAIXA LEGENDA
WRITE(900,*) 'set key right'					! POSIÇÃO HORIZONTAL
WRITE(900,*) 'set key top'					! POSIÇÃO VERTICAL
WRITE(900,*) 'set key outside'					! POSIÇÃO

WRITE(900,*) '### FONTE'
WRITE(900,*) 'set xlabel font "arial,12"'
WRITE(900,*) 'set ylabel font "arial,12"'
WRITE(900,*) 'set xtics font "arial,12"'
WRITE(900,*) 'set ytics font "arial,12"'
WRITE(900,*) 'set title font "arial,20"'
WRITE(900,*) 'set xtics font "arial,12"'

WRITE(900,*) '### PLOTAGEM'
WRITE(900,*) 'plot "TEMP100.dat" using 1:2 w p lc 1 lw 3 title "Numerica" &
              & , "TEMP101.dat" using 1:2 w l lc 3 lw 3 title "Analitica"'

END SUBROUTINE GEXEMPLOCC2
!----------------------------------------------
!EXEMPLOCONVEC CASO1 UPWIND
SUBROUTINE EXEMPLOCCU1()
	!VALORES
l=1
q=0.1
x=l/n
rho=1 
ux=0.1
D=q/x
F=rho*ux
phia=1
phib=0	
	!DECLARAÇÕES PARA DETERMIAÇÃO DA MATRIZ

aw1=0
aw2=D+F
aw5=D+F

ae1=D
ae2=D
ae5=0

sp1=-(2*D+F)
sp2=0
sp5=-(2*D)

ap1=aw1+ae1-sp1
ap2=aw2+ae2-sp2
ap5=aw5+ae5-sp5

su1=(2*D+F)*phia
su2=0
su5=(2*D)*phib
END SUBROUTINE EXEMPLOCCU1
!----------------------------------------------
	!GEXEMPLOCONVEC CASO1 UPWIND	
SUBROUTINE GEXEMPLOCCU1()
!PRECISÃO DO GRÁFICO
oi = 0
op = 0.01

OPEN(900,FILE='SCRIPT.gnu')
OPEN(901,FILE='TEMP100.dat')
OPEN(902,FILE='TEMP101.dat')

WRITE(901,*) 0, phia
DO i= 1,n,1
WRITE(901,*) (x*(i-1)+x/2), V(i)
ENDDO
WRITE(901,*) l, phib

DO c = 1,100,1
WRITE(902,*) oi + op*c, (2.7183-exp(oi + op*c))/(1.7183)

ENDDO

WRITE(900,*) '### GNUPLOT'
WRITE(900,*) 'set terminal pdf'					! FORMATO SAIDA
WRITE(900,*) 'set output "Grafico.pdf"'				! NOME SAIDA
WRITE(900,*) 'set title "Grafico"'				! TITULO

WRITE(900,*) '### EIXOS'
WRITE(900,*) 'set xlabel "Distancia" ; set ylabel "Propiedade"'
WRITE(900,*) 'set xrange [0:1] ; set yrange [0:1]'
WRITE(900,*) 'set xtics 0.2; set ytics 0.2'	
WRITE(900,*) 'set xzeroaxis'
WRITE(900,*) 'set yzeroaxis'
WRITE(900,*) 'set grid'

WRITE(900,*) '### LEGENDA'
WRITE(900,*) 'set key title "legenda"'				! TITULO LEGENDA
WRITE(900,*) 'set key box'					! CAIXA LEGENDA
WRITE(900,*) 'set key right'					! POSIÇÃO HORIZONTAL
WRITE(900,*) 'set key top'					! POSIÇÃO VERTICAL
WRITE(900,*) 'set key outside'					! POSIÇÃO

WRITE(900,*) '### FONTE'
WRITE(900,*) 'set xlabel font "arial,12"'
WRITE(900,*) 'set ylabel font "arial,12"'
WRITE(900,*) 'set xtics font "arial,12"'
WRITE(900,*) 'set ytics font "arial,12"'
WRITE(900,*) 'set title font "arial,20"'
WRITE(900,*) 'set xtics font "arial,12"'

WRITE(900,*) '### PLOTAGEM'
WRITE(900,*) 'plot "TEMP100.dat" using 1:2 w p lc 1 lw 3 title "Numerica" &
              & , "TEMP101.dat" using 1:2 w l lc 3 lw 3 title "Analitica"'

END SUBROUTINE GEXEMPLOCCU1
!----------------------------------------------
!EXEMPLOCONVEC CASO2 UPWIND
SUBROUTINE EXEMPLOCCU2()
	!VALORES
l=1
q=0.1
x=l/n
rho=1 
ux=2.5
D=q/x
F=rho*ux
phia=1
phib=0	
	!DECLARAÇÕES PARA DETERMIAÇÃO DA MATRIZ

aw1=0
aw2=D+F
aw5=D+F

ae1=D
ae2=D
ae5=0

sp1=-(2*D+F)
sp2=0
sp5=-(2*D)

ap1=aw1+ae1-sp1
ap2=aw2+ae2-sp2
ap5=aw5+ae5-sp5

su1=(2*D+F)*phia
su2=0
su5=(2*D)*phib
END SUBROUTINE EXEMPLOCCU2
!----------------------------------------------
	!GEXEMPLOCONVEC CASO2 UPWIND
SUBROUTINE GEXEMPLOCCU2
!PRECISÃO DO GRÁFICO
oi = 0
op = 0.01

OPEN(900,FILE='SCRIPT.gnu')
OPEN(901,FILE='TEMP100.dat')
OPEN(902,FILE='TEMP101.dat')

WRITE(901,*) 0, phia
DO i= 1,n,1
WRITE(901,*) (x*(i-1)+x/2), V(i)
ENDDO
WRITE(901,*) l, phib

DO c = 1,100,1
WRITE(902,*) oi + op*c, 1+(1-exp(25*(oi + op*c)))/(7.20*10000000000)

ENDDO

WRITE(900,*) '### GNUPLOT'
WRITE(900,*) 'set terminal pdf'					! FORMATO SAIDA
WRITE(900,*) 'set output "Grafico.pdf"'				! NOME SAIDA
WRITE(900,*) 'set title "Grafico"'				! TITULO

WRITE(900,*) '### EIXOS'
WRITE(900,*) 'set xlabel "Distancia" ; set ylabel "Propiedade"'
WRITE(900,*) 'set xrange [0:1] ; set yrange [0:2.5]'
WRITE(900,*) 'set xtics 0.2; set ytics 0.5'	
WRITE(900,*) 'set xzeroaxis'
WRITE(900,*) 'set yzeroaxis'
WRITE(900,*) 'set grid'

WRITE(900,*) '### LEGENDA'
WRITE(900,*) 'set key title "legenda"'				! TITULO LEGENDA
WRITE(900,*) 'set key box'					! CAIXA LEGENDA
WRITE(900,*) 'set key right'					! POSIÇÃO HORIZONTAL
WRITE(900,*) 'set key top'					! POSIÇÃO VERTICAL
WRITE(900,*) 'set key outside'					! POSIÇÃO

WRITE(900,*) '### FONTE'
WRITE(900,*) 'set xlabel font "arial,12"'
WRITE(900,*) 'set ylabel font "arial,12"'
WRITE(900,*) 'set xtics font "arial,12"'
WRITE(900,*) 'set ytics font "arial,12"'
WRITE(900,*) 'set title font "arial,20"'
WRITE(900,*) 'set xtics font "arial,12"'

WRITE(900,*) '### PLOTAGEM'
WRITE(900,*) 'plot "TEMP100.dat" using 1:2 w p lc 1 lw 3 title "Numerica" &
              & , "TEMP101.dat" using 1:2 w l lc 3 lw 3 title "Analitica"'

END SUBROUTINE GEXEMPLOCCU2
!----------------------------------------------
!EXEMPLOCONVEC CASO1 HIBRIDO
SUBROUTINE EXEMPLOCCH1()
	!VALORES
l=1
q=0.1
x=l/n
rho=1 
ux=0.1
D=q/x
F=rho*ux
phia=1
phib=0	
Pe=F/D
IF (Pe < 2) THEN
	!DECLARAÇÕES PARA DETERMIAÇÃO DA MATRIZ Pe < 2 
aw1=0
aw2=D+F/2
aw5=D+F/2

ae1=D-F/2
ae2=D-F/2
ae5=0

sp1=-(2*D+F)
sp2=0
sp5=-(2*D-F)

ap1=aw1+ae1-sp1
ap2=aw2+ae2-sp2
ap5=aw5+ae5-sp5

su1=(2*D+F)*phia
su2=0
su5=(2*D-F)*phib
ELSE
	!DECLARAÇÕES PARA DETERMIAÇÃO DA MATRIZ Pe >= 2 

aw1=0
aw2=F
aw5=F

ae1=0
ae2=0
ae5=0

sp1=-(2*D+F)
sp2=0
sp5=-(2*D)

ap1=aw1+ae1-sp1
ap2=aw2+ae2-sp2
ap5=aw5+ae5-sp5

su1=(2*D+F)*phia
su2=0
su5=(2*D)*phib
ENDIF
END SUBROUTINE EXEMPLOCCH1
!----------------------------------------------
	!GEXEMPLOCONVEC CASO1 HIBRIDO	
SUBROUTINE GEXEMPLOCCH1()
!PRECISÃO DO GRÁFICO
oi = 0
op = 0.01

OPEN(900,FILE='SCRIPT.gnu')
OPEN(901,FILE='TEMP100.dat')
OPEN(902,FILE='TEMP101.dat')

WRITE(901,*) 0, phia
DO i= 1,n,1
WRITE(901,*) (x*(i-1)+x/2), V(i)
ENDDO
WRITE(901,*) l, phib

DO c = 1,100,1
WRITE(902,*) oi + op*c, (2.7183-exp(oi + op*c))/(1.7183)

ENDDO

WRITE(900,*) '### GNUPLOT'
WRITE(900,*) 'set terminal pdf'					! FORMATO SAIDA
WRITE(900,*) 'set output "Grafico.pdf"'				! NOME SAIDA
WRITE(900,*) 'set title "Grafico"'				! TITULO

WRITE(900,*) '### EIXOS'
WRITE(900,*) 'set xlabel "Distancia" ; set ylabel "Propiedade"'
WRITE(900,*) 'set xrange [0:1] ; set yrange [0:1]'
WRITE(900,*) 'set xtics 0.2; set ytics 0.2'	
WRITE(900,*) 'set xzeroaxis'
WRITE(900,*) 'set yzeroaxis'
WRITE(900,*) 'set grid'

WRITE(900,*) '### LEGENDA'
WRITE(900,*) 'set key title "legenda"'				! TITULO LEGENDA
WRITE(900,*) 'set key box'					! CAIXA LEGENDA
WRITE(900,*) 'set key right'					! POSIÇÃO HORIZONTAL
WRITE(900,*) 'set key top'					! POSIÇÃO VERTICAL
WRITE(900,*) 'set key outside'					! POSIÇÃO

WRITE(900,*) '### FONTE'
WRITE(900,*) 'set xlabel font "arial,12"'
WRITE(900,*) 'set ylabel font "arial,12"'
WRITE(900,*) 'set xtics font "arial,12"'
WRITE(900,*) 'set ytics font "arial,12"'
WRITE(900,*) 'set title font "arial,20"'
WRITE(900,*) 'set xtics font "arial,12"'

WRITE(900,*) '### PLOTAGEM'
WRITE(900,*) 'plot "TEMP100.dat" using 1:2 w p lc 1 lw 3 title "Numerica" &
              & , "TEMP101.dat" using 1:2 w l lc 3 lw 3 title "Analitica"'

END SUBROUTINE GEXEMPLOCCH1
!----------------------------------------------
!EXEMPLOCONVEC CASO2 HIBRIDO
SUBROUTINE EXEMPLOCCH2()
	!VALORES
l=1
q=0.1
x=l/n
rho=1 
ux=2.5
D=q/x
F=rho*ux
phia=1
phib=0	
Pe=F/D
IF (Pe < 2) THEN
	!DECLARAÇÕES PARA DETERMIAÇÃO DA MATRIZ Pe < 2 
aw1=0
aw2=D+F/2
aw5=D+F/2

ae1=D-F/2
ae2=D-F/2
ae5=0

sp1=-(2*D+F)
sp2=0
sp5=-(2*D-F)

ap1=aw1+ae1-sp1
ap2=aw2+ae2-sp2
ap5=aw5+ae5-sp5

su1=(2*D+F)*phia
su2=0
su5=(2*D-F)*phib
ELSE
	!DECLARAÇÕES PARA DETERMIAÇÃO DA MATRIZ Pe >= 2 

aw1=0
aw2=F
aw5=F

ae1=0
ae2=0
ae5=0

sp1=-(2*D+F)
sp2=0
sp5=-(2*D)

ap1=aw1+ae1-sp1
ap2=aw2+ae2-sp2
ap5=aw5+ae5-sp5

su1=(2*D+F)*phia
su2=0
su5=(2*D)*phib
ENDIF
END SUBROUTINE EXEMPLOCCH2
!----------------------------------------------
	!GEXEMPLOCONVEC CASO2 HIBRIDO	
SUBROUTINE GEXEMPLOCCH2()



!PRECISÃO DO GRÁFICO
oi = 0

op = 0.01

OPEN(900,FILE='SCRIPT.gnu')
OPEN(901,FILE='TEMP100.dat')

OPEN(902,FILE='TEMP101.dat')

WRITE(901,*) 0, phia

DO i= 1,n,1
WRITE(901,*) (x*(i-1)+x/2), V(i)
ENDDO

WRITE(901,*) l, phib

DO c = 1,100,1

WRITE(902,*) oi + op*c, 1+(1-exp(25*(oi + op*c)))/(7.20*10000000000)

ENDDO


WRITE(900,*) '### GNUPLOT'
WRITE(900,*) 'set terminal pdf'					! FORMATO SAIDA
WRITE(900,*) 'set output "Grafico.pdf"'				! NOME SAIDA

WRITE(900,*) 'set title "Grafico"'				! TITULO

WRITE(900,*) '### EIXOS'

WRITE(900,*) 'set xlabel "Distancia" ; set ylabel "Propiedade"'
WRITE(900,*) 'set xrange [0:1] ; set yrange [0:2.5]'
WRITE(900,*) 'set xtics 0.2; set ytics 0.5'	
WRITE(900,*) 'set xzeroaxis'
WRITE(900,*) 'set yzeroaxis'
WRITE(900,*) 'set grid'

WRITE(900,*) '### LEGENDA'
WRITE(900,*) 'set key title "legenda"'				! TITULO LEGENDA
WRITE(900,*) 'set key box'					! CAIXA LEGENDA
WRITE(900,*) 'set key right'					! POSIÇÃO HORIZONTAL
WRITE(900,*) 'set key top'					! POSIÇÃO VERTICAL
WRITE(900,*) 'set key outside'					! POSIÇÃO

WRITE(900,*) '### FONTE'
WRITE(900,*) 'set xlabel font "arial,12"'
WRITE(900,*) 'set ylabel font "arial,12"'
WRITE(900,*) 'set xtics font "arial,12"'
WRITE(900,*) 'set ytics font "arial,12"'
WRITE(900,*) 'set title font "arial,20"'
WRITE(900,*) 'set xtics font "arial,12"'

WRITE(900,*) '### PLOTAGEM'
WRITE(900,*) 'plot "TEMP100.dat" using 1:2 w p lc 1 lw 3 title "Numerica" &
              & , "TEMP101.dat" using 1:2 w l lc 3 lw 3 title "Analitica"'


END SUBROUTINE GEXEMPLOCCH2
!----------------------------------------------
END
