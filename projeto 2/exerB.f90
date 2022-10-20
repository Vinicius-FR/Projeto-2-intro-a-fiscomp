!Definindo a função f(x) dada:
real(8) function f(x0, m, h)
  real(8), intent(in) :: x0, h
  integer, intent(in) :: m
  real(8) :: x
  
  x = x0 + m * h
  f = sin(x / 2d0)*cos(x)

end function f

program exerB
  implicit none

!Definindo as variáveis auxiliares:
  integer :: n, i, j, k
  integer, allocatable, dimension(:) :: N_list
  real(8), allocatable, dimension(:) :: h, desvio_trapezio, desvio_simpson, desvio_bode
  real(8) :: f, int_exata, int_trapezio, int_simpson, int_bode

!Integral calculada de forma direta/exata:
  int_exata = (cos(0.5d0) - (1d0/3d0)*cos(1.5d0)) - (cos(0d0) - (1d0/3d0)*cos(0d0))

!Lendo o número n de N's que serão testados e seus respectivos valores de cada teste
  open(1, file = 'tabB_in.dat', status='old', action = 'read')

    read(1,*) n
    
    allocate(N_list(n), h(n), desvio_trapezio(n), desvio_simpson(n), desvio_bode(n)) !Alocando o tamanho das listas

    read(1,*) N_list

  close(1)

  do i = 1, n
    h(i) = 1.0d0/N_list(i) !Calculando o espaçamento h para cada caso
  end do

  open(2, file = 'tabB_out.dat', status='replace')
    write(2,*) '      N              ', ' h                   ', ' Regra do Trapézio        ',&
    ' Regra de Simpson           ', ' Regra de Bode ' !Cabeçalho da tabela

  do i= 1, n
    !Setando os valores iniciais das integrais para cada loop:
    int_trapezio = 0.0d0
    int_simpson = 0.0d0
    int_bode = 0.0d0

    do j = 0, N_list(i), 2
      if (j /= N_list(i)) then

        int_trapezio = int_trapezio + ((h(i)/2)*(f(0.0d0,j+2,h(i)) + 2d0*f(0.0d0,j+1,h(i)) + f(0.0d0,j,h(i)))) !Método do trapézio
        int_simpson = int_simpson + ((h(i)/3)*(f(0.0d0,j,h(i)) + 4d0*f(0.0d0,j,h(i)) + f(0.0d0,j,h(i)))) !Método de Simpson

      else
        
        !Ao chegar no limite superior da integral (j = N), calcula-se o desvio e encerra-se o loop:
        desvio_trapezio(i) = abs(int_trapezio - int_exata)
        desvio_simpson(i) = abs(int_simpson - int_exata)
        exit

      end if
    end do

    do k = 0, N_list(i), 4 !Método de Bode
      if (k /= N_list(i)) then

        int_bode = int_bode + ((2d0*h(i)/45)*(7*f(0.0d0,k,h(i)) + 32d0*f(0.0d0,k+1,h(i)) + 12d0*f(0.0d0,k+2,h(i)) &
         + 32d0*f(0.0d0,k+3,h(i)) + 7d0*f(0.0d0,k+4,h(i))))

      else

        desvio_bode = abs(int_bode - int_exata)
        exit

      end if
    end do

    !Escrevendo as linhas da tabela a cada loop no arquivo de saída:
    write(2,*) N_list(i), h(i), desvio_trapezio(i), desvio_simpson(i), desvio_bode(i)
  end do

  !Conclusão
  close(2)
  print*, 'Ao analisar a tabela de saída, nota-se, inicialmente, que o valor de N que fornece o menor erro',&
  'para o método do trapézio e o de Simpson é N=4096, o que condiz com a ideia de separar a integral contínua',& 
  'em somas discretas de espaçamento, 2h cada vez menores para obter maior precisão. Porém, para o método de Bode',&
  'os erros obtidos no intervalo 256<=N<=2048 são equivalentes/muito próximos e menores que o caso N=4096, ou seja,',&
  'é válido considerar um número mais baixo de divisões, como N=256, para obter uma melhor precisão com menos loops.',&
  'para melhorar o desempenho e, assim, diminuir o tempo de processamento ao utilizar esse método.'

end program