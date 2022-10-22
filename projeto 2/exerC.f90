!Função para a qual descobriremos as 3 raízes:
real(8) function f(x)
  real(8), intent(in) :: x

  f = x**3 - x**2 - 2d0*x + 1

end function f

!Derivada da função:
real(8) function df(x)
  real(8), intent(in) :: x

  df = 3d0*(x**2) - 2d0*x - 2d0

end function df

program exerC
  implicit none

  integer :: i, iter
  real(8) :: f, df
  real(8) :: dir1, dir2, dir3, delta !Variáveis para a busca direta
  real(8) :: nr1, nr2, nr3 !Variáveis para o método de Newton-Raphson
  real(8) :: r1_x0, r1_x1, sec1, r2_x0, r2_x1, sec2, r3_x0, r3_x1, sec3 !Variáveis para o método da secante

  !Lendo o número de iterações desejado:
  print *, "Insira o número de iterações:"
  read(*,*) iter
  
  open(1, file = 'tabC_out.dat', status='replace')
  !Cabeçalho da tabela de saída:
  write(1,*) '         iter        ', 'dir1                      ', 'dir2                      ', 'dir3                      ',&
   'NR1                       ', 'NR2                       ', 'NR3                       ',&
    'sec1                     ', 'sec2                       ', 'sec3'

  !Analisando o gráfico da função, é possível "chutar" valores iniciais razoavelmente próximos das raízes
  !Condições iniciais para a busca direta:
  dir1 = -1.0d0
  dir2 = 0.3d0
  dir3 = 1.6d0
  delta = 0.04d0

  !Condições iniciais para o método de Newton-Raphson:
  nr1 = -1.0d0
  nr2 = 0.3d0
  nr3 = 1.6d0

  !Condições iniciais para o método da secante:
  r1_x0 = -1.5d0
  r1_x1 = -1.0d0

  r2_x0 = 0.1d0
  r2_x1 = 0.5d0

  r3_x0 = 1.3d0
  r3_x1 = 1.8d0

  !Iniciando então o loop que aplicará os cálculos para cada método, levando em conta os dados atribuidos anteriormente:
  do i = 1, iter, 1
    if (f(dir1 - delta) > 0) then !Quando a função inverte de sinal, encerra-se a busca
      dir1 = dir1 - delta
    end if

    if (f(dir2 + delta) > 0) then
      dir2 = dir2 + delta
    end if

    if (f(dir3 + delta) < 0) then
      dir3 = dir3 + delta
    end if

    !Método de Newton-Raphson:
    nr1 = nr1 - (f(nr1)/df(nr1))
    nr2 = nr2 - (f(nr2)/df(nr2))
    nr3 = nr3 - (f(nr3)/df(nr3))

    !Método da secante:
    sec1 = r1_x1 - f(r1_x1)*((r1_x1 - r1_x0)/(f(r1_x1) - f(r1_x0)))
    r1_x0 = r1_x1
    r1_x1 = sec1

    sec2 = r2_x1 - f(r2_x1)*((r2_x1 - r2_x0)/(f(r2_x1) - f(r2_x0)))
    r2_x0 = r2_x1
    r2_x1 = sec2

    sec3 = r3_x1 - f(r3_x1)*((r3_x1 - r3_x0)/(f(r3_x1) - f(r3_x0)))
    r3_x0 = r3_x1
    r3_x1 = sec3

    write(1,*) i, dir1, dir2, dir3, nr1, nr2, nr3, sec1, sec2, sec3
  end do
  close(1)

end program
