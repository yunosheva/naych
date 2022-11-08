clear();

h = 1
dt = 1
n = 50
teta = 1/3
u = zeros(n)
for i = 1 : n
    u(i) = 0
end
tau = 1
T = 0.8
k = 0.01
A = -0.152

//задаем начальную плотность
rho = zeros(n)
for i = 1 : n
    rho(i) = 1 + rand()/100
end

// задаем равновесную функцию
function F = f(c,u,w,rho)
    F = w*rho*(1 + c*u/teta + (c*u)^2/(2*teta^2) - u^2/2/teta)
endfunction

// решеточное уравнение Больцмана
function F_1 = F(f_k,f_ke)
    F_1 = f_k + (f_ke - f_k)/tau
endfunction

// уравнение состояния Ван-дер-Ваальса
function pl = P(rho)
    pl = 8*rho*T/(3 - rho) - 3*rho^2
endfunction


// задаем вектор возможных скоростей частиц (c_0 = c(1) и т.д.)
c = zeros(3)
c(2) = 1
c(3) = -1

// задаем коэффициенты w_k (w_0 = w(1) и т.д.)
w = zeros(3)
w(1) = 2/3
w(2) = 1/6
w(3) = 1/6

// посчитаем начальные одночастичные функции распределения f_k 
// и равновесные функции f_ke
f_0e = zeros(n)
f_2e = zeros(n)
f_1e = zeros(n)
f_0 = zeros(n)
f_2 = zeros(n)
f_1 = zeros(n)

buf = zeros(n)

for i = 1 : n
    f_0e(i) = f(c(1), u(i), w(1), rho(i))
    f_1e(i) = f(c(2), u(i), w(2), rho(i))
    f_2e(i) = f(c(3), u(i), w(3), rho(i))
    f_0(i) = f_0e(i)
    f_1(i) = f_1e(i)
    f_2(i) = f_2e(i)
end 

buf = f_1
f_1(1) = buf(n)
for i = 1 : n-1
    f_1(i+1) = buf(i)
end

buf = f_2
f_2(n) = buf(1)
for i = 1 : n-1
    f_2(i) = buf(i+1)
end

// посчитаем новую плотность
for i = 1 : n
    rho(i) = f_0(i) + f_1(i) + f_2(i)
end

// посчитаем новую скорость вещества в узле u
for i = 1 : n
    u(i) = f_1(i)*c(2)/rho(i) + f_2(i)*c(3)/rho(i)
end

// "эффективная" плотность
Fi = zeros(n)
for i = 1 : n
    Fi(i) = sqrt(rho(i)*teta - k*P(rho(i)))
end

// посчитаем силу
FF = zeros(n)
U = zeros(n)
U(n) = 0
FF(n) = 0
for i = 2 : n-1
    FF(i) = (A*Fi(i+1) + (1 - 2*A)*Fi(i) + A*Fi(i-1))*(Fi(i+1) - Fi(i-1))
    U(i) = FF(i)/rho(i)
end

// посчитаем новые равновесные функции и функции при другой скорости
f_01e = zeros(n)
f_11e = zeros(n)
f_21e = zeros(n)
for i = 1 : n
    f_0e(i) = f(c(1), u(i), w(1), rho(i))
    f_1e(i) = f(c(2), u(i), w(2), rho(i))
    f_2e(i) = f(c(3), u(i), w(3), rho(i))
    f_01e(i) = f(c(1), u(i) + U(i), w(1), rho(i))
    f_11e(i) = f(c(2), u(i) + U(i), w(2), rho(i))
    f_21e(i) = f(c(3), u(i) + U(i), w(3), rho(i))
end 

// сделаем одну итерацию
F_0 = zeros(n)
for i = 1 : n
    F_0(i) = F(f_0(i), f_0e(i)) + f_01e(i) - f_0e(i)
end

F_1 = zeros(n)
for i = 1 : n
    F_1(i) = F(f_1(i), f_1e(i)) + f_11e(i) - f_1e(i)
end

F_2 = zeros(n)
for i = 1 : n
    F_2(i) = F(f_2(i), f_2e(i)) + f_21e(i) - f_2e(i)
end

M = 50
rho_M = zeros(M,n)
for j = 1 : M
    // задаем одночастичные функции
    f_0 = F_0
    f_1 = F_1
    f_2 = F_2
    
    // учитываем движение
    buf = f_1
    f_1(1) = buf(n)
    for i = 1 : n-1
        f_1(i+1) = buf(i)
    end

    buf = f_2
    f_2(n) = buf(1)
    for i = 1 : n-1
        f_2(i) = buf(i+1)
    end
    
    // посчитаем новую плотность
    for i = 1 : n
        rho(i) = f_0(i) + f_1(i) + f_2(i)
        rho_M(j, i) = rho(i)
    end

    // посчитаем новую скорость вещества в узле u
    for i = 1 : n
        u(i) = f_1(i)*c(2)/rho(i) + f_2(i)*c(3)/rho(i)
    end
    
    Fi = zeros(n)
    for i = 1 : n
        Fi(i) = sqrt(rho(i)*teta - k*P(rho(i)))
    end

    // посчитаем силу
    FF = zeros(n)
    U = zeros(n)
    U(n) = 0
    FF(n) = 0
    for i = 2 : n-1
        FF(i) = (A*Fi(i+1) + (1 - 2*A)*Fi(i) + A*Fi(i-1))*(Fi(i+1) - Fi(i-1))
        U(i) = FF(i)/rho(i)
    end


    // посчитаем новые равновесные функции
    for i = 1 : n
        f_0e(i) = f(c(1), u(i), w(1), rho(i))
        f_1e(i) = f(c(2), u(i), w(2), rho(i))
        f_2e(i) = f(c(3), u(i), w(3), rho(i))
        f_01e(i) = f(c(1), U(i), w(1), rho(i))
        f_11e(i) = f(c(2), U(i), w(2), rho(i))
        f_21e(i) = f(c(3), U(i), w(3), rho(i))
    end 


    // сделаем одну итерацию
    F_0 = zeros(n)
    for i = 1 : n
        F_0(i) = F(f_0(i), f_0e(i)) + f_01e(i) - f_0e(i)
    end

    F_1 = zeros(n)
    for i = 1 : n
        F_1(i) = F(f_1(i), f_1e(i)) + f_11e(i) - f_1e(i)
    end

    F_2 = zeros(n)
    for i = 1 : n
        F_2(i) = F(f_2(i), f_2e(i)) + f_21e(i) - f_2e(i)
    end
end

// рисуем плотность
x=zeros(n)
x = 0 : h : n-1
rho_1 = zeros(n)
for j= 1:n
    rho_1(j)=rho_M(1, j)
end

rho_10 = zeros(n)
for j= 1:n
    rho_10(j)=rho_M(10, j)
end

rho_20 = zeros(n)
for j= 1:n
    rho_20(j)=rho_M(20, j)
end

rho_30 = zeros(n)
for j= 1:n
    rho_30(j)=rho_M(30, j)
end

rho_40 = zeros(n)
for j= 1:n
    rho_40(j)=rho_M(40, j)
end

rho_50 = zeros(n)
for j= 1:n
    rho_50(j)=rho_M(50, j)
end

subplot(1,1,1);
    plot2d(x, rho_1, style = color('blue'));
    plot2d(x, rho_10, style = color('yellow'));
    plot2d(x, rho_20, style = color('pink'));
    plot2d(x, rho_30, style = color('brown'));
    plot2d(x, rho_40, style = color('black'));
    plot2d(x, rho_50, style = color('orange'));
    xtitle("Плотность")
description(1) = "Плотность 1";
description(2) = "Плотность 10";
description(3) = "Плотность 20";
description(4) = "Плотность 30";
description(5) = "Плотность 40";
description(6) = "Плотность 50";
legend(description, "lower_caption");

disp(log2(0.147674/0.0415537))
disp(log2(0.0415537/0.010714 ))
disp(log2(0.010714/0.00270123))
disp(log2(0.00270123/0.00067735))
disp(log2(0.00067735/0.000169259))

m = 10
q = zeros(2,n,m)
zeros(q)








