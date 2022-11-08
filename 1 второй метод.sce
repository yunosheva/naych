clear();

h = 1
dt = 1
n = 20
teta = 1/3
u = zeros(n)
for i = 1 : n
    u(i) = 0
end
tau = 1
T = 0.8

//задаем начальную плотность
rho = zeros(n)
for i = 1 : n/2
    rho(i) = 2
end
for i = n/2 : n
    rho(i) = 1
end

// задаем равновесную функцию
function F = f(c,u,w,rho)
    F = w*rho*(1 + c*u/teta + (c*u)^2/(2*teta^2) - u^2/2/teta)
endfunction

// решеточное уравнение Больцмана
function F_1 = F(f_k,f_ke)
    F_1 = f_k + (f_ke - f_k)/tau
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

// посчитаем новые равновесные функции
for i = 1 : n
    f_0e(i) = f(c(1), u(i), w(1), rho(i))
    f_1e(i) = f(c(2), u(i), w(2), rho(i))
    f_2e(i) = f(c(3), u(i), w(3), rho(i))
end 


// сделаем одну итерацию
F_0 = zeros(n)
for i = 1 : n
    F_0(i) = F(f_0(i), f_0e(i))
end

F_1 = zeros(n)
for i = 1 : n
    F_1(i) = F(f_1(i), f_1e(i))
end

F_2 = zeros(n)
for i = 1 : n
    F_2(i) = F(f_2(i), f_2e(i))
end

M = 20
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
    sum = 0
    // посчитаем новую плотность
    for i = 1 : n
        rho(i) = f_0(i) + f_1(i) + f_2(i)
        sum = sum + rho(i)
        rho_M(j, i) = rho(i)
    end
    
    // посчитаем новую скорость вещества в узле u
    ux = 0
    for i = 1 : n
        u(i) = f_1(i)*c(2)/rho(i) + f_2(i)*c(3)/rho(i)
        ux = ux + u(i)
    end
    disp(ux)
    // посчитаем новые равновесные функции
    for i = 1 : n
        f_0e(i) = f(c(1), u(i), w(1), rho(i))
        f_1e(i) = f(c(2), u(i), w(2), rho(i))
        f_2e(i) = f(c(3), u(i), w(3), rho(i))
    end 


    // сделаем одну итерацию
    F_0 = zeros(n)
    for i = 1 : n
        F_0(i) = F(f_0(i), f_0e(i))
    end

    F_1 = zeros(n)
    for i = 1 : n
        F_1(i) = F(f_1(i), f_1e(i))
    end

    F_2 = zeros(n)
    for i = 1 : n
        F_2(i) = F(f_2(i), f_2e(i))
    end
end

// рисуем плотность
x=zeros(n)
x = 0 : h : 19
rho_1 = zeros(n)
for j= 1:n
    rho_1(j)=rho_M(1, j)
end

rho_4 = zeros(n)
for j= 1:n
    rho_4(j)=rho_M(4, j)
end

rho_7 = zeros(n)
for j= 1:n
    rho_7(j)=rho_M(7, j)
end

rho_10 = zeros(n)
for j= 1:n
    rho_10(j)=rho_M(10, j)
end

rho_15 = zeros(n)
for j= 1:n
    rho_15(j)=rho_M(15, j)
end

rho_20 = zeros(n)
for j= 1:n
    rho_20(j)=rho_M(20, j)
end

subplot(1,1,1);
    plot2d(x, rho_1, style = color('blue'));
    plot2d(x, rho_4, style = color('yellow'));
    plot2d(x, rho_7, style = color('pink'));
    plot2d(x, rho_10, style = color('brown'));
    plot2d(x, rho_15, style = color('black'));
    plot2d(x, rho_20, style = color('orange'));
    xtitle("Плотность")
description(1) = "Плотность 1";
description(2) = "Плотность 4";
description(3) = "Плотность 7";
description(4) = "Плотность 10";
description(5) = "Плотность 15";
description(6) = "Плотность 20";
legend(description, "lower_caption");
















