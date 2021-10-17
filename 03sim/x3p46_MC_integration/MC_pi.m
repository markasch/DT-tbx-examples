% MC integration
N = 1000000;
count = 0;

for i = 1:N
    x = rand(); y=rand();
    if (x^2 + y^2) < 1
        count = count +1;
    end
end

disp(['Approximation of pi = ', num2str(4*count/N)]);
