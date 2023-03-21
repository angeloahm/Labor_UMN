function root = bisection(f, a, b, max_iter, tol)

for iter=1:max_iter
    mean = (a+b)/2;
    if f(a)*f(b)>0
        root=mean;
        fprintf('Error:  No change in sign \n')
        break
    else
        if f(mean)*f(a)>0
            a = mean;
        else
            b = mean;
        end 
    end
    
    if abs(f(mean))<tol
        root=mean;
        fprintf('Bisection concluded in: %d iterations \n', iter)
        %fprintf('Root is:                           %d \n', root)
        fprintf('Error:                             %d \n', abs(f(root)))
        break
    else
        fprintf('Error is:                          %d \n', abs(f(mean)))
    end
    
end


end