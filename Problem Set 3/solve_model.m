%% Solve model 
% This code solves the wage equation for each value in the grid of the
% latent state variable 

fprintf('=============================\n')
fprintf('Solving the model...\n')
fprintf('=============================\n')

theta0 = ones(length(y_grid),1);
theta = zeros(length(y_grid),1);
wages = zeros(length(y_grid),1);


for iter=1:max_iter

    for i=1:length(y_grid)
        y = y_grid(i);
        p = p_grid(i);
        
        % Turn bargain shocks on/off
        if bargain==0
            beta = beta0;
        else
            beta = normcdf(-p)+0.5;
        end

        % Turn leisure shocks on 
        if leisure==0
            z=z0;
        else
            z=p-0.0001;
        end
        
        if p<0
            keyboard
        end
        
        %Define wage equation to be solved
        if i==1
            G = @(t) ((r+lambda+s)/mu)*t^alpha + beta*t - (1-beta)*((p-z)/c) - (lambda/mu)*(P(i,i+1)*theta0(i+1)^alpha);
        elseif i==length(y_grid)
            G = @(t) ((r+lambda+s)/mu)*t^alpha + beta*t - (1-beta)*((p-z)/c) - (lambda/mu)*(P(i,i-1)*theta0(i-1)^alpha);
        else
            G = @(t) ((r+lambda+s)/mu)*t^alpha + beta*t - (1-beta)*((p-z)/c) - (lambda/mu)*(P(i,i+1)*theta0(i+1)^alpha + P(i,i-1)*theta0(i-1)^alpha);
        end
        
        %Find theta that solves the wage equation
        t_star = fsolve(G, theta0(i), options);
        %Update theta
        theta(i) = t_star;
        wages(i) = (1-beta)*z + beta*(p + c*t_star);

    end
    
    %Compute norm
    norm = max(abs(theta(:)-theta0(:)));
    
    %Check convergence
    if norm>eps
        %Update theta
        theta0=theta;
        if mod(iter, print)==0
            fprintf('Iter = %.4f \n',iter)
            fprintf('Error = %.8f \n',norm)
        end
    else
        fprintf('Solution has been found! \n')
        fprintf('Iter = %.4f \n',iter)
        fprintf('Error = %.8f \n',norm)
        break
    end



end









