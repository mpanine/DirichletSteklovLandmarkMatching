function vf = DumDumOptimizer(funcObj, v0, numsteps)
%DUMDUMOPTIMIZER Minimizes funcObj starting at v0 following a
%gradient descent. Used for rough testing.

funcObj = @(x)autoGrad(x,1,funcObj); % Numeric differentiator


% numsteps = 5000;
v_cur = v0; % current v

stepsize = 0.0001;
stepup = 1.1;
stepdown = 0.7;

[E_cur, D_cur] = funcObj(v0);

Energies = zeros(1,numsteps);
StepSizes = zeros(1,numsteps);

for i = 1:numsteps
    
    StepSizes(i) = stepsize;
    
    v_next = v_cur - stepsize * D_cur;
    [E_next, D_next] = funcObj(v_next);

    if E_next < E_cur
       
        v_cur = v_next;
        E_cur = E_next;
        D_cur = D_next;
        stepsize = stepup*stepsize;
        
    else
           
        stepsize = stepdown*stepsize;
        
    end

    Energies(i) = E_cur;
    
    
    subplot(1,2,1)
    plot(Energies)
    title('Energy')
    subplot(1,2,2)
    plot(StepSizes)
    title('Step Size')
    pause(0.0000000001)
    
end    
    
vf = v_cur;
end

