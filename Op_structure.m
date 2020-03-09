% Optimization structure
classdef Op_structure
   properties
      Optim_Method = OptimizeConstant.LINE_SEARCH_METHOD; % or 'BTR'
      CombinedType = 'NONE';
      Hessian_approx = OptimizeConstant.BFGS;     
      n = 4;          % size of vector
      alpha = 0.0001; % Paramters for linesearch algorithm
      gamma = 0.9 ;   % Paramters for linesearch algorithm
      maxIter = 200;
      nFev = 0;       % Number of function evaluation
      tol = 0.00001;
      x = [];
      grad = [];
      deltaGrad = [];
      step = [];
      stepLength = 0;
      k = 0;
      value = 0;
      deltaValue = 0;
      H = [];   
      Ak = [];% Hessian approximation.
      Hi;
      radius = 1;   
      % trust region radius.   
      delta = 1;
      ETA1 = 0.2;
      ETA2 = 0.75;
   end
   
   methods
       
   function [] = PrintOut(instance)
      global resultsTXT;
      fprintf('[Iteration]: %d\n', instance.k);
      fprintf('     LL = %f\n', instance.value);
      fprintf('     x = \n');
      fprintf('         %i\n', instance.x');
      fprintf('     norm of step = %f\n', norm(instance.step));
      fprintf('     radius = %f\n', instance.delta);  
      fprintf('     Norm of grad = %f\n', norm(instance.grad));
      relatice_grad = relative_gradient(instance.value, instance.x, instance.grad, 1.0);
      fprintf('     Norm of relative gradient = %f\n', relatice_grad);
      fprintf('     Number of function evaluation = %f\n', instance.nFev);
      
      % To string
      resultsTXT = [resultsTXT sprintf('[Iteration]: %d\n', instance.k)];
      resultsTXT = [resultsTXT sprintf('     LL = %f\n', instance.value)];
      resultsTXT = [resultsTXT sprintf('     x = \n')];
      resultsTXT = [resultsTXT sprintf('         %i\n', instance.x')];
      resultsTXT = [resultsTXT sprintf('     norm of step = %f\n', norm(instance.step))];
      resultsTXT = [resultsTXT sprintf('     radius = %f\n', instance.delta)];  
      resultsTXT = [resultsTXT sprintf('     Norm of grad = %f\n', norm(instance.grad))];
      resultsTXT = [resultsTXT sprintf('     Norm of relative gradient = %f\n', relatice_grad)];    
      resultsTXT = [resultsTXT sprintf('     Number of function evaluation = %d\n', instance.nFev)]; 
       
   end
   end
end