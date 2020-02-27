%   Get MUtility
% returns exp(v(a|k)) as a matrix for all (a,k)
%pretty sure this evaluates the matrix v(a|k)
% Op.x is the parameter vector of betas
% AttsLc are the the matrices of Travel time, turn angles, leftturn and
% uturn - as a k-a matrix

%%
function Mfull = getM(x,AttsLc)   
    global incidenceFull;
    global Op;
%     fprintf("x, AttsLc = ");
%     disp(x)
%     disp(AttsLc)
%     fprintf("================\n")
    % short term value function matrix
    u = AttsLc(1).value * x(1); 
    for i = 2:Op.n
        u = u + AttsLc(i).value * x(i);
    end
    expM = u;
    expM(find(incidenceFull)) = exp(u(find(incidenceFull)));
    % This product appears to be redundant - we are multiplying 
    % element wise the non zero entries with the values that they already
    % were
    Mfull = incidenceFull .* expM;
end