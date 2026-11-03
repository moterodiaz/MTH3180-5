%computes the force exerted by the spring at one of its ends
%INPUTS:
%PA: the position of the first end of the spring
%PB: the position of the second end of the spring
%k: the stiffness of the spring
%l0: the natural length of the spring
%OUTPUTS:
%F: the force exerted by the spring at end B
function F = compute_spring_force(k,l0,PA,PB)
    %current length of the spring
    dP = PB(:)-PA(:);
    l = norm(dP);
    %unit vector pointing from PA to PB
    if l == 0
        error('Spring endpoints must not coincide.')
    end
    e_s = dP/l;

    %Force exerted by spring at point B
    F = -k*(l-l0)*e_s;
end
