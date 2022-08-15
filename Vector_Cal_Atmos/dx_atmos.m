function dx=dx_atmos(longitude, latitude)
    R=6.3781e6; % earth's radius
    [dx, ~]=gradient(longitude);
    dx=dx.*(pi./180).*R.*cos(latitude*pi./180);
end