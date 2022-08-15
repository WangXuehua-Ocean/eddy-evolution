function dy=dy_atmos(latitude)
    R=6.3781e6; % earth's radius
    [~, dy]=gradient(latitude);
    dy=dy.*(pi./180).*R;
end