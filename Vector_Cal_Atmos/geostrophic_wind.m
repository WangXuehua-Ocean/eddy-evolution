function [ug vg]=geostrophic_wind(longitude, latitude, H)
    omega=7.292*1e-5; % rotational velocity
    g=9.8; % m/s^2 % gravity
    f=2*omega*sin(latitude*pi./180); % geostrophic parameter
    [grad_x grad_y]=grad_atmos(longitude, latitude, H);
    ug=-(g./f).*grad_y;
    vg=(g./f).*grad_x;
    ug(abs(latitude)==0)=NaN;
    vg(abs(latitude)==0)=NaN;
end