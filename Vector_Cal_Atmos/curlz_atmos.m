function curlz=curlz_atmos(longitude, latitude, u, v)
    R=6.3781e6; % earth's radius
    dx=dx_atmos(longitude, latitude);
    dy=dy_atmos(latitude);
    [~, du]=gradient(u);
    [dv, ~]=gradient(v);
    curlz=dv./dx-du./dy+u.*tan(latitude.*pi./180)./R;
    curlz(abs(latitude)==90)=NaN;
end