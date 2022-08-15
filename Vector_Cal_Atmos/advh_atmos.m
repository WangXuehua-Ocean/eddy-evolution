function advh=advh_atmos(longitude, latitude, u, v, T)
    R=6.3781e6; % earth's radius
    dx=dx_atmos(longitude, latitude);
    dy=dy_atmos(latitude);
    [dTx dTy]=gradient(T);
    advh=-(u.*dTx./dx+v.*dTy./dy)./R;
end