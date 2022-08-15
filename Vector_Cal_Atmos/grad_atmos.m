function [grad_x grad_y]=grad_atmos(longitude, latitude, H)
    dx=dx_atmos(longitude, latitude);
    dy=dy_atmos(latitude);
    [grad_x grad_y]=gradient(H);
    grad_x=grad_x./dx;
    grad_y=grad_y./dy;
end