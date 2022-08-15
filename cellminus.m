% cell 转成 矩阵 并计算矩阵最后一个数减去第一个数
% Created by Wang Xuehua 2022/06/29
function data=cellminus(cemat)
n=length(cemat);
data=[];
a=0;
b=0;
for i=1:n
    a=cemat{i};
    b=a(length(a))-a(1)+1;
    data=[data;b];
end
end