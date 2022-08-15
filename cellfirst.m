% cell 转成 矩阵 并取出第一个数
% Created by Wang Xuehua 2022/06/29
function data=cellfirst(cemat)
n=length(cemat);
data=[];
for i=1:n
    a=cemat{i};
    b=a(1);
    data=[data;b];
end
end
