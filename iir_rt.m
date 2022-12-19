%iir 的实时实现
function  y = iir_rt(x)
   
    global iir
    
    len = length(x);
    y = zeros(len,1);
    for j = 1:len
        y(j) = sum(iir.mem.*iir.b);
        for i = length(iir.mem)-1:-1:1
            iir.mem(i+1) = iir.mem(i);
        end
        iir.mem(1) = x(j) - sum(iir.mem(2:end).*iir.a(2:end));
    end
end