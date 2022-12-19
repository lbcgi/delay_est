clc
clear
close all

global iir

frame_len = 512;
Len = frame_len*4;
lamda=0.998;
Wg=64;
h = zeros(1,Wg);
h(Wg/2) = 1;
NN=Len+Wg;
delta=0.001;

x=sqrt(10)*randn(NN,1);
m=sqrt(10);%signal power
y=filter(h,1,x);
vn=sqrt(0.1)*randn(NN,1);
dnoise=y+vn;
%%%%%%%initialization of PNLMS
en=0;
ev=0;
ex=0;
ed=0;

mu_PNLMS=1;
red=0;
%%%%%%%%%%%
%%%%%%%%%%%initialization of mu=1 NLMS
e_NLMS_1=zeros(NN,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%initialization of NPVSS
en_1=0;
%%%%%%%%%%%
%低通滤波和降采样的信号用于延迟测量，使用自适应滤波vss-nlms算法测量延迟
fs = 48000;
D = 4;
% Frames = 100;
n = 12;
Wn = 0.2;
ftype = 'low';
[b,a] = butter(n,Wn,ftype);

iir.b = b;
iir.a = a;
iir.mem = zeros(1,length(b));
% xk_lp = filter(b,a,xk);
% xk_lp_dec = xk_lp(1:D:end);
% x = xk_lp_dec;
% 
% dnoise_lp = filter(b,a,dnoise);
% dnoise_lp_dec = dnoise_lp(1:D:end);
% dnoise = dnoise_lp_dec;

h_NPVSS = zeros(1,round(Wg/D)).';
e_PNLMS=zeros(NN,1);
rex=zeros(round(Wg/D),1);
%分帧
% 遍历各帧
frms = Len/frame_len;
x_frame_lp_dec_pre = [];
dnoise_frame_lp_dec_pre = [];
for fcnt = 1:frms
   
    x_frame = x(1+(fcnt-1)*frame_len:fcnt*frame_len);
%     x_frame_lp = filter(b,a,x_frame);
    x_frame_lp = iir_rt(x_frame);
    x_frame_lp_dec = x_frame_lp(1:D:end);
    dnoise_frame = dnoise(1+(fcnt-1)*frame_len:fcnt*frame_len);
%     dnoise_frame_lp = filter(b,a,dnoise_frame);
    dnoise_frame_lp = iir_rt(dnoise_frame);
    dnoise_frame_lp_dec = dnoise_frame_lp(1:D:end);
    
    if(~isempty(x_frame_lp_dec_pre))
        x_frame_lp_dec_cur = [x_frame_lp_dec_pre(end-round(Wg/D)+1:end); x_frame_lp_dec];
        dnoise_frame_lp_dec_cur = [dnoise_frame_lp_dec_pre(end-round(Wg/D)+1:end); dnoise_frame_lp_dec];
    else
        x_frame_lp_dec_cur = x_frame_lp_dec;
        dnoise_frame_lp_dec_cur = dnoise_frame_lp_dec;
    end
    
    for k=round(Wg/D):length(x_frame_lp_dec_cur)
            X = x_frame_lp_dec_cur(k:-1:k-round(Wg/D)+1);
             %%NPVSS1
            e_NPVSS=dnoise_frame_lp_dec_cur(k)-(h_NPVSS')*X;
            en_1=lamda*en_1+(1-lamda)*e_NPVSS^2;
            en_2=sqrt(en_1);

            rex=lamda*rex+(1-lamda)*X*e_NPVSS;%互相关
            ex=lamda*ex+(1-lamda)*x(k)^2;%参考信号能量
            ev=en-(1/ex)*(rex'*rex);%什么东西

            mu_NPVSS=0.0;
            if en_2>sqrt(abs(ev))
                mu_NPVSS=1/(m+X'*X)*(1-sqrt(abs(ev))/(delta+en_2));
            end
            h_NPVSS=h_NPVSS+mu_NPVSS*X*e_NPVSS;
    
    end
    x_frame_lp_dec_pre = x_frame_lp_dec_cur;
    dnoise_frame_lp_dec_pre = dnoise_frame_lp_dec_cur;
%     if(congerged)
%     end
%     plot(h_NPVSS)
%     pause(1)
end
        

plot(h,'g')
hold on
plot(h_NPVSS,'b')
legend('原始标准延迟','低通滤波后4倍抽取后自适应滤波器的延迟估计');
