# DNN
```matlab
%固定信噪比用此方法生成数据集
%信噪比为-19dB，采样点数为1000，协作用户数为5，检测概率与虚警概率之间的关系
clear
clc
%snrdb=-21;
snrdb=[-20,-19,-18,-17];
% SNR=(10).^(snrdb./10);
N_Simulation=11;%定义模拟的数量为11。
N_GradIterative=30;
% S=2;
N_Pulse=5;%协作用户数
fs=1500;%采样频率Hz
f=100;%信号频率Hz
Sigma=1;   %噪声方差
ns=5000;
e=0;

for iteration=1:1:600  
    e=e+1;
    e
%     SNR=(10).^(snrdb./10);
    %产生参考的噪声协方差矩阵
    N_NoiseMatrices=30;
    for i=1:N_NoiseMatrices
         for j=1:N_Pulse
             %[ signal_noise,noise ] = hunhexinghao(N_Pulse, snrdb);
             [ signal_noise,noise ] = hunhexinghao(N_Pulse, snrdb(j));
             a(:,j)=noise';
        end
        B(:,:,i)=cov(a);
    end
    %利用梯度下降算法计算噪声协方差阵的均值
    dt=0.1;     %迭代步长
    Sum_Matrices=zeros(N_Pulse,N_Pulse); %定义矩阵的算术和
    for i=1:N_NoiseMatrices     %求和
        Sum_Matrices=Sum_Matrices+B(:,:,i);
    end
    Mean_Matrices=Sum_Matrices/N_NoiseMatrices;      %矩阵的算术均值
    A=Mean_Matrices;
    for k=1:N_GradIterative    %迭代次数循环
        Sum_Log_DistanceMatrices=zeros(N_Pulse,N_Pulse); %定义对数距离矩阵的和
        for i=1:N_NoiseMatrices    %参考单元数16
            DistanceMatrices=A^(-0.5)*(B(:,:,i))*A^(-0.5); %两个矩阵的距离
            Log_DistanceMatrices=logm(DistanceMatrices);   %第i个距离矩阵的对数
            Sum_Log_DistanceMatrices=Sum_Log_DistanceMatrices+Log_DistanceMatrices;  %对数距离矩阵的和
        end
        Mean_Log_DistanceMatrices=Sum_Log_DistanceMatrices/N_NoiseMatrices;  %对数距离矩阵的均值
        A=A^0.5*expm(dt*Mean_Log_DistanceMatrices)*A^0.5;  %迭代结束：噪声的距离特征
    end
%     As=sqrt(2*Sigma^2*SNR);
    for m=1:N_Simulation
        for j=1:N_Pulse 
           % [ signal_noise,noise ] = hunhexinghao(N_Pulse, snrdb);
             [ signal_noise,noise ] = hunhexinghao(N_Pulse, snrdb(j));
            ray(:,j)=noise';
            b(:,j)= signal_noise';
        end
        sig_noi_cov=cov(b);
        noi_cov=cov(ray);

        %计算噪声和信号距离
        sig_noi_noi_Distance=(A^-0.5)*sig_noi_cov*(A^-0.5);
        sig_noi_E=eig(sig_noi_noi_Distance);
        sig_noi_D=sqrt(1/2*sum(log(sig_noi_E).^2));
        signal_noise_distance(m)=sig_noi_D;
        
        %计算噪声和噪声距离
        noi_noi_Distance=(A^-0.5)*noi_cov*(A^-0.5);
        noi_E=eig(noi_noi_Distance);
        noi_D=sqrt(1/2*sum(log(noi_E).^2));
        noise_distance(m)=noi_D;
    end%已知噪声和
%     SU(1,1:500)=signal_noise(1:500);
%     SU(1,501:1000)=noise(1:500);
%     SU(2,1:500)=signal_noise(501:1000);
%     SU(2,501:1000)=noise(501:1000);
    SU(:,2*e-1)=(signal_noise_distance)';
    SU(:,2*e)=(noise_distance)';
end

SU(11,:)=repmat([1;0],600,1)';
save SU17_21_CYS1500.mat;
```
```matlab
%信噪比在一个区间内，用此方法生成数据集
clc
clear all
tic
N_Simulation=11;%模拟的数据量为11
N_GradIterative=30;%定义梯度下降的迭代次数为30
N_Pulse=6;%次用户数
fs=1000;%采样频率Hz
f=100;%信号频率Hz
Sigma=1;   %噪声方差
%ns=500;
e=0;%初始化为0，用于之后的循环计数
for snrdb=-20:0.01:-17.01  
        e=e+1;
        e
%产生参考的噪声协方差矩阵
    N_NoiseMatrices=30;%定义噪声矩阵的数量为30
    for i=1:N_NoiseMatrices%对噪声矩阵进行遍历
         for j=1:N_Pulse
              [ signal_noise,noise ] = hunhexinghao(N_Pulse, snrdb);
            a(:,j)=noise';%将噪声矩阵存储到a矩阵的第j列
        end
        B(:,:,i)=cov(a);%计算噪声协方差并将结果存到矩阵B的第i页
    end
%利用梯度下降算法计算噪声协方差阵的均值
    dt=0.1;     %迭代步长
    Sum_Matrices=zeros(N_Pulse,N_Pulse); %定义矩阵的算术和
    for i=1:N_NoiseMatrices     %求和
        Sum_Matrices=Sum_Matrices+B(:,:,i);
    end
    Mean_Matrices=Sum_Matrices/N_NoiseMatrices;      %矩阵的算术均值
    A=Mean_Matrices;
    for k=1:N_GradIterative    %迭代次数循环
        Sum_Log_DistanceMatrices=zeros(N_Pulse,N_Pulse); %定义对数距离矩阵的和
        for i=1:N_NoiseMatrices    %参考单元数16
            DistanceMatrices=A^(-0.5)*(B(:,:,i))*A^(-0.5); %两个矩阵的距离
            Log_DistanceMatrices=logm(DistanceMatrices);   %第i个距离矩阵的对数
            Sum_Log_DistanceMatrices=Sum_Log_DistanceMatrices+Log_DistanceMatrices;  %对数距离矩阵的和
        end
         Mean_Log_DistanceMatrices=Sum_Log_DistanceMatrices/N_NoiseMatrices;  %对数距离矩阵的均值
         A=A^0.5*expm(dt*Mean_Log_DistanceMatrices)*A^0.5;  %迭代结束：噪声的距离特征
    end
  % As=sqrt(2*Sigma^2*SNR);
    
    for m=1:N_Simulation
        for j=1:N_Pulse 
             [ signal_noise,noise ] = hunhexinghao(N_Pulse, snrdb);
             ray(:,j)=noise';
             b(:,j)= signal_noise';
       
       end
        sig_noi_cov=cov(b);
        noi_cov=cov(ray);

        %计算噪声和信号距离
        sig_noi_noi_Distance=(A^-0.5)*sig_noi_cov*(A^-0.5);
        sig_noi_E=eig(sig_noi_noi_Distance);
        sig_noi_D=sqrt(sum(log(sig_noi_E).^2));
        signal_noise_distance(m)=sig_noi_D;
        
        %计算噪声和噪声距离
        noi_noi_Distance=(A^-0.5)*noi_cov*(A^-0.5);
        noi_E=eig(noi_noi_Distance);
        noi_D=sqrt(sum(log(noi_E).^2));
        noise_distance(m)=noi_D;
    end%已知噪声和
% SU(1,1:500)=signal_noise(1:500);
% SU(1,501:1000)=noise(1:500);
% SU(2,1:500)=signal_noise(501:1000);
% SU(2,501:1000)=noise(501:1000);
SU(:,2*e-1)=(signal_noise_distance)';
SU(:,2*e)=(noise_distance)';
end
SU(11,:)=repmat([1;0],600,1)';
toc
disp(['运行时间',num2str(toc)]);
save SU_20_11;
```
```matlab
clc
clear
%%数据导入
%DATA=importdata('matlab19-13.mat');%matlab为数据格式假设样本个数200，特征个数为30，200*31，第31列表示分类标签，这里分为0-3共4类
load SU17_21_CYS1500.mat
DATA = SU;
data=DATA';
k=rand(1,1000);         %生成样本个数的随机数，以200为例
[m,n]=sort(k);         %若k为二维矩阵，则sort（k）表示对每列进行升序排序,m为移动后矩阵，n为移动顺序。

%上两步的操作是为了生成随机数，之后按比例随机挑选样本；
%输入输出数据
input=data(:,1:10);
output=data(:,11);           %第101列为标签列
%随机抽取1000个样本为训练样本，200个样本为预测样本
input_train=input(n(1:800),:)';  %取打乱顺序后的前1000个数字，行为样本列为特征点
label_train=output(n(1:800),:)';%类似于标签
input_test=input(n(801:1000),:)';
label_test=output(n(801:1000),:)';

%输入数据归一化
 [inputn,inputps]=mapminmax(input_train);

%% BP网络训练
% %初始化网络结构
%    net=newff(inputn,label_train,10);
%      net=newff(minmax(inputn),[10,7,1],{'tansig','tansig','purelin'},'trainlm');
%     net=newff(minmax(inputn),[10,9,7,1],{'tansig','tansig','tansig','purelin'},'trainlm');
   net=newff(minmax(inputn),[10,9,7,5,1],{'tansig','tansig','tansig','tansig','purelin'},'trainlm');%网络结构设置
%  net=newff(minmax(inputn),[10,9,7,7,5,1],{'tansig','tansig','tansig','tansig','tansig','purelin'},'trainlm');
net.trainParam.epochs=1000;%网络迭代次数epochs为1000次
net.trainParam.lr=0.1;%学习速率lr为0. 1
net.trainParam.goal=0.0001;%期望误差goal为0.00000004
% net.trainParam.max_fail=100; 
% net.divideFcn = '';
%% 网络训练
net=train(net,inputn,label_train);

%% BP网络预测
%预测数据归一化
tic
inputn_test=mapminmax('apply',input_test,inputps);
 
%网络预测输出
BPoutput=sim(net,inputn_test);
 aa=mean(BPoutput);

%% 结果分析
%根据网络输出找出数据属于哪类
BPoutput(find(BPoutput<0.5))=0;
BPoutput(find(BPoutput>=0.5))=1;
toc
%% 结果分析
% figure(1)
% plot(BPoutput,'og')
% ylim([-0.5 1.5])
% legend('目标标签')
% title('BP网络预测类别分布','fontsize',12)
% %画出预测种类和实际种类的分类图
% figure(2)
% plot(BPoutput,'og')
% hold on
% plot(label_test,'r*');
% legend('目标标签','实际输出')
% title('BP网络预测分类与实际类别比对','fontsize',12)
% ylabel('类别标签','fontsize',12)
% xlabel('样本数目','fontsize',12)
% ylim([-0.5 1.5])




%预测正确率
rightnumber=0;
for i=1:size(label_test,2)%size(output_test,2)指output_test的列
    if BPoutput(i)==label_test(i)
        rightnumber=rightnumber+1;
    end
end
rightratio=rightnumber/size(label_test,2)*100;

sprintf('测试准确率=%0.2f',rightratio)

%画ROC曲线
 figure(3);
 BPoutput1= mapminmax(BPoutput,0,1);
 plot_roc(BPoutput1,label_test);
%  plotroc(label_test,BPoutput);
grid on
```
```matlab
function  auc = plot_roc( predict, ground_truth )  
% INPUTS  
%  predict       - 分类器对测试集的分类结果  
%  ground_truth - 测试集的正确标签,这里只考虑二分类，即0和1  
% OUTPUTS  
%  auc            - 返回ROC曲线的曲线下的面积  
  
%初始点为（1.0, 1.0）  
x = 1.0;  
y = 1.0;  
%计算出ground_truth中正样本的数目pos_num和负样本的数目neg_num  
pos_num = sum(ground_truth==1);  
neg_num = sum(ground_truth==0);  
%根据该数目可以计算出沿x轴或者y轴的步长  
x_step = 1.0/neg_num;  
y_step = 1.0/pos_num;  
%首先对predict中的分类器输出值按照从小到大排列  
[predict,index] = sort(predict);  
ground_truth = ground_truth(index);  
%对predict中的每个样本分别判断他们是FP或者是TP  
%遍历ground_truth的元素，  
%若ground_truth[i]=1,则TP减少了1，往y轴方向下降y_step  
%若ground_truth[i]=0,则FP减少了1，往x轴方向下降x_step  
for i=1:length(ground_truth)  
    if ground_truth(i) == 1  
        y = y - y_step;  
    else  
        x = x - x_step;  
    end  
    X(i)=x;  
    Y(i)=y;  
end  
%画出图像  
% c = polyfit(X, Y, 6);  %进行拟合，c为2次拟合后的系数
% d = polyval(c, X, 1);  %拟合后，每一个横坐标对应的值即为d
% plot(X,d,'*-','LineWidth',1.5,'MarkerSize',1.5);  

% plot(X,Y,'*-','LineWidth',1.5,'MarkerSize',1.5);  
%画出图像  
X1=zeros(1,21);
Y1=zeros(1,21);
X1=X(1);
Y1=Y(1);
for n = 1:20
    X1(n+1)=X(:,10*n);
    Y1(n+1)=Y(:,10*n);
end
plot(X1,Y1,'*-','linewidth',1.5);  
xlabel('Pf');  
ylabel('Pd');  
title('ROC曲线图');
axis([0 1 0 1])

% saveas(gcf, 'IGDNN_SNR21', 'fig')

%计算小矩形的面积,返回auc  
auc = -trapz(X,Y);

end 
```
```matlab
function [ signal_noise,noise ] = hunhexinghao(N_Pulse, snrdb )
Dn=1;       %噪声方差N_Pulse
% fs=1000;%采样频率
% n=0:N-1;
% t=n/fs;  
t= 0:0.05:15.9*pi;
SNR_in=10^(snrdb/10);  %将db转化为十进制
Ap=sqrt(2*Dn*SNR_in);
% signal=cos(2*pi*10*t);
signal=cos(t)+cos(4*t+0.2*t.^2);
%signal=cos(2*pi*t)+2*cos((2*pi-4)*t);
noise=sqrt(Dn)*randn(size(t));
signal_noise=Ap*signal+noise;
end
```
