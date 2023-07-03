close all;clear all

Path='E:\3Dimage\POCM_System\Human Imaging\IR_card\scan\pic\cross_unflat\';
%Path='C:\\3Dimage\\Sync_Calibration\\16_microsecond_1\\pic_1 frame missed\\';
%Path='C:\Users\Patrice TANKAM\Documents\Rochester\Program\1400-June_11_after_correction\';
%Path='C:\Users\Patrice TANKAM\Documents\Rochester\Program\MEMS-1400_frames-corrected\';

 FF1=dicomread([Path,'corl109.DCM']);
 FF2=dicomread([Path,'corl110.DCM']);
%  FF1=rot90(F1);
%  FF2=rot90(F2);

size = size(FF1);
 
FF1 = FF1(:, 1:(size(1, 2)));
FF2 = FF2(:, 1:(size(1, 2)));
 
%  F1=imread([Path,'frame951.gif']);
%  F2=imread([Path,'frame952.gif']);
%  FF1=F1;%flipud(rot90(F1));
%  FF2=F2;%flipud(rot90(F2));


figure;
subplot(2,1,1);imagesc(FF1);colormap(gray);title('Forward Frame');
subplot(2,1,2);imagesc(FF2);colormap(gray);title('Backward Frame');

l1=sum(FF1);
l2=sum(FF2);
len=length(l1);

figure; plot(l1,'b'); hold on; plot(l2, 'r');grid

%corel=xcorr(l1,l2);

%figure; plot (corel);grid
%[M, shift]=max(corel)

Mcorel=xcorr2(FF1,FF2);
figure; imagesc(Mcorel);colormap(gray);grid;title('Images cross correlation');
[m, pos]=max(max(Mcorel))
shift=pos-len