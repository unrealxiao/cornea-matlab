%Intensity=zeros(400,100,100);
%Intensity=zeros(1000,600);  % (x,z)
%Intensity=zeros(1000,550);  % (x,z)

%delete(gcp)
%parpool(20)
tic
path='E:\\3Dimage\\POCM_System\\Human Imaging\\IR_card\\scan\\';

newpath=fullfile([path,'pic']);
mkdir(newpath);
picpath = fullfile([path, 'image_2.asc']);
pic_size = size(load(picpath));
pic_h = pic_size(1, 1);%obtain the size of the picture
pic_w = pic_size(1, 2);
parfor ii=1:500  % y
Intensity=zeros(pic_h, pic_w);
%path='E:\\3Dimage\\Intensity\\Sync_Testing\\scan1\\';
%path='C:\\Users\\vwaller\\Desktop\\3D_vwaller\\';
%path='E:\\3Dimage\\Intensity\\OCT_Pupillary_Test\\Pig_Eyeball-2021-3-30\\Eye2\\iris4-1.5V\\';
filepath=sprintf([path,'image_%d.asc'],ii);
filepathsave=sprintf([path,'pic\\frame%d.dcm'],ii);
Intensity(:,:)=load(filepath);
Intensity=uint16(rot90(Intensity));
dicomwrite(Intensity,filepathsave);
%ii
end

toc
%filepath=sprintf('C:\\Users\\Kye\\Documents\\Highspeed_OCT_Labview\\Highspeed\\3Dimage\\Loop0\\Kye_5162011_forearm_100 by 500 50 us 6 zones\\fused\\image%d%d.asc',jj,ii);
%filepath2=sprintf('C:\\Users\\Kye\\Documents\\Highspeed_OCT_Labview\\Highspeed\\3Dimage\\Loop0\\Kye_5162011_forearm_100 by 500 50 us 6 zones\\fused\\pic\\frame%d%d.gif',jj,ii);


%imwrite(Intensity,filepath2,'tif', 'WriteMode', 'append', 'Compression', 'none');
%imagesc(Intensity,[0 255]); colormap(gray);
%imwrite(Intensity,'frame0.gif');
%siz = size(Intensity);
%imshow(B,'DisplayRange',[0 255]);

