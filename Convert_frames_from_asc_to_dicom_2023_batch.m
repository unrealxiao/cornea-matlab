
clear all;
tic

path=uigetdir('E:\\3Dimage\\POCM_System\Human_Imaging\OK1_V1_2023_2_7\\',' Choose the MAIN FOLDER');
nb_scan =34;

for scan =22:nb_scan
folder=['\scan',num2str(scan),'\'];
pathfile=[path,folder];
newpath=fullfile([pathfile,'pic']);%create a new folder to store the dicom files
mkdir(newpath);  
file0=[pathfile,'image_1.asc'];
Int0(:,:)=load([file0]);% get the dimension of the images
fprintf(['Processing...scan%d ...'],scan);
parfor ii=1:500 %987   % y number of frames per scan
Intensity=zeros(size(Int0)); %For OCMFM systemnewpath=fullfile([path,'pic1']);
filepath=strcat([pathfile,'\image_',num2str(ii),'.asc']);
filepathsave=strcat([newpath,'\frame',num2str(ii),'.dcm']);
Intensity(:,:)=load(filepath);
Intensity=uint16(rot90(Intensity)); 
dicomwrite(Intensity,filepathsave);
% ii
end
fprintf('Completed!!! ');
toc

end