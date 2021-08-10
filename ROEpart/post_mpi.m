clear
clc
filepath='/public/lapc/Wangtao/ROE_V1.0/calculation/emission_hour/test_mpi6/';
outputpath='/public/lapc/Wangtao/ROE_V1.0/calculation/emission_hour/post_mpi/';
species={'CO','HC','NOx','PM25','PM10'};
rank=120; %100
%timestart='2021-02-05';
%timeend='2021-02-16';
%timestart='2020-07-30';
%timeend='2020-08-15';
%timestart='2020-08-16';
%timeend='2021-05-31';
timestart='2021-06-01';
timeend='2021-06-30';
timelength=datenum(timeend,'yyyy-mm-dd')-datenum(timestart,'yyyy-mm-dd')+1;

for i_time=1:timelength
    for i_spe=1:length(species)
        clear needdata1
        needdata1=[];
        for i_rank=0:rank-1
    filename=strcat(filepath,datestr(datenum(timestart,'yyyy-mm-dd')+i_time-1,'yyyy-mm-dd'),'_hour_emission_',species{i_spe},num2str(i_rank,'%.3d'),'.txt');
clear data needdata
data=importdata(filename);
for i=2:length(data.textdata)
    needdata(i-1,1)=str2num(data.textdata{i,1});
end
needdata(:,2:24)=data.data(2:end,:);
needdata1=[needdata1;needdata];
        end
outputfilename=strcat(outputpath,datestr(datenum(timestart,'yyyy-mm-dd')+i_time-1,'yyyy-mm-dd'),'_hour_emission_',species{i_spe},'.txt');
fid=fopen(outputfilename,'wt');
fprintf(fid,'%s \n','# 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23');
%dlmwrite(strcat(outputpath,'2020-10-08_hour_emission_',species{3},'.txt'),'\n\r','newline','pc','-append');
dlmwrite(outputfilename,needdata1,'delimiter',',','precision','%.6f','-append');%strcat(outputpath,'2020-10-08_hour_emission_',species{3},'.txt')
fclose(fid);
    end
end


