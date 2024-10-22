% Jacob Kerr, 2019

%% Processing raw data and producing results

dataDir =  '/home/labadmin/Documents/MATLAB/VarTests/Data/newdata/';
cd(dataDir);
dirlist = dir('Data*');
parfor idx = 1:length(dirlist)
    tic
    disp(strcat('Data Set ',int2str(idx)));
    pca_process_v1(strcat(dataDir,dirlist(idx).name),20);
    toc
end

%% Graphing and statistics

close all

%Variable legend
%var_scaled = compact matrix of all scaled variance readings
%sum_var = sum of all scaled variance readings
%time_scaled = compact matrix of all scaled time readings
%sum_time = sum of all scaled time readings

var_scaled = [];
sum_var = [0];
time_scaled = [];
sum_time = [0];

%Loading results and packing into matrices
for idx = 1:85
    cd(strcat('/home/labadmin/Documents/MATLAB/VarTests/Data/newdata/Data',int2str(idx)));
    load('pca_process_output.mat');
    var_scaled = [var_scaled,pca_struct.var_scale(:,1:5)];
    sum_var = sum_var + pca_struct.var_scale(:,1:5);
    time_scaled = [time_scaled;pca_struct.time_scale];
    sum_time = sum_time + pca_struct.time_scale;
end

%Fixing dimentions 
time_scaled = transpose(time_scaled);

%Calculating mean of data sets
mean_var = sum_var / 85;
mean_time = sum_time / 85;

%Plotting data
f1 = figure;
f2 = figure;
bar(mean_var,'stacked')
hold on
%plot(var_scaled)
hold off
figure(f1)
bar(mean_time)
hold on
plot(time_scaled)
hold off

%TTest for difference between PCA and gPCA
ttest_var = var_scaled(4,:) - var_scaled(3,:);
h_var = ttest(var_scaled(3,:),var_scaled(4,:));
h_var_dif = ttest(ttest_var);
ttest_time = time_scaled(4,:) - time_scaled(3,:);
h_time = ttest(time_scaled(3,:),time_scaled(4,:));
h_time_dif = ttest(ttest_time);
