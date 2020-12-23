function param_selection

% This script is for selecting the values of A and B based on 

% coarse_sd3(); % 1 coarse the dataset so that we sample at 0.01 s

% calculate_distmat(); % 3 this takes time! go for a coffee break!

% calculate_nn_dist()

plot_nn_dist()

function coarse_sd3()
% subsample the trajectories at 0.01 s intervals

files=dir('./PanicPackage/dump/sd_A*.dat3');
for ff=1:size(files,1)
    data=dlmread(['./PanicPackage/dump/', files(ff).name]);
    fprintf('processing %s ...\n', files(ff).name);
    % sort time into 1/100 dt increments
    dt1=1/100;
    simtime=unique(data(:,1));

    keeptime=simtime(1);
    for ii=2:numel(simtime)
        simtime1=simtime(simtime-keeptime(end)-dt1>0);
        [val, idx]=min(simtime1-keeptime(end)-dt1);
        keeptime=[keeptime; simtime1(idx)];
    end

    data2=[];
    for ii=1:numel(keeptime)
        data2=[data2; data(data(:,1)==keeptime(ii),:)];
    end
    
    % remove any 0,0 points
    data2=data2(data2(:,3)~=0,:);
    
    % write to csv
    csvwrite(['./PanicPackage/dump/', files(ff).name(1:end-4), 'csv'], data2);
end


function calculate_distmat()

inf1=100000;
files=dir('./PanicPackage/dump/sd_A*.csv');

for ff=1:size(files,1)

    fprintf('processing %s ...\n', files(ff).name);
    
    traj=csvread(['./PanicPackage/dump/', files(ff).name]);
    ids=unique(traj(:,2));
    
    k=sort(unique(traj(:,1)));
    dt=k(2)-k(1);
    
    D=zeros(numel(ids), numel(ids), numel(k));
    
    for kk=1:numel(k)
        for ii=1:numel(ids)
            pos_ii=traj(traj(:,1)==k(kk) & traj(:,2)==ids(ii),3:4);
            for jj=1:numel(ids)
                pos_jj=traj(traj(:,1)==k(kk) & traj(:,2)==ids(jj),3:4);
                if ~isempty(pos_ii) && ~isempty(pos_jj)
                    if pos_ii(1)<10 && pos_jj(1)<10
                        D(ii,jj,kk)=norm(pos_ii-pos_jj);
                    else
                        D(ii,jj,kk)=inf1;
                    end
                else
                    D(ii,jj,kk)=inf1;
                end
                if jj==ii
                    D(ii,jj,kk)=inf1;
                end
            end
        end
    end
    
    
    
    save(sprintf('./PanicPackage/dump/distmat_%s.mat', ...
        files(ff).name(1:end-4)), 'D', 'ids', 'dt', 'k');
end



function calculate_nn_dist()


files=dir('./PanicPackage/dump/distmat_sd_A*.mat');

iter=10;

for ff=1:size(files,1)

    fprintf('processing %s ...\n', files(ff).name);
    load(['./PanicPackage/dump/', files(ff).name], 'D', 'ids', 'k', 'dt');
        
    nn_dist=zeros(numel(ids),numel(k));
    for kk=1:numel(k)
        D1=D(:,:,kk);
        nn_dist(:,kk)=min(D1,[],2);
        nn_dist(nn_dist(:,kk)>10*sqrt(2),kk)=nan;
    end
    
    csvwrite(sprintf('./PanicPackage/dump/nn_dist_%s.csv', files(ff).name(9:end-4)),...
        nn_dist);
end



function plot_nn_dist()

simfiles=dir(['./PanicPackage/dump/nn_dist_sd_A*.*']);

for ss=1:size(simfiles,1)
    nndist=csvread(['./PanicPackage/dump/', simfiles(ss).name]);
    mu(ss)=mean(nanmean(nndist,2));
    st(ss)=std(nanmean(nndist,2));
    A(ss)=str2double(simfiles(ss).name(14:16));
    B(ss)=str2double(simfiles(ss).name(20:23));
end

AA=reshape(A,9,10);
BB=reshape(B,9,10);
MU=reshape(mu,9,10);
ST=reshape(st,9,10);
figure(1); gcf; clf;
subplot(1,2,1);
surf(AA,BB,MU);
colorbar
view(2);
axis tight
title('mean NN distance');
xlabel('A'); ylabel('B');
subplot(1,2,2);
surf(AA,BB,ST);
colorbar
view(2);
axis tight
title('std NN distance');
xlabel('A'); ylabel('B');


keyboard